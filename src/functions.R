require(devtools)
load_all("~/git/rmaize")
require(ape)
require(ggtree)
require(ggforce)
require(Rtsne)
dirp = '~/projects/rnaseq'
dird = file.path(dirp, 'data')
dirc = '/scratch.global/zhoux379/rnaseq'
t_cfg = read_gspread_master(lib='rnaseq')
#f_yml = file.path(dird, '10.cfg.yaml')
#Sys.setenv("R_CONFIG_FILE" = f_yml)

read_rnaseq <- function(yid) {
    #{{{
    res = rnaseq_cpm(yid)
    th = res$th; tm = res$tm
    th = th %>% replace_na(list(Tissue='',Genotype='B73',Treatment='')) %>%
        mutate(Tissue=as.character(Tissue)) %>%
        mutate(Genotype=as.character(Genotype)) %>%
        mutate(Treatment=as.character(Treatment))
    yids_dev = c('rn10a','rn11a','rn13b','rn14b','rn14c','rn14e',"rn16b","rn16c","rn18g")
    if(yid == 'rn12a') {
        th = th %>% filter(Treatment == 'WT')
    } else if(yid == 'rn17c') {
        th = th %>% filter(Treatment == 'con')
    } else if(yid %in% c(yids_dev,'rn19c')) {
        if(yid == 'rn13b') th = th %>% filter(!str_detect(Treatment, "^ET"))
        if(yid == 'rn18g') th = th %>% filter(Genotype == 'B73')
        th = th %>% mutate(Tissue=str_c(Tissue,Treatment, sep="_")) %>%
            mutate(Treatment=yid)
    }
    th = th %>% mutate(study = yid) %>%
        mutate(SampleID = str_c(study, SampleID, sep='_')) %>%
        replace_na(list(Treatment='')) %>%
        select(SampleID, Tissue, Genotype, Treatment, Replicate, study)
    tm = tm %>% mutate(SampleID = str_c(yid, SampleID, sep='_')) %>%
        filter(SampleID %in% th$SampleID)
    list(th=th, tm=tm)
    #}}}
}
read_multiqc_trimmomatic <- function(fi, paired = T) {
    #{{{
    ti = read_tsv(fi)
    types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
    if (paired == F) {
        nd = ti %>% mutate(nd = input_reads - surviving - dropped) %>%
            group_by(1) %>% summarise(nd = sum(nd)) %>% pull(nd)
        stopifnot(nd == 0)
        to = ti %>% mutate(SampleID = Sample, total = input_reads, 
                           surviving_f=0, surviving_r=0)
    } else if(paired == T | paired == 'both') {
        ti2 = ti %>% 
            separate(Sample, c("SampleID",'suf'), sep="_", fill='right', extra='merge') %>% select(-suf) %>%
            mutate(surviving_f = forward_only_surviving,
                   surviving_r = reverse_only_surviving)
        if(paired == 'both') 
            ti2 = ti2 %>% 
                replace_na(list('input_reads'=0, 'input_read_pairs'=0,
                                'surviving_f'=0, 'surviving_r'=0)) %>%
                mutate(input_read_pairs = 
                    ifelse(input_read_pairs == 0, input_reads, input_read_pairs))
        nd = ti2 %>% mutate(nd = input_read_pairs - surviving - 
                            surviving_f - surviving_r - dropped) %>%
            group_by(1) %>% summarise(nd = sum(nd)) %>% pull(nd)
        stopifnot(nd == 0)
        to = ti2 %>% mutate(total = input_read_pairs)
    } else {
        stop(sprintf("unsupported option: %s", paired))
    }
    to %>%
        select(SampleID, total, surviving, surviving_f, surviving_r, dropped)
    #}}}
}
read_multiqc_star <- function(fi, paired = T) {
    #{{{
    ti = read_tsv(fi)
    if(paired == F) {
        ti2 = ti %>% mutate(SampleID = Sample) %>% select(-Sample)
    } else {
        ti2 = ti %>% separate(Sample, c("SampleID",'suf'), sep="_", fill='right', extra='merge') %>% 
            select(-suf) 
    }
    ti2 = ti2 %>%
        transmute(SampleID = SampleID, total = total_reads,
                  uniquely_mapped = uniquely_mapped,
                  multimapped = multimapped + multimapped_toomany,
                  unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
                  nd = total - uniquely_mapped - multimapped - unmapped) 
    stopifnot(sum(ti2$nd) < 1000)
    ti2 = ti2 %>% group_by(SampleID) %>%
        summarise(uniquely_mapped = sum(uniquely_mapped),
                  multimapped = sum(multimapped),
                  unmapped = sum(unmapped))
    types = c("uniquely_mapped", "multimapped", "unmapped")
    to = ti2 %>% select(SampleID, uniquely_mapped, multimapped, unmapped)
    to
    #}}}
}
read_multiqc_hisat2 <- function(fi, paired = T) {
    #{{{
    ti = read_tsv(fi)
    if(paired == F) {
        ti2 = ti %>%
            transmute(SampleID = Sample,
                      total = unpaired_total,
                      uniquely_mapped = unpaired_aligned_one,
                      multimapped = unpaired_aligned_multi,
                      unmapped = unpaired_aligned_none)
    } else {
        ti2 = ti %>%
            transmute(SampleID = Sample,
                      #total = paired_total + unpaired_total,
                      #uniquely_mapped = paired_aligned_one+paired_aligned_discord_one+unpaired_aligned_one,
                      #multimapped = paired_aligned_multi+unpaired_aligned_multi,
                      #unmapped = paired_aligned_none+unpaired_aligned_none)
                      total = paired_total,
                      uniquely_mapped = paired_aligned_one+paired_aligned_discord_one,
                      multimapped = paired_aligned_multi,
                      unmapped = paired_aligned_none)
    }
    ti2 = ti2 %>% mutate(nd = total - uniquely_mapped - multimapped - unmapped)
    cat(sum(ti2$nd),"\n")
    stopifnot(sum(ti2$nd) < 1000)
    to = ti2 %>% group_by(SampleID) %>%
        summarise(uniquely_mapped = sum(uniquely_mapped),
                  multimapped = sum(multimapped),
                  unmapped = sum(unmapped))
    types = c("uniquely_mapped", "multimapped", "unmapped")
    to
    #}}}
}
read_multiqc_featurecounts <- function(fi) {
    #{{{
    ti = read_tsv(fi)
    ti2 = ti %>% mutate(SampleID = Sample) %>%
        select(SampleID, Total, Assigned, Unassigned_Unmapped,
               Unassigned_MultiMapping,
               Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
        mutate(nd = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
    stopifnot(sum(as.numeric(ti2$nd)) == 0)
    #
    types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
              "Unassigned_NoFeatures", "Unassigned_Ambiguity")
    to = ti2 %>% select(SampleID, Assigned, Unassigned_MultiMapping,
                        Unassigned_NoFeatures, Unassigned_Ambiguity,
                        Unassigned_Unmapped)
    to
    #}}}
}
plot_pca0 <- function(tp, fo, opt = 'col=tis', labsize = 2.5, wd = 8, ht = 8) {
    #{{{
    if(opt == 'col=tis,sha=rep') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, color = Tissue, shape = Replicate)) +
            geom_point(size = 1.5) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'lab=tis,sha=rep') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, label = Tissue, shape = Replicate)) +
            geom_point(size = 1.5) +
            geom_text_repel(size = labsize) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'lab=tre,sha=rep') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, label = Treatment, shape = Replicate)) +
            geom_point(size = 1.5) +
            geom_text_repel(size = labsize) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'lab=tre') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, label = Treatment)) +
            geom_point(size = 1.5) +
            geom_text_repel(size = labsize) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'lab=gen,sha=tis') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, label = Genotype, shape = Tissue)) +
            geom_point(size = 1.5) +
            geom_text_repel(size = labsize) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'lab=tis,sha=gen') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, label = Tissue, shape = Genotype)) +
            geom_point(size = 1.5) +
            geom_text_repel(size = labsize) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'lab=tis,col=sid,sha=sid') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, label = Tissue, color = sid, shape = sid)) +
            geom_point(size = 1.5) +
            geom_text_repel(size = labsize) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'col=tis,sha=tre') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, color = Tissue, shape = Treatment)) +
            geom_point(size = 1.5) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', fill = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'col=tis') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, color = Tissue)) +
            geom_point(size = 1.5) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'col=tis,sha=tis') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, color = Tissue, shape = Tissue)) +
            geom_point(size = 1.5) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'col=tre,sha=rep') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, color = Treatment, shape = Replicate)) +
            geom_point(size = 1.5) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else if(opt == 'sha=rep') {
        #{{{
        p1 = ggplot(tp, aes(x = PC1, y = PC2, shape = Replicate)) +
            geom_point(size = 1.5) +
            scale_x_continuous(name = xlab) + scale_y_continuous(name = ylab) +
            scale_color_d3() +
            scale_shape_manual(values = shapes) +
            guides(direction = 'vertical', color = guide_legend(ncol = 1)) +
            guides(shape = guide_legend(ncol = 1, byrow = T)) +
            otheme(legend.pos = 'top.left', xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T)
        #}}}
    } else {
        stop(sprintf("unknown opt: %s", opt))
    }
    ggsave(p1, filename = fo, width = wd, height = ht)
    #}}}
}

