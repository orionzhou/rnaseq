#{{{ load required libraries, define common variables
require(grid)
require(tidyverse)
require(gtable)
require(ggtree)
require(RColorBrewer)
require(viridis)
require(cluster)
require(Hmisc)
require(ggsignif)
require(cowplot)
require(GGally)
require(ggridges)
require(ggpubr)
require(ggsci)
require(ggrepel)
require(scales)
require(pheatmap)
require(yaml)
options(stringsAsFactors = FALSE)
dirr = '~/git/luffy/r'
source(file.path(dirr, 'plot.R'))
#
dirp = '~/projects/maize.expression/data'
dird = '~/projects/maize.expression/data'
dirc = '/scratch.global/zhoux379/maize.expression'
f_cfg = file.path(dird, '01.cfg.tsv')
t_cfg = read_tsv(f_cfg)
f_yml = file.path(dird, '11.cfg.yaml')
Sys.setenv("R_CONFIG_FILE" = f_yml)
#}}}

readcount_norm <- function(t_rc, t_gs) {
#{{{ normalize
    tm = t_rc
    tw = tm %>% 
        select(SampleID, gid, ReadCount) %>%
        spread(SampleID, ReadCount)
    gids = tw$gid
    twd = data.frame(tw[,-1])
    rownames(twd) = tw$gid
    # nRC with DESeq2
    require(DESeq2)
    th = tm %>% distinct(SampleID) %>% arrange(SampleID)
    thd = data.frame(th)
    rownames(thd) = th$SampleID
    stopifnot(identical(thd$SampleID, colnames(twd)))
    dds = DESeqDataSetFromMatrix(countData=twd, colData=thd, design = ~SampleID)
    dds = estimateSizeFactors(dds)
    sf = sizeFactors(dds)
    t_sf = tibble(SampleID = names(sf), sizeFactor = as.numeric(sf))
    t_nrc = counts(dds, normalized = T) %>% as_tibble() %>%
        mutate(gid = names(dds)) %>% gather(SampleID, nRC, -gid)
    # rCPM and CPM with edgeR
    require(edgeR)
    y = DGEList(counts = twd)
    y = calcNormFactors(y, method = 'TMM')
    t_nf = y$samples %>% as_tibble() %>% 
        mutate(SampleID = rownames(y$samples)) %>%
        transmute(SampleID = SampleID, libSize = lib.size, normFactor = norm.factors)
    t_cpm = cpm(y) %>% as_tibble() %>% mutate(gid = rownames(cpm(y))) %>%
        select(gid, everything()) %>%
        gather(SampleID, CPM, -gid)
    t_rcpm = cpm(y, normalized = F) %>% as_tibble() %>% 
        mutate(gid = rownames(cpm(y))) %>%
        gather(SampleID, rCPM, -gid)
    # rFPKM & FPKM
    t_cpm = t_cpm %>% left_join(t_gs, by = 'gid') %>% 
        mutate(FPKM = CPM / (size / 1000)) %>%
        select(-size)
    t_rcpm = t_rcpm %>% left_join(t_gs, by = 'gid') %>% 
        mutate(rFPKM = rCPM / (size / 1000)) %>%
        select(-size)
#
    tl = th %>% inner_join(t_sf, by = 'SampleID') %>%
        inner_join(t_nf, by = 'SampleID')
    tm = tm %>% 
        left_join(t_nrc, by = c('SampleID','gid')) %>%
        left_join(t_rcpm, by = c('SampleID','gid')) %>%
        left_join(t_cpm, by = c('SampleID','gid')) 
    stopifnot(nrow(tm) == length(gids) * nrow(th))
    list(tl = tl, tm = tm)
#}}}
}

get_read_list <- function(dird, sid) {
    #{{{
    fh1 = sprintf("%s/05_read_list/%s.tsv", dird, sid)
    fh2 = sprintf("%s/05_read_list/%s.c.tsv", dird, sid)
    fh = ifelse(file.exists(fh2), fh2, fh1) 
    read_tsv(fh)
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
            separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
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
        ti2 = ti %>% separate(Sample, c("SampleID", 'suf'), sep = "_") %>% 
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

read_multiqc <- function(diri, th) {
    #{{{
    paired = unique(th$paired)
    if(length(paired) == 2) paired = 'both'
    fi = file.path(diri, "multiqc_trimmomatic.txt")
    tt1 = read_multiqc_trimmomatic(fi, paired = paired)
    fi = file.path(diri, 'multiqc_star.txt')
    tt2 = read_multiqc_star(fi, paired = paired)
    fi = file.path(diri, 'multiqc_featureCounts.txt')
    tt3 = read_multiqc_featurecounts(fi)
    tt = th %>% select(-paired) %>% 
        left_join(tt1, by = 'SampleID') %>%
        left_join(tt2, by = 'SampleID') %>%
        left_join(tt3, by = 'SampleID')
    tt %>% group_by(Tissue, Genotype, Treatment) %>%
        summarise(total = sum(total), Assigned = sum(Assigned)) %>%
        ungroup() %>% group_by(1) %>%
        summarise(total_median = median(total/1000000),
                  total_mean = mean(total/1000000),
                  assigned_median = median(Assigned/1000000),
                  assigned_mean = mean(Assigned/1000000)) %>% print(n=1)
    tt
    #}}}
}

plot_hclust_tree <- function(tree, tp, fo, labsize = 3, x.expand = .2, x.off = .05, wd = 6, ht = 8) {
    #{{{
    cols1 = c('gray80','black','red','seagreen3', pal_d3()(5))
    p1 = ggtree(tree) +
        #geom_tiplab(size = labsize, color = 'black') +
        scale_x_continuous(expand = expand_scale(mult=c(.02,x.expand))) +
        scale_y_discrete(expand = c(.01,0)) +
        theme_tree2()
    p1 = p1 %<+% tp +
        geom_tiplab(aes(label = lab), size = labsize, offset = x.off, family = 'mono')
        #geom_text(aes(label = SampleID), size = 2, nudge_x = .001, hjust = 0) +
        #geom_text(aes(label = Genotype), size = 2, nudge_x= .015, hjust = 0) +
        #geom_text(aes(label = Rep), size = 2, nudge_x = .022, hjust = 0)
    ggsave(p1, filename = fo, width = wd, height = ht)
    #}}}
}

plot_pca <- function(tp, fo, opt = 'col=tis', labsize = 2.5, wd = 8, ht = 8) {
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


