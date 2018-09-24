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
options(stringsAsFactors = FALSE)
dirr = '~/git/luffy/r'
source(file.path(dirr, 'plot.R'))

dirp = '~/projects/maize.expression/data'
dird = '~/projects/maize.expression/data'
dirc = '/scratch.global/zhoux379/maize.expression'
f_cfg = file.path(dirp, '01.cfg.tsv')
t_cfg = read_tsv(f_cfg)
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
        ti2 = ti %>% rename(SampleID = Sample)
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

merge_me_datasets <- function(sids, t_cfg, dird, group = 'Tissue') {
    #{{{
    th = tibble(); t_rc = tibble()
    for (sid1 in sids) {
        study1 = t_cfg %>% filter(sid == !!sid1) %>% pull(author)
        diri = file.path(dird, sid1, 'data')
        fh1 = file.path(diri, '01.reads.tsv')
        fh2 = file.path(diri, '02.reads.corrected.tsv')
        fh = ifelse(file.exists(fh2), fh2, fh1) 
        th1 = read_tsv(fh) 
        if(group == 'Tissue') {
            th1 = th1 %>% mutate(Tissue = sprintf("%s_%s", study1, Tissue))
        } else if(group == 'Genotype') {
            th1 = th1 %>% mutate(Genotype = sprintf("%s_%s", study1, Genotype))
        } else if(group == 'Treatment') {
            th1 = th1 %>% mutate(Treatment = sprintf("%s_%s", study1, Treatment))
        } else {
            stop("unsupported group")
        }
        th = rbind(th, th1)
        fi = file.path(diri, 'raw/featurecounts.tsv')
        t_rc1 = read_tsv(fi) %>% select(one_of(c('gid', th1$SampleID)))
        stopifnot(ncol(t_rc1) == nrow(th1) + 1)
        if(nrow(t_rc) == 0)
            t_rc = t_rc1
        else {
            t_rc = t_rc %>% inner_join(t_rc1, by = 'gid')
        }
    }
    list(th = th, t_rc = t_rc)
    #}}}
}

create_cache_dir <- function(sid, dird, dirc) {
    #{{{
    dirw = file.path(dird, sid)
    cmd = sprintf("mkdir -p %s", dirw)
    system(cmd)
    dirc1 = file.path(dirc, sid)
    cmd = sprintf("mkdir -p %s", dirc1)
    system(cmd)
    if(file.exists(file.path(dirw, 'cache'))) system(sprintf("rm %s/cache", dirw))
    cmd = sprintf("ln -sf %s/ %s/cache", dirc1, dirw)
    system(cmd)
    fread = sprintf('%s/05_read_list/%s.tsv', dird, sid)
    stopifnot(file.exists(fread))
    cmd = sprintf("ln -sf %s %s/01.reads.tsv", fread, dirc1)
    system(cmd)
    cmd = sprintf("mkdir -p %s/raw", dirw)
    system(cmd)
    cmd = sprintf("ln -sf %s/raw/ %s/data", dirw, dirc1)
    system(cmd)
    #}}}
}


