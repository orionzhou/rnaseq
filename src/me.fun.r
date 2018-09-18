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
    if (paired == F) {
        types = c("surviving", "dropped")
        nd = ti %>% mutate(nd = input_reads - surviving - dropped) %>%
            group_by(1) %>% summarise(nd = sum(nd)) %>% pull(nd)
        stopifnot(nd == 0)
        to = ti %>% mutate(SampleID = Sample) %>%
            select(SampleID, surviving, dropped) %>%
            gather(type, nseq, -SampleID)
    } else if(paired == T) {
        types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
        ti2 = ti %>% 
            separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
            select(SampleID, input_read_pairs, surviving, forward_only_surviving,
                   reverse_only_surviving, dropped) 
        nd = ti2 %>% mutate(nd = input_read_pairs - surviving - 
                            forward_only_surviving - reverse_only_surviving -
                            dropped) %>%
            group_by(1) %>% summarise(nd = sum(nd)) %>% pull(nd)
        stopifnot(nd == 0)
        to = ti2 %>% select(-input_read_pairs) %>%
            gather(type, nseq, -SampleID) 
    } else if (paired == 'both') {
        types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
        ti2 = ti %>% replace_na(list('input_reads'=0, 'input_read_pairs'=0,
                    'forward_only_surviving'=0, 'reverse_only_surviving'=0)) %>%
            mutate(input_read_pairs = ifelse(input_read_pairs == 0, input_reads, input_read_pairs))
        ti2 = ti2 %>% 
            separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
            select(SampleID, input_read_pairs, surviving, forward_only_surviving,
                   reverse_only_surviving, dropped) 
        nd = ti2 %>% mutate(nd = input_read_pairs - surviving - 
                            forward_only_surviving - reverse_only_surviving -
                            dropped) %>%
            group_by(1) %>% summarise(nd = sum(nd)) %>% pull(nd)
        stopifnot(nd == 0)
        to = ti2 %>% select(-input_read_pairs) %>%
            gather(type, nseq, -SampleID) 
    } else {
        stop(sprintf("unsupported option: %s", paired))
    }
    to %>% mutate(nseq = nseq/1000000,
                  type = factor(type, levels = types))
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
                  n.diff = total - uniquely_mapped - multimapped - unmapped) 
    stopifnot(sum(ti2$n.diff) < 1000)
    types = c("uniquely_mapped", "multimapped", "unmapped")
    to = ti2 %>% select(-n.diff, -total) %>%
        gather(type, rc, -SampleID) %>%
        group_by(SampleID, type) %>% 
        summarise(rc = sum(rc) / 1000000) %>%
        mutate(type = factor(type, levels = types))
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
    to = ti2 %>% select(-Total, -nd) %>%
        gather(type, rc, -SampleID) %>%
        mutate(rc = rc / 1000000) %>%
        mutate(type = factor(type, levels = types))
    to
    #}}}
}
