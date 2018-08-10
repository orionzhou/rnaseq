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

dirp = '~/projects/maize.expression'
studies = c(
    'li2013',
    'hirsch2014',
    'leiboff2015',
    'stelpflug2016',
    'walley2016',
    'jin2016', 
    'lin2017',
    'baldauf2018',
    'kremling2018',
    'kaeppler2018'
)
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


