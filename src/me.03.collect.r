#{{{ header
source("me.fun.r")
t_cfg
fi = file.path('~/data/genome/B73', "v37/t2.tsv")
t_gs = read_tsv(fi, col_types = 'ccccciic') %>% 
    filter(etype == 'exon') %>% 
    group_by(gid, tid) %>% 
    summarise(size = sum(end - beg + 1)) %>%
    group_by(gid) %>%
    summarise(size = max(size))
#}}}

sid = 'me99c'
#{{{ mapping stats
Sys.setenv(R_CONFIG_ACTIVE = sid)
diri = file.path(dird, '08_raw_output', sid, 'multiqc_data')
dirw = file.path(dird, '11_qc', sid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))
fh1 = sprintf("%s/05_read_list/%s.tsv", dird, sid)
fh2 = sprintf("%s/05_read_list/%s.c.tsv", dird, sid)
fh = ifelse(file.exists(fh2), fh2, fh1)
th = read_tsv(fh)
tiss = unique(th$Tissue); genos = unique(th$Genotype); treas = unique(th$Treatment)
reps = unique(th$Replicate)
#
tt = read_multiqc(diri, th)
fo = file.path(dirw, '10.mapping.stat.tsv')
write_tsv(tt, fo)
#}}}

#{{{ obtain raw read counts, normalize and save
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)
#
tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm
#
fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read from 20.rc.norm.rda
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x
#}}}

#{{{ prepare for hclust and pca 
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)
#
tp = th %>% mutate(taxa = SampleID, lab = SampleID) 
if(length(tiss)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Tissue), lab)
if(length(genos)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Genotype), lab)
if(length(treas)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Treatment), lab)
if(length(reps)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Replicate), lab)
tp = tp %>% select(taxa, everything())
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
plot_hclust_tree(tree, tp, fo, 
                 labsize = config::get("hc.labsize"), 
                 x.expand = config::get("hc.x.expand"),
                 x.off = config::get("hc.x.off"), 
                 wd = config::get("hc.wd"), ht = config::get("hc.ht"))
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(Treatment = factor(Treatment), Replicate = factor(Replicate))
fo = sprintf("%s/22.pca.pdf", dirw)
plot_pca(tp, fo, opt = config::get("pca.opt"), labsize = config::get("pca.labsize"),
         wd = config::get("pca.wd"), ht = config::get("pca.ht"))
#}}}

#{{{ #identify mis-labelled replicate
cls = cutree(hc, h = .01)
tcl = tibble(SampleID = names(cls), grp = as.integer(cls))
th2 = th %>% inner_join(tcl, by = 'SampleID')
th3 = th2 %>% group_by(Genotype) %>% 
    summarise(nrep = length(Replicate), ngrp = length(unique(grp))) %>%
    ungroup() %>%
    filter(nrep > 1, ngrp > 1)
th2 %>% filter(Genotype %in% th3$Genotype) %>% print(n=40)
#}}}

#{{{ generate corrected read list me??.c.tsv
if(sid == 'me14c') {
    th = th %>% filter(! SampleID %in% c("SRR254169"))
} else if(sid == 'me14d') {
    th = th %>% filter(! SampleID %in% c("SRR1573518", 'SRR1573513'))
} else if(sid == 'me17a') {
    th = th %>%
        mutate(Tissue = ifelse(SampleID == 'SRR445601', 'tassel', Tissue)) %>%
        mutate(Tissue = ifelse(SampleID == 'SRR445416', 'tassel', Tissue)) %>%
        mutate(Genotype = ifelse(SampleID == 'SRR426798', 'Mo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID == 'SRR426814', 'M37W', Genotype))
} else if(sid == 'me99b') {
    gts = c("B73", "Mo17", "B73xMo17")
    tissues = sort(unique(th$Tissue))
    th = th %>% 
        filter(Genotype %in% gts, ! SampleID %in% c('BR207', 'BR230', "BR235"))
}
write_tsv(th, fh2, na = '')
#}}}


#{{{ ##me99d - enders stress response 3' RNA-Seq
#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           lab = sprintf("%s %s %s", Tissue, Genotype, Treatment)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(-.02,3.3)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = lab), size = 2.5, nudge_x = .01, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 10)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
p1 = ggplot(tp, aes(x = PC1, y = PC2, shape = Genotype, color = Tissue)) +
    geom_point(size = 1.5) +
    #geom_label_repel() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_d3() +
    #scale_shape_manual(values = c(16, 4)) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(1,1), legend.justification = c(1,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}
#}}}

#{{{ ##me99e - settles endosperm
#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           lab = sprintf("%s %s %s", Tissue, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) +
    #geom_tiplab(size = 4, color = 'black') +
    scale_x_continuous(expand = c(0,0), limits=c(-.02,5.5)) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) +
    geom_text(aes(label = lab), size = 2.5, nudge_x = .01, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 6, height = 8)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
p1 = ggplot(tp, aes(x = PC1, y = PC2, shape = Genotype, color = Tissue)) +
    geom_point(size = 1.5) +
    #geom_label_repel() +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_color_d3() +
    #scale_shape_manual(values = c(16, 4)) +
    theme_bw() +
    #theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = c(0,1), legend.justification = c(0,1)) +
    theme(legend.direction = "vertical", legend.background = element_blank()) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.key.width = unit(.8, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}
#}}}

# combined datasets
sid = 'mec03'
#{{{ collect featurecounts data & normalize
Sys.setenv(R_CONFIG_ACTIVE = sid)
sids = str_split(config::get("sids"), "[\\+]")[[1]] 
sids
dirw = file.path(dird, '11_qc', sid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))
#
th = tibble(); t_rc = tibble()
for (sid1 in sids) {
    diri = file.path(dird, '08_raw_output', sid1, 'multiqc_data')
    fh1 = sprintf("%s/05_read_list/%s.tsv", dird, sid1)
    fh2 = sprintf("%s/05_read_list/%s.c.tsv", dird, sid1)
    fh = ifelse(file.exists(fh2), fh2, fh1) 
    th1 = read_tsv(fh) %>% mutate(sid = sid1) %>% select(sid, everything())
    if(sid1 == 'me99b') {
        th1 = th1 %>% filter(Genotype == 'B73')
    } else if(sid1 == 'me13b') {
        th1 = th1 %>% filter(!str_detect(Treatment, "ET"))
    } else if(sid1 == 'me12a') {
        th1 = th1 %>% filter(Treatment == 'WT')
    }
    th = rbind(th, th1)
    fi = file.path(diri, '../featurecounts.tsv')
    t_rc1 = read_tsv(fi) %>% select(one_of(c('gid', th1$SampleID)))
    stopifnot(ncol(t_rc1) == nrow(th1) + 1)
    if(nrow(t_rc) == 0)
        t_rc = t_rc1
    else {
        t_rc = t_rc %>% inner_join(t_rc1, by = 'gid')
    }
}
dim(th); dim(t_rc)
th %>% dplyr::count(sid)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm
#}}}

#{{{ merge reps and save to 20.rc.norm.rda
ths = th %>% distinct(sid, Tissue, Genotype, Treatment) %>%
    mutate(nSampleID = sprintf("%s_%03d", !!sid, 1:length(Tissue)))
th = th %>% inner_join(ths, by = c("sid", "Tissue", "Genotype", "Treatment"))
t_map = th %>% select(SampleID, nSampleID)
th = ths %>% select(SampleID=nSampleID, sid, Tissue, Genotype, Treatment)
tiss = unique(th$Tissue); genos = unique(th$Genotype); treas = unique(th$Treatment)
#
tm = tm %>% inner_join(t_map, by = 'SampleID') %>%
    mutate(SampleID = nSampleID) %>%
    group_by(gid, SampleID) %>%
    #summarise(ReadCount = sum(ReadCount), nRC = sum(nRC), rCPM = mean(rCPM),
    #          rFPKM = mean(rFPKM), CPM = mean(CPM), FPKM = mean(FPKM)) %>%
    summarise(CPM = mean(CPM), FPKM = mean(FPKM)) %>%
    ungroup()
#
fo = file.path(dirw, '20.rc.norm.rda')
save(th, tm, file = fo)
#}}}

#{{{ read from 20.rc.norm.rda
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x
#}}}

#{{{ prepare for hclust and pca 
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust
cor_opt = "pearson"
hc_opt = "ward.D"
plot_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)
#
tp = th %>% mutate(taxa = SampleID, lab = sid)
if(length(tiss)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Tissue), lab)
if(length(genos)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Genotype), lab)
if(length(treas)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Treatment), lab)
tp = tp %>% select(taxa, everything())
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
plot_hclust_tree(tree, tp, fo, 
                 labsize = config::get("hc.labsize"), 
                 x.expand = config::get("hc.x.expand"),
                 x.off = config::get("hc.x.off"), 
                 wd = config::get("hc.wd"), ht = config::get("hc.ht"))
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(Treatment = factor(Treatment))
fo = sprintf("%s/22.pca.pdf", dirw)
plot_pca(tp, fo, opt = config::get("pca.opt"), labsize = config::get("pca.labsize"),
         wd = config::get("pca.wd"), ht = config::get("pca.ht"))
#}}}

# write output
#{{{ 
diro = file.path(dird, '15_output')
for(i in 1:nrow(t_cfg)) {
    sid = t_cfg$sid[i]; study = t_cfg$study[i]
    #
    diri = file.path(dird, '11_qc', sid)
    fi = file.path(diri, '20.rc.norm.rda')
    if(!file.exists(fi)) next
    cat(sid, study, "\n")
    x = load(fi)
    if(!str_detect(sid, 'me[ct]')) {
        th = get_read_list(dird, sid)
        tm = tm %>% filter(SampleID %in% th$SampleID)
    }
    #stopifnot(nrow(th) * ngene == nrow(tm))
    #
    tl = th
    fo = sprintf("%s/%s.rda", diro, sid)
    save(tl, tm, file = fo)
}
#}}}

