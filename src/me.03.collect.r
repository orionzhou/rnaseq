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

sid = 'me99b'
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
study = 'me99d'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
gts = sort(unique(th$Genotype))
tissues = sort(unique(th$Tissue))

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>% 
    transmute(SampleID = Sample, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~Tissue, ncol = 1, scale = 'free') +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)
#
types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Genotype, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Tissue, ncol = 1, scale = 'free') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

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
study = 'me99e'
dirw = file.path(dirp, study, 'data')
diri = file.path(dirp, study, 'data/raw/multiqc_data')
fh = file.path(dirw, '01.reads.tsv')
th = read_tsv(fh)
th = th %>% mutate(Genotype = Treatment)
tissues = sort(unique(th$Tissue))
gts = sort(unique(th$Genotype))

#{{{ trim - number read pairs - w. both SE and PE
fi = file.path(diri, "multiqc_trimmomatic.txt")
ti = read_tsv(fi)
ti2 = ti %>% separate(Sample, c("SampleID", "suf"), sep = "_") %>% select(-suf) %>%
    replace_na(list('input_reads'=0, 'input_read_pairs'=0,
                    'forward_only_surviving'=0, 'reverse_only_surviving'=0)) %>%
    mutate(input_read_pairs = ifelse(input_read_pairs == 0, input_reads, input_read_pairs))
ti2 %>% mutate(ndiff = input_read_pairs - surviving - forward_only_surviving - reverse_only_surviving - dropped) %>%
    group_by(1) %>% summarise(ndiff = sum(ndiff))
ti3 = ti2 %>% 
    select(SampleID, surviving, forward_only_surviving,
           reverse_only_surviving, dropped) %>%
    gather(type, nseq, -SampleID)

types = c("surviving", "forward_only_surviving", "reverse_only_surviving", "dropped")
tp = th %>% inner_join(ti3, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(nseq = sum(nseq)/1000000) %>% ungroup() %>% 
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = rev(tissues)),
           Genotype = factor(Genotype, levels = gts))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = nseq, fill = type), stat = 'identity', position = position_stack(reverse = T), width = .7) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 1, strip.position = 'top', scales = 'free_x') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 7, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4)) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/11.reads.pdf", dirw)
ggsave(p1, filename = fp, width = 6, height = 6)
#}}}

#{{{ star - mapped reads
fi = file.path(diri, 'multiqc_star.txt')
ti = read_tsv(fi)
ti2 = ti %>%
    separate(Sample, c("SampleID", 'suf'), sep = "_") %>% select(-suf) %>%
    transmute(SampleID = SampleID, total = total_reads,
              uniquely_mapped = uniquely_mapped,
              multimapped = multimapped + multimapped_toomany,
              unmapped = unmapped_mismatches + unmapped_tooshort + unmapped_other,
              n.diff = total - uniquely_mapped - multimapped - unmapped)
sum(ti2$n.diff)
ti2 = ti2 %>% gather(type, rc, -SampleID) %>%
    group_by(SampleID, type) %>% summarise(rc = sum(rc)) %>%
    spread(type, rc)
#
types = c("uniquely_mapped", "multimapped", "unmapped")
tp = ti2 %>% select(-total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Tissue, Genotype, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    theme_bw() +
    facet_wrap(~Genotype, ncol = 1, scale = 'free') +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/12.mapping.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ featurecounts
fi = file.path(diri, 'multiqc_featureCounts.txt')
ti = read_tsv(fi)
ti2 = ti %>% mutate(SampleID = Sample) %>%
    select(SampleID, Total, Assigned, Unassigned_Unmapped, Unassigned_MultiMapping,
           Unassigned_NoFeatures, Unassigned_Ambiguity) %>%
    mutate(n.diff = Total-Assigned-Unassigned_Unmapped-Unassigned_MultiMapping-Unassigned_NoFeatures-Unassigned_Ambiguity)
sum(ti2$n.diff)
#
types = c("Assigned", "Unassigned_MultiMapping", "Unassigned_Unmapped",
          "Unassigned_NoFeatures", "Unassigned_Ambiguity")
tp = ti2 %>% select(-Total, -n.diff) %>%
    gather(type, rc, -SampleID) %>%
    inner_join(th, by = 'SampleID') %>%
    group_by(Genotype, Tissue, type) %>%
    summarise(rc = sum(rc)/1000000) %>% ungroup() %>%
    mutate(type = factor(type, levels = types),
           Tissue = factor(Tissue, levels = tissues),
           Genotype = factor(Genotype, levels = rev(gts)))
p1 = ggplot(tp) +
    geom_bar(mapping = aes(x = Tissue, y = rc, fill = type), stat = 'identity', position = position_stack(reverse = T), width=.8) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(name = 'Num. Read Pairs (/million)', expand = expand_scale(mult=c(0,.03))) +
    scale_fill_simpsons() +
    coord_flip() +
    facet_wrap(~Genotype, ncol = 1, scale = 'free') +
    theme_bw() +
    theme(strip.background = element_blank(), strip.text = element_text(size = 8, margin = margin(0,0,0,0,'lines'))) + 
    theme(legend.position = c(.5,1), legend.justification = c(.5,-.4), legend.background = element_blank()) +
    guides(direction = 'horizontal', fill = guide_legend(nrow = 1, byrow = T)) +
    theme(legend.title = element_blank(), legend.key.size = unit(.8, 'lines'), legend.text = element_text(size = 8)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(1.5,.5,.5,.5), "lines")) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank()) +
    theme(axis.title.x = element_text(size = 9), axis.title.y = element_blank()) +
    theme(axis.text = element_text(size=8))
fp = sprintf("%s/13.assigned.pdf", dirw)
ggarrange(p1, nrow = 1, ncol = 1)  %>%
    ggexport(filename = fp, width = 6, height = 6)
#}}}

#{{{ collect featurecounts data & normalize
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

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
#{{{ mec01: stelpflug2016 + walley2016
sid = 'mec01'
dirw = file.path(dird, sid, 'data')
t_cfg %>% filter(sid == !!sid)
sids = str_split(t_cfg %>% filter(sid == !!sid) %>% pull(author), "[\\+]")[[1]] 
sids

#{{{ collect featurecounts data & normalize
res = merge_me_datasets(sids, t_cfg, dird, group = 'Tissue')
th = res$th; t_rc = res$t_rc
dim(th); dim(t_rc)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fh = file.path(dirw, '01.reads.tsv')
write_tsv(th, fh)
fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) + 
    scale_x_continuous(expand = c(0,0), limits=c(-.2,12)) +
    scale_y_discrete(expand = c(.02,0)) +
    theme_tree2()
p1 = p1 %<+% tp + 
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) + 
    geom_text(aes(label = Tissue), size = 2, nudge_x = .1, hjust = 0) +
    geom_text(aes(label = Rep), size = 2, nudge_x = 4, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 12)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]

xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
cols6 = pal_d3()(6)
shapes7 = c(15:18, 7,8,9)
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
tps = tp %>% distinct(Tissue) %>%
    mutate(colors = rep(cols6, each = 7)[1:length(tissues)],
           shapes = rep(shapes7, 6)[1:length(tissues)])
p1 = ggplot(tp) +
    #geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Tissue, group = Tissue), size = 2) +
    geom_point(aes(x = PC1, y = PC2), color = 'red', size = 2) +
    geom_text_repel(aes(x = PC1, y = PC2, label = Tissue), size = 2) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    guides(direction = 'vertical', fill = guide_legend(ncol = 1)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) +
    otheme(xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T) 
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 12, height = 12)
#}}}
#}}}

#{{{ mec02: stelpflug2016 + walley2016 + briggs B73
sid = 'mec02'
dirw = file.path(dird, sid, 'data')
t_cfg %>% filter(sid == !!sid)
sids = str_split(t_cfg %>% filter(sid == !!sid) %>% pull(author), "[\\+]")[[1]] 
sids

#{{{ collect featurecounts data & normalize
res = merge_me_datasets(sids, t_cfg, dird, group = 'Tissue')
th = res$th; t_rc = res$t_rc
th = th %>% filter(Genotype == 'B73')
t_rc = t_rc %>% select(one_of(c('gid', th$SampleID)))
dim(th); dim(t_rc)
tissues = unique(th$Tissue)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fh = file.path(dirw, '01.reads.tsv')
write_tsv(th, fh)
fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           Rep = as.character(Replicate),
           lab = sprintf("%s %s %s", SampleID, Genotype, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) + 
    scale_x_continuous(expand = c(0,0), limits=c(-.5,17)) +
    scale_y_discrete(expand = c(.02,0)) +
    theme_tree2()
p1 = p1 %<+% tp + 
    #geom_tiplab(aes(label = lab), size = 2, offset = 0.04) + 
    geom_text(aes(label = Tissue), size = 2, nudge_x = .1, hjust = 0) +
    geom_text(aes(label = Rep), size = 2, nudge_x = 5, hjust = 0)
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width = 8, height = 18)
#}}}

#{{{ PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]

xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
cols6 = pal_d3()(6)
shapes7 = c(15:18, 7,8,9)
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
tps = tp %>% distinct(Tissue) %>%
    mutate(colors = rep(cols6, each = 7)[1:length(tissues)],
           shapes = rep(shapes7, 6)[1:length(tissues)])
p1 = ggplot(tp) +
    #geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Tissue, group = Tissue), size = 2) +
    geom_point(aes(x = PC1, y = PC2), color = 'red', size = 2) +
    geom_text_repel(aes(x = PC1, y = PC2, label = Tissue), size = 2) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    guides(direction = 'vertical', fill = guide_legend(ncol = 1)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) +
    otheme(xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T) 
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 12, height = 12)
#}}}
#}}}

#{{{ met01: liu2013 + yu2015
sid = 'met01'
dirw = file.path(dird, sid, 'data')
t_cfg %>% filter(sid == !!sid)
sids = str_split(t_cfg %>% filter(sid == !!sid) %>% pull(author), "[\\+]")[[1]] 
sids

#{{{ collect featurecounts data & normalize
res = merge_me_datasets(sids, t_cfg, dird, group = 'Treatment')
th = res$th; t_rc = res$t_rc
th = th %>% filter(!str_detect(Treatment, "Liu2013_ET"))
t_rc = t_rc %>% select(one_of('gid', th$SampleID))
dim(th); dim(t_rc)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

fh = file.path(dirw, '01.reads.tsv')
write_tsv(th, fh)
fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}

#{{{ read data for hclust and pca
fi = file.path(dirw, '20.rc.norm.rda')
x = load(fi)
x

tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
#}}}

#{{{ hclust tree
cor_opt = "pearson"
hc_opt = "ward.D"
e.c.dist <- as.dist(1-cor(e, method = cor_opt))
e.c.hc <- hclust(e.c.dist, method = hc_opt)
hc = e.c.hc
tree = as.phylo(e.c.hc)

tp = tl %>% inner_join(th, by = 'SampleID') %>%
    mutate(taxa = SampleID,
           lab = sprintf("%s %s %s", SampleID, Treatment, Replicate)) %>%
    select(taxa, everything())
p1 = ggtree(tree) + 
    scale_x_continuous(expand = expand_scale(mult=c(0,.5))) +
    scale_y_discrete(expand = c(.02,0)) +
    theme_tree2()
p1 = p1 %<+% tp + 
    geom_tiplab(aes(label = lab), size = 3, offset = 0.04)
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
cols6 = pal_d3()(6)
shapes7 = c(15:18, 7,8,9)
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID')
p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2), color = 'red', size = 2) +
    geom_text_repel(aes(x = PC1, y = PC2, label = Treatment), size = 2.5) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    guides(direction = 'vertical', fill = guide_legend(ncol = 1)) +
    guides(shape = guide_legend(ncol = 1, byrow = T)) +
    otheme(xgrid = T, ygrid = T, xtitle = T, ytitle = T, xtext = T, ytext = T) 
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5))
fp = sprintf("%s/22.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}
#}}}

# write output
#{{{ 
diro = file.path(dird, '09_output')
ngene = 46117

for(i in 1:nrow(t_cfg)) {
    sid = t_cfg$sid[i]; study = t_cfg$study[i]
    cat(sid, study, "\n")
    #
    diri = file.path(dird, sid, 'data')
    fi = file.path(diri, '20.rc.norm.rda')
    if(!file.exists(fi)) next
    x = load(fi)
    fh1 = file.path(diri, '01.reads.tsv')
    fh2 = file.path(diri, '02.reads.corrected.tsv')
    fh = ifelse(file.exists(fh2), fh2, fh1) 
    th = read_tsv(fh)
    tm = tm %>% filter(SampleID %in% th$SampleID)
    #stopifnot(nrow(th) * ngene == nrow(tm))
    #
    tl = th
    fo = sprintf("%s/%s.rda", diro, sid)
    save(tl, tm, file = fo)
}
#}}}

