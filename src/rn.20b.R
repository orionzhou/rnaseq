source("functions.R")
gcfg = read_genome_conf()
tg = gcfg$gene %>% select(gid, ttype)

yid = 'rn20b'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in, filter/fix samples
ref = t_cfg %>% filter(yid == !!yid) %>% pull(ref)
th = rnaseq_sample_meta(yid)
tt = rnaseq_mapping_stat(yid)
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
sum_stat_tibble(tt)

#{{{ plot read number
tt2 = tm %>% select(gid, SampleID, ReadCount) %>%
    left_join(tg, by='gid') %>%
    mutate(ttype = ifelse(ttype == 'mRNA', ttype, 'non_mRNA')) %>%
    group_by(SampleID, ttype) %>%
    summarise(ReadCount = sum(ReadCount)) %>% ungroup() %>%
    spread(ttype, ReadCount) %>% print(n=50)

tt3 = tt %>% select(SampleID, Total=total, Trimmed=passed, Mapped=mapped, UniquelyMapped=mappedu) %>%
    inner_join(tt2, by='SampleID') %>%
    mutate(mRNA = mRNA/1e6, non_mRNA=non_mRNA/1e6) %>%
    mutate(noFeature = UniquelyMapped - mRNA - non_mRNA,
        multiMapped = Mapped - UniquelyMapped,
        unMapped = Trimmed - Mapped,
        failedQC = Total - Trimmed) %>%
    select(SampleID, failedQC, unMapped, multiMapped, noFeature, non_mRNA, mRNA) %>%
    gather(opt, rc, -SampleID)

#{{{ bar plot
gts = th %>% distinct(Genotype) %>% pull(Genotype)
opts = tt3 %>% distinct(opt) %>% pull(opt)
tp = tt3 %>%
    inner_join(th, by='SampleID') %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, opt, rc) %>%
    mutate(txt=str_c(Tissue,Genotype,Replicate, sep='-')) %>%
    mutate(Genotype = factor(Genotype, levels=gts)) %>%
    mutate(opt = factor(opt, levels=opts)) %>%
    mutate(Treatment = ifelse(Treatment=='m', 'poly-A mRNA', 'RiboZero RNA'))

ytit = 'Million Read Pairs'
p = ggplot(tp) +
    geom_bar(aes(txt, rc, fill=opt), stat='identity', position='stack', width=.8) +
    scale_x_discrete(expand=expand_scale(mult=c(.01,.01))) +
    scale_y_continuous(name=ytit, expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_npg() +
    facet_wrap(~Treatment, ncol = 1, scale='free') +
    otheme(ytitle=T, ytick=T, xtext=T, ytext=T, xtick=T, legend.pos='right') +
    theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1))
fo = file.path(dirw, 'nread.pdf')
ggsave(fo, p, width=13, height=10)
#}}}

#{{{ heatmap - obsolete
gts = th %>% distinct(Genotype) %>% pull(Genotype)
opts = tt3 %>% distinct(opt) %>% pull(opt)
tp = tt3 %>%
    inner_join(th, by='SampleID') %>%
    mutate(Tissue=str_c(Tissue,Replicate, sep='-')) %>%
    select(-Treatment, -Replicate, -MergeID, -paired, -spots, -avgLength) %>%
    mutate(Genotype = factor(Genotype, levels=gts)) %>%
    mutate(opt = factor(opt, levels=opts))
tps = tp %>% group_by(opt) %>% summarise(rc_max = max(rc)) %>% ungroup()
tp = tp %>% inner_join(tps, by=c('opt')) %>%
    mutate(prop = rc / rc_max)

p = ggplot(tp) +
    geom_tile(aes(Tissue,Genotype,fill=prop), color='black') +
    geom_text(aes(Tissue,Genotype,label=number(rc,accuracy=.1), color=rc<.5), hjust=1, size=3, nudge_x=.3) +
    scale_x_discrete(expand=expand_scale(mult=c(.01,.01))) +
    scale_y_discrete(expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    facet_wrap(~opt, ncol = 2) +
    otheme(xtext=T, ytext=T, xtick=T, legend.pos='none') +
    theme(axis.text.x=element_text(angle=30, hjust=1, vjust=1)) +
    theme(panel.border = element_blank()) +
    ggtitle('Million Read Pairs') +
    theme(plot.title=element_text(hjust=.5))
fo = file.path(dirw, 'nread.m.pdf')
ggsave(fo, p, width=8, height=7)
#}}}
#}}}

min_mappedu = 5
tt %>% filter(mappedu<min_mappedu)
sids_keep = tt %>% filter(mappedu>=min_mappedu) %>% pull(SampleID)
sum_stat_tibble(tt %>% filter(SampleID %in% sids_keep))

# fix th
th2 = th %>% filter(SampleID %in% sids_keep) %>%
    mutate(Treatment=ifelse(SampleID=='mBR570', 't',Treatment)) %>%
    mutate(Replicate=ifelse(SampleID=='mBR570', 5, Replicate))

th = th2
tt = tt %>% filter(SampleID %in% th$SampleID)

fh = file.path(dirw, 'meta.tsv')
write_tsv(th, fh, na='')
# run snakemake again
#}}}

res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
thl = th %>% group_by(Tissue, Genotype, Treatment) %>%
    summarise(groupLab = Genotype[1], Replicate = Replicate[1]) %>% ungroup()
th = th %>%
    mutate(group = str_c(Tissue, Genotype, Treatment, sep="_")) %>%
    left_join(thl, by=c('Tissue','Genotype','Treatment','Replicate')) %>%
    replace_na(list(groupLab='')) %>%
    mutate(lab = sprintf("%s_%s_%s_%d", Treatment, Tissue, Genotype, Replicate))

#{{{ hclust
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)

cor_opt = "spearman"
cor_opt = "pearson"
hc_opt = "ward.D"
hc_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
edist <- as.dist(1-cor(e, method = cor_opt))
ehc <- hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
lnames = ehc$labels[ehc$order]

tp = th %>% mutate(taxa = SampleID) %>%
    select(taxa, everything())
p1 = ggtree(tree, layout = 'rectangular') +
    scale_x_continuous(expand = expand_scale(mult=c(0,.2))) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+%
    tp + geom_tiplab(aes(label=lab, color=Genotype), size=2) +
    scale_color_aaas()
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width=8, height=12)
#}}}

suf='t'
sids = th %>% filter(Treatment==suf) %>% pull(SampleID)
#{{{ tSNE
require(Rtsne)
tw = tm %>% filter(SampleID %in% sids) %>%
    select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=4,
              pca = T, max_iter = 1000)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp, aes(x=V1,y=V2)) +
    geom_mark_ellipse(aes(fill=group,label=Genotype),
        expand=unit(2,'mm'), alpha=0, size = .2,
        con.type='none',label.fontsize=0,label.minwidth=unit(0,'mm'),
        label.buffer=unit(0,'mm'),label.margin = margin(0,0,0,0,"mm")) +
    geom_point(aes(x=V1, y=V2, color=Genotype, shape=Tissue), size=2) +
    geom_text_repel(aes(label=groupLab), size=3) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(0,1,4)) +
    scale_color_aaas() +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(color=F, fill=F)
fp = sprintf("%s/25.tsne.%s.pdf", dirw, suf)
ggsave(p_tsne, filename = fp, width=6, height=6)
#}}}

#{{{ PCA
tw = tm %>% select(SampleID, gid, rFPKM) %>% spread(SampleID, rFPKM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(rFPKM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .5) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]

xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
tismap = LETTERS[1:length(tissues23)]
names(tismap) = tissues23
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    inner_join(th, by = 'SampleID') %>% 
    mutate(lab = tismap[Tissue])
cols17 = c(brewer.pal(8, 'Dark2'), brewer.pal(9, 'Set1'))
cols4 = c(brewer.pal(4, 'Set1'))
p1 = ggplot(tp) +
    geom_point(aes(x = PC1, y = PC2, shape = Tissue, color = Genotype), size = 3) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    scale_shape_manual(values = as.character(tismap)) +
    scale_color_manual(values = cols4) +
    theme_bw() +
    theme(axis.ticks.length = unit(0, 'lines')) +
    theme(plot.margin = unit(c(.2,.2,.2,.2), "lines")) +
    theme(legend.position = 'right', legend.direction = "vertical", legend.justification = c(0,.5)) +
    theme(legend.title = element_blank(), legend.key.size = unit(1, 'lines'), legend.key.width = unit(1, 'lines'), legend.text = element_text(size = 9)) +
    guides(shape = guide_legend(ncol = 1, byrow = T))
fp = sprintf("%s/09.pca.pdf", dirw)
ggsave(p1, filename = fp, width = 8, height = 8)
#}}}

#{{{ ase gene
fi = file.path(dird, 'raw', yid, 'ase.rds')
ti = readRDS(fi)

tp = ti %>% filter(allele1 + allele2 >= 20) %>%
    mutate(af = allele1/(allele1 + allele2)) %>%
    inner_join(th, by=c('sid'='SampleID'))
tp %>% group_by(lab) %>%
    summarise(q50=median(af), m50=sum(allele1)/sum(allele1+allele2)) %>%
    ungroup() %>% print(n=40)
p = ggplot(tp) +
    geom_histogram(aes(af), binwidth=.02) +
    geom_vline(xintercept = .5, color='red') +
    scale_y_continuous(expand=expand_scale(mult=c(0,.03))) +
    facet_wrap(~lab, ncol=8, scale='free_y') +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T)
fo = file.path(dirw, 'afs_gene.pdf')
ggsave(fo, p, width=10, height=10)
#}}}


