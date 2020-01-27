source("functions.R")

yid = 'rn15d'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in, filter/fix samples
ref = t_cfg %>% filter(yid == !!yid) %>% pull(ref)
th = rnaseq_sample_meta(yid)
tt = rnaseq_mapping_stat(yid)
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
sum_stat_tibble(tt)

sids_keep = tt %>% filter(mapped>5) %>% pull(SampleID)
sum_stat_tibble(tt %>% filter(SampleID %in% sids_keep))

# fix th
th2 = th %>% filter(SampleID %in% sids_keep)


th = th2
tt = tt %>% filter(SampleID %in% th$SampleID)

fh = file.path(dirw, 'meta.tsv')
write_tsv(th, fh, na='')
# run snakemake again
#}}}

res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
th = th %>% mutate(lab = sprintf("%s_%s", Genotype, Replicate))

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
#
tp = th %>% mutate(taxa = SampleID) %>%
    select(taxa, everything())
p1 = ggtree(tree, layout = 'rectangular') +
    scale_x_continuous(expand = expand_scale(0,.1)) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+%
    tp + geom_tiplab(aes(label=lab, color=Genotype), size=2.5) +
    scale_color_viridis_d()
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width=8, height=12)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=10,
              pca = T, max_iter = 1500)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp) +
    geom_text_repel(aes(x=V1,y=V2,label=Genotype), size=2.5) +
    geom_point(aes(x=V1, y=V2, color=Genotype), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_color_viridis_d(name = 'Genotype',direction=-1) +
    otheme(legend.pos='top.left', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(color = F)
fp = file.path(dirw, "25.tsne.pdf")
ggsave(p_tsne, filename = fp, width=8, height=8)
#}}}

#{{{ RIL
fw = '~/projects/genome/data/Zmays_B73/15_intervals/20.win11.tsv'
tw = read_tsv(fw)
offs = c(0, cumsum(tw$size)[-nrow(tw)]) + 0:10 * 10e6
tx = tw %>% mutate(off = offs) %>%
    mutate(gstart=start+off, gend=end+off, gpos=(gstart+gend)/2) %>%
    select(rid,chrom,gstart,gend,gpos,off) %>% filter(chrom!='B99')
fi = file.path(dird, 'raw', yid, 'ril.rds')
res = readRDS(fi)

ti = res$cp
tps = th %>% select(sid=SampleID,Genotype,Replicate) %>% arrange(Genotype) %>%
    mutate(y = 1:n())
tp = ti %>% inner_join(tx, by='rid') %>%
    mutate(start=start+off, end=end+off) %>%
    inner_join(tps, by='sid')
tz = ti %>% mutate(size=end-start) %>% group_by(sid, gt) %>%
    summarise(size = sum(size)) %>%
    mutate(total_size = sum(size)) %>%
    mutate(prop = size / total_size) %>% ungroup() %>%
    select(sid, gt, prop) %>% spread(gt, prop) %>%
    replace_na(list(a=0,b=0,h=0)) %>%
    inner_join(tps, by='sid') %>% arrange(Genotype) %>%
    filter(Replicate==1) %>% select(-y)

xmax = max(tp$end)
gts = c("B73","H99")
tcol = pal_startrek()(1)
ty = tps %>% group_by(Genotype) %>% summarise(y=mean(y)) %>% ungroup() %>%
    mutate(col=ifelse(Genotype %in% gts, tcol, 'black')) %>%
    inner_join(tz, by='Genotype') %>% mutate(lab=sprintf("%.1f", a*100))
p = ggplot(tp) +
    geom_rect(aes(xmin=start,xmax=end,ymin=y-.3,ymax=y+.4, fill=gt)) +
    geom_text(data=ty, aes(x=xmax+5e6,y=y, label=lab), hjust=0, size=2.5, color=ty$col) +
    scale_x_continuous(breaks=tx$gpos, labels=tx$chrom, expand=expand_scale(mult=c(.001,.03))) +
    scale_y_continuous(breaks=ty$y, labels=ty$Genotype, expand=expand_scale(mult=c(.001,.001))) +
    scale_fill_manual(values=pal_simpsons()(8)[c(1,2,5)], labels=c('B73','H99','het')) +
    otheme(legend.pos='top.center.out', legend.dir='h',
           xtext=T, xtick=T, ytext=T, ytick=T) +
    theme(axis.text.y = element_text(color=ty$col))
fo = file.path(dirw, 'ril_genotype.pdf')
ggsave(fo, p, width=10, height=10)
#}}}


