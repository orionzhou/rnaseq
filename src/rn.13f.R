source("functions.R")

yid = 'rn13f'
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
th = th %>% mutate(lab = str_c(Genotype, Treatment, sep='_'))

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
    scale_x_continuous(expand = expand_scale(0,.2)) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+%
    tp + geom_tiplab(aes(label=lab, color=Genotype), size=2.5) +
    scale_color_aaas()
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width=6, height=6)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=3,
              pca = T, max_iter = 1000)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp) +
    geom_text_repel(aes(x=V1,y=V2,label=Treatment), size=2.5) +
    geom_point(aes(x=V1, y=V2, color=Genotype, shape=Genotype), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(15,0)) +
    scale_color_aaas() +
    otheme(legend.pos='top.left', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines'))
fp = file.path(dirw, "25.tsne.pdf")
ggsave(p_tsne, filename = fp, width=6, height=6)
#}}}


#{{{ ase gene
fi = file.path(dird, 'raw', yid, 'ase.rds')
ti = readRDS(fi)

tp = ti %>% filter(allele1 + allele2 >= 20) %>%
    mutate(af = allele1/(allele1 + allele2)) %>%
    inner_join(th, by=c('sid'='SampleID'))
tp %>% group_by(lab) %>%
    summarise(q50=median(af), m50=sum(allele1)/sum(allele1+allele2)) %>%
    ungroup() %>% print(n=70)
p = ggplot(tp) +
    geom_histogram(aes(af), binwidth=.02) +
    geom_vline(xintercept = .5, color='red') +
    scale_y_continuous(expand=expand_scale(mult=c(0,.03))) +
    facet_wrap(~lab, ncol=3, scale='free_y') +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T)
fo = file.path(dirw, 'afs_gene.pdf')
ggsave(fo, p, width=8, height=8)
#}}}

#{{{ ase SNP
fi = file.path(dird, 'raw', yid, 'ase2.rds')
ti2 = readRDS(fi)

tp2 = ti2 %>% filter(allele1 + allele2 >= 20) %>%
    mutate(af = allele1/(allele1 + allele2)) %>%
    inner_join(th, by=c('sid'='SampleID'))
tp2 %>% group_by(Treatment,Genotype) %>%
    summarise(q50=median(af), m50=sum(allele1)/sum(allele1+allele2)) %>% ungroup()
p = ggplot(tp2) +
    geom_histogram(aes(af), binwidth=.02) +
    geom_vline(xintercept = .5, color='red') +
    scale_y_continuous(expand=expand_scale(mult=c(0,.03))) +
    facet_wrap(~lab, ncol=3, scale='free') +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T)
fo = file.path(dirw, 'afs_site.pdf')
ggsave(fo, p, width=8, height=6)
#}}}






