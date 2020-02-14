source("functions.R")

yid = 'rn14f'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in, filter/fix samples
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(Genotype, Treatment, Replicate, sep='_'))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.25)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=6, height=8)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.25)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=6, height=8)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Treatment', var.col='Genotype', var.lab='Replicate',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=3,iter=1200, seed=2,
    var.shape='Treatment', var.col='Genotype', var.lab='Replicate',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), width=6, height=6)
#}}}

th2 = res$th %>% filter(SampleID != 'SRR1238724')
th2 = complete_sample_list(th2)

fh = file.path(dirw, 'meta.tsv')
write_tsv(th2, fh, na='')
#}}}

res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = th %>% mutate(lab = str_c(Genotype, Treatment, Replicate, sep='_'))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.25)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=6, height=8)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.25)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=6, height=8)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Treatment', var.col='Genotype', var.lab='Replicate',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=3,iter=1200, seed=42,
    var.shape='Treatment', var.col='Genotype', var.lab='Replicate',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=6, height=6)
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
    facet_wrap(~lab, ncol=5, scale='free_y') +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T)
fo = file.path(dirw, 'afs_gene.pdf')
ggsave(fo, p, width=6, height=6)
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
    facet_grid(Treatment ~ Genotype) +
    otheme(xtext=T, ytext=T, xtick=T, ytick=T)
fo = file.path(dirw, 'afs_site.pdf')
ggsave(fo, p, width=8, height=6)
#}}}


