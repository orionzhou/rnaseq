source("functions.R")

yid = 'rn13c'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) dir.create(dirw)

#{{{ raw: read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(Replicate = as.character(Replicate)) %>%
    mutate(lab = str_c(Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Replicate==1, Genotype, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Replicate',
    expand.x=.1)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=6, height=10)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Replicate',
    expand.x=.1)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=6, height=10)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Replicate', var.col='Replicate', var.ellipse='Genotype', var.lab='clab',
    legend.pos='none', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=7,iter=1500, seed=2,
    var.shape='Replicate', var.col='Replicate', var.ellipse='Genotype', var.lab='clab',
    legend.pos='none', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), p3, width=6, height=6)
#}}}


#{{{ raw: filter/fix samples
th2 = res$th %>%
    mutate(Genotype = ifelse(SampleID=='SRR767691','Oh43',Genotype)) %>%
    mutate(Genotype = ifelse(SampleID=='SRR651079','Oh7b',Genotype))

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

#{{{ read in
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(Replicate = as.character(Replicate)) %>%
    mutate(lab = str_c(Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Replicate==1, Genotype, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Replicate',
    expand.x=.1)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=6, height=10)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Replicate',
    expand.x=.1)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=6, height=10)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Replicate', var.col='Replicate', var.ellipse='Genotype', var.lab='clab',
    legend.pos='none', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=7,iter=1500, seed=2,
    var.shape='Replicate', var.col='Replicate', var.ellipse='Genotype', var.lab='clab',
    legend.pos='none', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=6, height=6)
#}}}



