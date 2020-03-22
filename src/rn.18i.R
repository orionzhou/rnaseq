source("functions.R")

yid = 'rn18i'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in, filter/fix samples
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(Genotype, Tissue, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Genotype=='B73' & Replicate==1, Tissue, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    pal.col='aaas', expand.x=.2)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=6, height=12)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    pal.col='aaas', expand.x=.2)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=6, height=12)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
   var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='Tissue',
   legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=4,iter=800, seed=2,
   var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='Tissue',
   legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), p3, width=6, height=6)
#}}}

#{{{ fix
th2 = res$th
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = th %>%
    mutate(lab = str_c(Genotype, Tissue, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Genotype=='B73' & Replicate==1, Tissue, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    pal.col='aaas', expand.x=.2)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=6, height=12)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    pal.col='aaas', expand.x=.2)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=6, height=12)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
   var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='Tissue',
   legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=4,iter=800, seed=2,
   var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='Tissue',
   legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=6, height=6)
#}}}


