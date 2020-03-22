source("functions.R")

yid = 'rn20d'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ raw: read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(SampleID, Tissue, Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Genotype=='B73' & Replicate==1, Tissue, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    expand.x=.3, lab.size=2)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=8, height=35)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    expand.x=.3, lab.size=2)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=8, height=35)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='',var.col='Tissue',var.lab='clab',var.ellipse='Tissue',leg.col=F,
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=8,iter=1200, seed=2,
    var.shape='',var.col='Tissue',var.lab='clab',var.ellipse='Tissue',leg.col=F,
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), p3, width=8, height=8)
#}}}

#{{{ raw: filter/fix samples
th2 = res$th %>%
    mutate(Genotype=str_replace(Genotype, 'xSelf$', '')) %>%
    mutate(Tissue=ifelse(Tissue=='unk', "U", Tissue)) %>%
    mutate(Tissue=ifelse(SampleID=='SRR8043188','L',Tissue)) %>%
    mutate(Tissue=ifelse(SampleID=='SRR8043090','L', Tissue))
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

#{{{ read in
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(Tissue, Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Genotype=='B73' & Replicate==1, Tissue, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    expand.x=.3, lab.size=2)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=8, height=35)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    expand.x=.3, lab.size=2)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=8, height=35)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='',var.col='Tissue',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=8,iter=1500, seed=2,
    var.shape='',var.col='Tissue',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=8, height=8)
#}}}


