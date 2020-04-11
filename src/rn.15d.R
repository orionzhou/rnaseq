source("functions.R")

yid = 'rn15d'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) dir.create(dirw)

#{{{ read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(rep1 = Replicate==1) %>%
    mutate(lab = sprintf("%s_%s", Genotype, Replicate))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='rep1',
    pal.col='aaas', expand.x=.1)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=8, height=12)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='rep1',
    pal.col='aaas', expand.x=.1)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=8, height=12)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='', var.col='rep1', var.lab='Genotype', var.ellipse='Genotype',
    legend.pos='none', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=8,iter=1500, seed=2,
    var.shape='', var.col='rep1', var.lab='Genotype', var.ellipse='Genotype',
    legend.pos='none', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), width=8, height=8)
#}}}

#{{{ read in
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(rep1 = Replicate==1) %>%
    mutate(lab = sprintf("%s_%s", Genotype, Replicate))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='rep1',
    pal.col='aaas', expand.x=.1)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=8, height=12)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='rep1',
    pal.col='aaas', expand.x=.1)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=8, height=12)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='', var.col='rep1', var.lab='Genotype', var.ellipse='Genotype',
    legend.pos='none', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=8,iter=1500, seed=2,
    var.shape='', var.col='rep1', var.lab='Genotype', var.ellipse='Genotype',
    legend.pos='none', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), width=8, height=8)
#}}}

gts = c("B73","H99")
sids_red = th %>% filter(Genotype %in% c(gts)) %>% pull(SampleID)
#{{{ RIL haplotype blocks
cp = res$ril$cp
p = plot_ril_genotype(cp, th, sids_red=sids_red, gts=c("B73",'H99','het'))
fo = file.path(dirw, '41.ril.genotype.pdf')
ggsave(fo, p, width=8, height=12)
#}}}

