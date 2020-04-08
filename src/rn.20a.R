source("functions.R")

yid = 'rn20a'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) dir.create(dirw)

#{{{ read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    arrange(Experiment, Genotype, Treatment, Timepoint) %>%
    mutate(lab = str_c(Experiment, Genotype, Treatment, Timepoint, Replicate, sep='_'))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Experiment',
    expand.x=.2)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=8, height=24)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Experiment',
    expand.x=.2)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=8, height=24)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Treatment', var.col='Experiment', var.lab='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=8,iter=1500, seed=2,
    var.shape='Treatment', var.col='Experiment', var.lab='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), width=8, height=8)
#}}}

#{{{ switch TC19 and TC20
th2 = res$th %>%
    mutate(Genotype=ifelse(SampleID=='TC19', 'Mo17', Genotype)) %>%
    mutate(Genotype=ifelse(SampleID=='TC20', 'B73', Genotype))
th2 %>% filter(SampleID %in% c("TC19",'TC20')) %>% print(width=Inf)
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

#{{{ read in again
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    arrange(Experiment, Genotype, Treatment, Timepoint) %>%
    mutate(lab = str_c(Experiment, Genotype, Treatment, Timepoint, Replicate, sep='_'))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Experiment',
    expand.x=.2)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=8, height=24)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Experiment',
    expand.x=.2)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=8, height=24)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Treatment', var.col='Experiment', var.lab='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=8,iter=1500, seed=2,
    var.shape='Treatment', var.col='Experiment', var.lab='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), width=8, height=8)
#}}}

#{{{ ase
pa1 = plot_ase(res$ase_gene, th, val.col='Genotype', pal.col='')
fo = file.path(dirw, '31.afs_gene.pdf')
ggsave(fo, pa1, width=7, height=25)

pa2 = plot_ase(res$ase_snp, th, val.col='Genotype', pal.col='')
fo = file.path(dirw, '32.afs_site.pdf')
ggsave(fo, pa2, width=7, height=25)
#}}}



