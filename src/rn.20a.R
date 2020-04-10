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

#{{{ sample clustering
#{{{ TC
ex = 'TC'
th1 = th %>% filter(Experiment==ex) %>%
    mutate(Genotype = factor(Genotype, levels=gts3)) %>%
    mutate(Treatment = factor(Treatment, levels=c("Control",'Cold','Heat'))) %>%
    mutate(has = sprintf("h%03d", Timepoint*10)) %>%
    mutate(grp = sprintf("%s_%s_%s", Treatment, has, Genotype)) %>%
    mutate(cond = sprintf("%s_%s", Treatment, has)) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(Genotype, Treatment, Timepoint) %>%
    mutate(lab = cond)
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1, width=6, height=8)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1, width=6, height=8)

p2 = plot_tsne(tm1,th1,pct.exp=.5,perp=4,iter=900, seed=12,
    var.shape='Treatment',var.col='Genotype',var.lab='has',#var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.pdf", dirw, ex), p2, width=6, height=6)

p3 = plot_umap(tm1,th1,pct.exp=.5,seed=42,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    pal.col='aaas')
ggsave(sprintf("%s/22.umap.%s.pdf", dirw, ex), p3, width=6, height=6)
#}}}

#{{{ HY
ex = 'HY'
th1 = th %>% filter(Experiment==ex) %>%
    mutate(Genotype = factor(Genotype, levels=gts6)) %>%
    mutate(grp = sprintf("%s_%dh_%s", Treatment, Timepoint, Genotype)) %>%
    mutate(cond = sprintf("%s_%dh", Treatment, Timepoint)) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(Genotype, Treatment, Timepoint) %>%
    mutate(Replicate = 1:n()) %>% ungroup() %>%
    mutate(lab = str_c(SampleID, grp, sep=' ')) %>%
    mutate(clab = ifelse(Replicate==1, cond, ''))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Experiment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1, width=8, height=10)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Experiment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1, width=8, height=10)

p2 = plot_tsne(tm1,th1,pct.exp=.5,perp=4,iter=900, seed=12,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0:3,5,6), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.pdf", dirw, ex), p2, width=6, height=6)

p3 = plot_umap(tm1,th1,pct.exp=.5,seed=42,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    pal.col='aaas')
ggsave(sprintf("%s/22.umap.%s.pdf", dirw, ex), p3, width=6, height=6)
#}}}

#{{{ NM
ex = 'NM'
th1 = th %>% filter(Experiment==ex) %>%
    mutate(Genotype = factor(Genotype, levels=gts25)) %>%
    mutate(Treatment = factor(Treatment, levels=c("Control",'Cold'))) %>%
    mutate(grp = sprintf("%s_%dh_%s", Treatment, Timepoint, Genotype)) %>%
    mutate(cond = sprintf("%s_%dh", Treatment, Timepoint)) %>%
    mutate(cond = factor(cond, levels=c("Control_1h",'Control_25h','Cold_1h','Cold_25h'))) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(Genotype, Treatment, Timepoint) %>%
    mutate(lab = str_c(cond, Genotype, sep="_"))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1, width=6, height=10)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1, width=6, height=10)

p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=4,iter=1200, seed=12,
    var.shape='cond',var.col='cond',var.lab='Genotype',#var.ellipse='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0:1,15:16), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.pdf", dirw, ex), p2, width=6, height=6)
#}}}
#}}}

#{{{ ase
pa1 = plot_ase(res$ase_gene, th, val.col='Genotype', pal.col='')
fo = file.path(dirw, '31.afs_gene.pdf')
ggsave(fo, pa1, width=7, height=25)

pa2 = plot_ase(res$ase_snp, th, val.col='Genotype', pal.col='')
fo = file.path(dirw, '32.afs_site.pdf')
ggsave(fo, pa2, width=7, height=25)
#}}}



