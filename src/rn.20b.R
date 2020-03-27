source("functions.R")

yid = 'rn20b'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

thl = th %>% arrange(Tissue,Genotype,Treatment) %>%
    select(Tissue, Genotype, Treatment, Replicate) %>%
    group_by(Tissue, Genotype, Treatment) %>% slice(1) %>% ungroup() %>%
    mutate(elab = Genotype)
th = res$th %>%
    mutate(grp = str_c(Tissue, Genotype, Treatment, sep="_")) %>%
    left_join(thl, by=c('Tissue','Genotype','Treatment','Replicate')) %>%
    replace_na(list(elab='')) %>%
    mutate(lab = str_c(Treatment, Tissue, Genotype, Replicate, sep='_'))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    expand.x=.25)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=8, height=15)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    expand.x=.25)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=8, height=15)

p2 = plot_pca(tm,th[th$Treatment=='m',],pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Tissue', var.col='Genotype', var.lab='elab', var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.m.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th[th$Treatment=='m',],pct.exp=.7,perp=4,iter=1000, seed=42,
    var.shape='Tissue', var.col='Genotype', var.lab='elab', var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.m.pdf'), width=6, height=6)

p2 = plot_pca(tm,th[th$Treatment=='t',],pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Tissue', var.col='Genotype', var.lab='elab', var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.t.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th[th$Treatment=='t',],pct.exp=.7,perp=4,iter=1000, seed=42,
    var.shape='Tissue', var.col='Genotype', var.lab='elab', var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.t.pdf'), width=6, height=6)
#}}}

#{{{ fix
th2 = res$th %>%
    mutate(Treatment=ifelse(SampleID=='mBR570', 't',Treatment)) %>%
    mutate(Replicate=ifelse(SampleID=='mBR570', 5, Replicate))
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}


#{{{ read in
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

thl = th %>% arrange(Tissue,Genotype,Treatment) %>%
    select(Tissue, Genotype, Treatment, Replicate) %>%
    group_by(Tissue, Genotype, Treatment) %>% slice(1) %>% ungroup() %>%
    mutate(elab = Genotype)
th = res$th %>%
    mutate(grp = str_c(Tissue, Genotype, Treatment, sep="_")) %>%
    left_join(thl, by=c('Tissue','Genotype','Treatment','Replicate')) %>%
    replace_na(list(elab='')) %>%
    arrange(Tissue, Genotype) %>%
    mutate(lab = str_c(Treatment, Tissue, Genotype, Replicate, sep='_'))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.25)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=8, height=15)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.25)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=8, height=15)

p2 = plot_pca(tm,th[th$Treatment=='m',],pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Tissue', var.col='Genotype', var.lab='elab', var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.m.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th[th$Treatment=='m',],pct.exp=.7,perp=4,iter=1000, seed=42,
    var.shape='Tissue', var.col='Genotype', var.lab='elab', var.ellipse='grp',
    legend.pos='bottom.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.m.pdf'), width=6, height=6)

p2 = plot_pca(tm,th[th$Treatment=='t',],pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Tissue', var.col='Genotype', var.lab='elab', var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.t.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th[th$Treatment=='t',],pct.exp=.7,perp=4,iter=1000, seed=42,
    var.shape='Tissue', var.col='Genotype', var.lab='elab', var.ellipse='grp',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.t.pdf'), width=6, height=6)
#}}}

#{{{ ase
pa1 = plot_ase(res$ase_gene, th, val.col='Genotype', pal.col='aaas')
fo = file.path(dirw, '31.afs_gene.pdf')
ggsave(fo, pa1, width=7, height=15)

pa2 = plot_ase(res$ase_snp, th, val.col='Genotype', pal.col='aaas')
fo = file.path(dirw, '32.afs_site.pdf')
ggsave(fo, pa2, width=7, height=15)
#}}}

