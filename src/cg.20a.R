source("functions.R")
genome = 'Zmays_B73'
t_cfg = read_projects(genome)

yid = 'cg20a'
dirw = glue("{dird}/11_qc/{genome}/{yid}")
if(!dir.exists(dirw)) dir.create(dirw)

#{{{ raw: read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(Treatment=ifelse(avgLength<50, '2015', Treatment)) %>%
    mutate(Treatment=ifelse(Tissue %in% c("husk",'stem'), 'Maike', Treatment)) %>%
    replace_na(list(Treatment='normal')) %>%
    mutate(grp = str_c(Tissue,Genotype,Treatment, sep='_')) %>%
    mutate(clab = ifelse(Replicate==1, Treatment, '')) %>%
    mutate(lab = str_c(SampleID, Tissue, Genotype, Treatment, Replicate, sep=' '))
tm1 = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.35)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=6, height=8)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.35)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=6, height=8)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Genotype', var.col='Tissue', var.ellipse='grp', var.lab='clab',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=4,iter=1200, seed=2,
    var.shape='Genotype', var.col='Tissue', var.ellipse='grp', var.lab='clab',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), p3, width=6, height=6)
#}}}

#{{{ raw: filter/fix samples
th2 = res$th
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

#{{{ read in
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(Treatment=ifelse(avgLength<50, '2015', Treatment)) %>%
    mutate(Treatment=ifelse(Tissue %in% c("husk",'stem'), 'Maike', Treatment)) %>%
    replace_na(list(Treatment='normal')) %>%
    mutate(grp = str_c(Tissue,Genotype,Treatment, sep='_')) %>%
    mutate(clab = ifelse(Replicate==1, Treatment, '')) %>%
    mutate(lab = str_c(SampleID, Tissue, Genotype, Treatment, Replicate, sep=' '))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.35)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=6, height=8)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.35)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=6, height=8)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Genotype', var.col='Tissue', var.ellipse='grp', var.lab='clab',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.6,perp=3,iter=1500, seed=2,
    var.shape='Genotype', var.col='Tissue', var.ellipse='grp', var.lab='clab',
    legend.pos='bottom.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=6, height=6)
#}}}


#{{{ ase
pa1 = plot_ase(res$ase_gene, th, min_rc=10, val.col='Genotype', pal.col='aaas')
fo = file.path(dirw, '31.afs_gene.pdf')
ggsave(fo, pa1, width=7, height=10)

pa2 = plot_ase(res$ase_snp, th, min_rc=10, val.col='Genotype', pal.col='aaas')
fo = file.path(dirw, '32.afs_site.pdf')
ggsave(fo, pa2, width=7, height=10)
#}}}


yid = 'cg20a2'

#{{{ raw: read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(SampleID, Tissue, Genotype, sep=' '))
tm1 = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ ase
pa1 = plot_ase(res$ase_gene, th, min_rc=10, val.col='SampleID', pal.col='aaas')
fo = file.path(dirw, '36.afs_gene.pdf')
ggsave(fo, pa1, width=7, height=4)

pa2 = plot_ase(res$ase_snp, th, min_rc=10, val.col='SampleID', pal.col='aaas')
fo = file.path(dirw, '36.afs_site.pdf')
ggsave(fo, pa2, width=7, height=4)
#}}}
