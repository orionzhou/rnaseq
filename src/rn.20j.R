source("functions.R")

yid = 'rn20j'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) dir.create(dirw)

#{{{ raw: read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(Tissue, Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Genotype=='B73' & Replicate==1, Tissue, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    expand.x=.25)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=6, height=8)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    expand.x=.25)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=6, height=8)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=3,iter=1200, seed=2,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
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
    mutate(lab = str_c(Tissue, Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Genotype=='B73' & Replicate==1, Tissue, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    expand.x=.25)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=6, height=8)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    expand.x=.25)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=6, height=8)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=3,iter=1200, seed=2,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=6, height=6)
#}}}


#{{{ ase
gt = 'A632d_D'
pa1 = plot_ase(res$ase_gene, th, val.col='Genotype', pal.col='aaas')
fo = sprintf("%s/31.afs_gene.%s.pdf", dirw, gt)
ggsave(fo, pa1, width=7, height=7)

pa2 = plot_ase(res$ase_snp, th, val.col='Genotype', pal.col='aaas')
fo = sprintf("%s/32.afs_site.%s.pdf", dirw, gt)
ggsave(fo, pa2, width=7, height=7)
#}}}

#{{{ RIL haplotype blocks
gts = c("B73","Mo17")
gts = c("B73","LH145")
gts = c("B73","A632")
gts = c("B73","A635")
gts = c("B73","LH143")
gts = c("B73","A632d")
sids_red = th %>% filter(Genotype %in% c(gts)) %>% pull(SampleID)
cp = res$ril$cp
p = plot_ril_genotype(cp, th, sids_red=sids_red, gts=c(gts,'het'))
fo = sprintf("%s/41.ril.genotype.%s.pdf", dirw, gts[2])
ggsave(fo, p, width=8, height=7)
#}}}


