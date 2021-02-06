source("functions.R")
genome = 'Osativa'
t_cfg = read_projects(genome)
t_cfg %>% print(n=50)

yid = 'mo20a3'
dirw = glue("{dird}/11_qc/{genome}/{yid}")
if(!dir.exists(dirw)) dir.create(dirw)

#{{{ raw: read in
res = rnaseq_cpm_raw(yid, genome)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th2 = read_xlsx(glue("~/projects/cold/data/samples.xlsx")) %>%
    replace_na(list(Name='?', ColdTolerant='?'))
th = res$th %>% left_join(th2, by="Genotype") %>%
    mutate(lab = str_c(Genotype, Name, Treatment, ColdTolerant, sep='_'))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.4)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=6, height=8)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.4)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=6, height=8)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='ColdTolerant', var.col='Treatment', var.lab='Name',
    legend.pos='bottom.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.8,perp=3,iter=1200, seed=2,
    var.shape='ColdTolerant', var.col='Treatment', var.lab='Name',
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

th2 = read_xlsx(glue("~/projects/cold/data/samples.xlsx")) %>%
    replace_na(list(Name='?', ColdTolerant='?'))
th = res$th %>% left_join(th2, by="Genotype") %>%
    mutate(lab = str_c(Genotype, Name, Treatment, ColdTolerant, sep='_'))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.4)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=6, height=8)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.4)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=6, height=8)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='ColdTolerant', var.col='Treatment', var.lab='Name',
    legend.pos='bottom.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.8,perp=3,iter=1200, seed=2,
    var.shape='ColdTolerant', var.col='Treatment', var.lab='Name',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=6, height=6)
#}}}



