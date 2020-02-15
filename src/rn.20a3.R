source("functions.R")

yid = 'rn20a3'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in, filter/fix samples
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

hr_map = c("11am"=11,'12pm'=12,'2pm'=14,'4pm'=16)
th = res$th %>%
    separate(Treatment,c('day','hr','Treatment'), sep='_') %>%
    mutate(day = as.integer(str_sub(day,4))) %>%
    mutate(hr = hr_map[hr]) %>%
    mutate(Timepoint = (day-9)*24 + hr - 11) %>%
    select(SampleID,Tissue,Genotype,Treatment,Timepoint) %>%
    mutate(grp=str_c(Genotype,Treatment,sep='_')) %>%
    mutate(tis=ifelse(Tissue=='leaf_2', '', Tissue)) %>%
    mutate(lab=str_c(Genotype,Treatment,Timepoint,tis,sep='_')) %>%
    mutate(Timepoint = as.character(Timepoint))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='grp',
    expand.x=.3)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=6, height=10)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='grp',
    expand.x=.3)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=6, height=10)

p2 = plot_pca(tm,th,pct.exp=.5, pca.center=T, pca.scale=F,
    var.shape='grp', var.col='grp', var.lab='Timepoint',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.5,perp=3,iter=1200, seed=2,
    var.shape='grp', var.col='grp', var.lab='Timepoint',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), p3, width=6, height=6)
#}}}

th2 = th %>% filter(!SampleID %in% c("B03",'C03','C04','D03','E06','G02','H04'))
th2 = complete_sample_list(th2)

fh = file.path(dirw, 'meta.tsv')
write_tsv(th2, fh, na='')
#}}}

res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='grp',
    expand.x=.3)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=6, height=10)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='grp',
    expand.x=.3)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=6, height=10)

p2 = plot_pca(tm,th,pct.exp=.5, pca.center=T, pca.scale=F,
    var.shape='grp', var.col='grp', var.lab='Timepoint',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=6, height=6)

p3 = plot_tsne(tm,th,pct.exp=.5,perp=3,iter=1200, seed=42,
    var.shape='grp', var.col='grp', var.lab='Timepoint',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=6, height=6)
#}}}


