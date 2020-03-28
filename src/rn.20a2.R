source("functions.R")

yid = 'rn20a2'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in, filter/fix samples
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

require(lubridate)
fh = '~/projects/stress/data/samples.xlsx'
th0 = read_xlsx(fh, sheet='merged') %>%
    mutate(Time = sprintf("%02d:%02d", hour(Time), minute(Time))) %>%
    mutate(Tissue = 'leaf') %>% rename(Replicate=Rep) %>%
    select(SampleID,Tissue,Genotype,Treatment,Timepoint,Experiment,Replicate)

th = res$th %>% separate(SampleID, c("batch", "coord"), sep="-", remove=F) %>%
    select(SampleID,batch,coord,ExpID=Treatment,paired,spots,avgLength) %>%
    left_join(th0, by=c('ExpID'='SampleID')) %>%
    mutate(batch = str_replace(batch, '[ab]$', '')) %>%
    mutate(batch = str_replace(batch, 'batch', 'b')) %>%
    mutate(lab = str_c(batch, Genotype, Treatment, Timepoint, sep='_')) %>%
    mutate(Treatment = as.character(Treatment)) %>%
    select(SampleID,Tissue,Genotype,Treatment,Timepoint,lab,Replicate,everything())
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.4,cor.opt='pearson',var.col='batch',
    expand.x=.2)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=8, height=24)

p1 = plot_hclust(tm,th,pct.exp=.5,cor.opt='spearman',var.col='batch',
    expand.x=.2)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=8, height=24)

p2 = plot_pca(tm,th,pct.exp=.5, pca.center=T, pca.scale=F,
    var.shape='Treatment', var.col='batch', var.lab='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.5,perp=8,iter=1500, seed=2,
    var.shape='Treatment', var.col='batch', var.lab='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), width=8, height=8)
#}}}

#{{{ fix
th2 = th %>%
    mutate(Genotype=ifelse(ExpID=='TC19', 'Mo17', Genotype)) %>%
    mutate(Genotype=ifelse(ExpID=='TC20', 'B73', Genotype)) %>%
    mutate(lab = str_c(batch, Genotype, Treatment, Timepoint, sep='_'))
th2 %>% filter(ExpID %in% c("TC19",'TC20')) %>% print(width=Inf)
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

#{{{ read
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.4,cor.opt='pearson',var.col='batch',
    expand.x=.2)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=8, height=24)

p1 = plot_hclust(tm,th,pct.exp=.5,cor.opt='spearman',var.col='batch',
    expand.x=.2)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=8, height=24)

p2 = plot_pca(tm,th,pct.exp=.5, pca.center=T, pca.scale=F,
    var.shape='Treatment', var.col='batch', var.lab='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.5,perp=8,iter=1500, seed=42,
    var.shape='Treatment', var.col='batch', var.lab='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=8, height=8)
#}}}




