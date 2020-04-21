source("functions.R")
gts3 = c("B73",'Mo17','W22')
gts6 = c("B73",'Mo17','W22','B73xMo17','W22xB73','W22xMo17')
gts25 = c("B73", "B97", "CML322", "CML333", "CML52", "CML69", "DK105",
    "EP1", "F7", "Il14H", "Ki11", "Ki3", "M162W", "M37W",
    "Mo17", "Mo18W", "MS71", "NC350", "NC358", "Oh43", "Oh7B",
    "P39", "PH207", "Tx303", "W22")

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
    mutate(lab = str_c(SampleID, grp, sep=' '))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/11.hclust.%s.p.pdf",dirw,ex), p1, width=6, height=8)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/11.hclust.%s.s.pdf",dirw,ex), p1, width=6, height=8)

p2 = plot_tsne(tm1,th1,pct.exp=.7,perp=5,iter=900, seed=12,
    var.shape='Treatment',var.col='Genotype',var.lab='has',#var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/12.tsne.%s.pdf", dirw, ex), p2, width=6, height=6)
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

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/11.hclust.%s.p.pdf",dirw,ex), p1, width=6, height=10)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/11.hclust.%s.s.pdf",dirw,ex), p1, width=6, height=10)

p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=6,iter=800, seed=12,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0:3,5,6), pal.col='aaas')
ggsave(sprintf("%s/12.tsne.%s.pdf", dirw, ex), p2, width=6, height=6)
#}}}

#{{{ NM
ex = 'NM'
th1 = th %>% filter(Experiment==ex) %>%
    mutate(Genotype = factor(Genotype, levels=gts25)) %>%
    mutate(Treatment = factor(Treatment, levels=c("Control",'Cold'))) %>%
    mutate(cond = sprintf("%s_%dh", Treatment, Timepoint)) %>%
    mutate(cond = factor(cond, levels=c("Control_1h",'Control_25h','Cold_1h','Cold_25h'))) %>%
    arrange(Genotype,Treatment,Timepoint) %>%
    group_by(Genotype, Treatment, Timepoint) %>%
    mutate(lab = str_c(SampleID, cond, Genotype, sep=' ')) %>%
    mutate(clab = ifelse(cond=='Control_1h',as.character(Genotype),''))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.25)
ggsave(sprintf("%s/11.hclust.%s.p.pdf",dirw,ex), p1, width=6, height=10)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.25)
ggsave(sprintf("%s/11.hclust.%s.s.pdf",dirw,ex), p1, width=6, height=10)

p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=4,iter=1200, seed=12,
    var.shape='cond',var.col='cond',var.lab='Genotype', #var.ellipse='Genotype',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0:1,15:16), pal.col='aaas')
ggsave(sprintf("%s/12.tsne.%s.pdf", dirw, ex), p2, width=6, height=6)
#}}}
#}}}

#{{{ fix samples
# TC19 and TC20 genotype
th2 = res$th %>%
    mutate(Genotype=ifelse(SampleID=='TC19', 'Mo17', Genotype)) %>%
    mutate(Genotype=ifelse(SampleID=='TC20', 'B73', Genotype))
th2 %>% filter(SampleID %in% c("TC19",'TC20')) %>% print(width=Inf)

# fix HY samples
sids_y = c('NM76','NM100','TC64','TC66')
sids_n = c('HY93','HY91')
th3b = th2 %>% filter(SampleID %in% sids_y) %>% mutate(Experiment="HY")
th3 = th2 %>% filter(!SampleID %in% sids_n) %>%
    bind_rows(th3b)

# remove TC - time8 samples
sids8 = th3 %>% filter(Experiment=='TC', Timepoint == 8) %>% pull(SampleID)
sids8
th4 = th3 %>% filter(!SampleID %in% sids8)

th4 %>% count(Experiment)
#th4 = complete_sample_list(th4)
fh = file.path(dirw, '01.meta.tsv')
write_tsv(th4, fh, na='')
#}}}

#{{{ read in again
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th
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
    mutate(lab = str_c(SampleID, grp, sep=' '))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1, width=6, height=8)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    expand.x=.3)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1, width=6, height=8)

p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=4,iter=1000, seed=12,
    var.shape='Treatment',var.col='Genotype',var.lab='has',#var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.pdf", dirw, ex), p2, width=6, height=6)

th2 = th1 %>% filter(Genotype=='B73')
tm2 = tm1 %>% filter(SampleID %in% th2$SampleID)
p2 = plot_tsne(tm2,th2,pct.exp=.7,perp=3,iter=800, seed=12,
    var.shape='Treatment',var.col='Treatment',var.lab='has',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.B73.pdf", dirw, ex), p2, width=6, height=6)

th2 = th1 %>% filter(Genotype=='Mo17')
tm2 = tm1 %>% filter(SampleID %in% th2$SampleID)
p2 = plot_tsne(tm2,th2,pct.exp=.7,perp=3,iter=800, seed=12,
    var.shape='Treatment',var.col='Treatment',var.lab='has',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.Mo17.pdf", dirw, ex), p2, width=6, height=6)

th2 = th1 %>% filter(Genotype=='W22')
tm2 = tm1 %>% filter(SampleID %in% th2$SampleID)
p2 = plot_tsne(tm2,th2,pct.exp=.8,perp=4,iter=800, seed=12,
    var.shape='Treatment',var.col='Treatment',var.lab='has',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0,1,2), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.W22.pdf", dirw, ex), p2, width=6, height=6)
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

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1, width=8, height=10)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.2)
ggsave(sprintf("%s/21.hclust.%s.s.pdf",dirw,ex), p1, width=8, height=10)

p2 = plot_tsne(tm1,th1,pct.exp=.6,perp=6,iter=800, seed=12,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0:3,5,6), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.pdf", dirw, ex), p2, width=6, height=6)

th2 = th1 %>% filter(Genotype %in% gts3)
p2 = plot_tsne(tm1,th2,pct.exp=.7,perp=4,iter=1000, seed=12,
    var.shape='Genotype',var.col='Genotype',var.lab='clab',var.ellipse='grp',
    legend.pos='top.center.out', legend.dir='h', legend.box='h',legend.title=F,
    shapes=c(0:3,5,6), pal.col='aaas')
ggsave(sprintf("%s/22.tsne.%s.inbreds.pdf", dirw, ex), p2, width=6, height=6)
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
    mutate(lab = str_c(SampleID, cond, Genotype, sep=' '))
tm1 = tm %>% filter(SampleID %in% th1$SampleID) %>%
    mutate(value=asinh(CPM))

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='pearson',var.col='Treatment',
    expand.x=.25)
ggsave(sprintf("%s/21.hclust.%s.p.pdf",dirw,ex), p1, width=6, height=10)

p1 = plot_hclust(tm1,th1,pct.exp=.7,cor.opt='spearman',var.col='Treatment',
    expand.x=.25)
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



