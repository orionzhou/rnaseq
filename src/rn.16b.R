source("functions.R")

yid = 'rn16b'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>% replace_na(list(Treatment='')) %>%
    mutate(lab=str_c(Tissue,Treatment,Replicate,sep='_')) %>%
    mutate(grp=str_c(Tissue,Treatment,sep='_')) %>%
    mutate(grp = str_replace(grp, '_$', '')) %>%
    mutate(clab = ifelse(Replicate==1, grp, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    pal.col='simpsons', expand.x=.3)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=8, height=25)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    pal.col='simpsons', expand.x=.3)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=8, height=25)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='', var.col='Tissue', var.lab='clab', var.ellipse='grp', leg.col=F,
    legend.pos='top.left', legend.dir='v', pal.col='simpsons')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=10,iter=1500, seed=2,
    var.shape='', var.col='Tissue', var.lab='clab', var.ellipse='grp', leg.col=F,
    legend.pos='top.left', legend.dir='v', pal.col='simpsons')
ggsave(file.path(dirw, '11.tsne.pdf'), p3, width=8, height=8)
#}}}

#{{{ filter/fix
sids_keep = res$bamstat %>% filter(pair_map+unpair_map>2e6) %>% pull(sid)
sids_rm = c("SRR1620908","SRR1620913",'SRR1620838','SRR1620927')
#
th2 = res$th %>% filter(SampleID %in% sids_keep) %>%
    filter(!SampleID %in% sids_rm) %>%
    mutate(Treatment=ifelse(str_detect(Treatment,"^3[\\.\\-]7\\.1"), 'TC_control', Treatment))
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

#{{{ read in
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>% replace_na(list(Treatment='')) %>%
    mutate(lab=str_c(Tissue,Treatment,Replicate,sep='_')) %>%
    mutate(grp=str_c(Tissue,Treatment,sep='_')) %>%
    mutate(grp = str_replace(grp, '_$', '')) %>%
    mutate(clab = ifelse(Replicate==1, grp, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    pal.col='simpsons', expand.x=.3)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=8, height=25)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    pal.col='simpsons', expand.x=.3)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=8, height=25)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='', var.col='Tissue', var.lab='clab', var.ellipse='grp', leg.col=F,
    legend.pos='top.left', legend.dir='v', pal.col='simpsons')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=10,iter=1500, seed=2,
    var.shape='', var.col='Tissue', var.lab='clab', var.ellipse='grp', leg.col=F,
    legend.pos='top.left', legend.dir='v', pal.col='simpsons')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=8, height=8)
#}}}

