source("functions.R")

yid = 'rn13a'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) dir.create(dirw)

#{{{ read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Replicate==1, Genotype, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    pal.col='viridis_d', expand.x=.1)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=8, height=12)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    pal.col='viridis_d', expand.x=.1)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=8, height=12)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='', var.col='', var.lab='clab', var.ellipse='Genotype',
    legend.pos='top.left', legend.dir='v', pal.col='')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=10,iter=1500, seed=42,
    var.shape='', var.col='', var.lab='clab', var.ellipse='Genotype',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), width=8, height=8)
#}}}

#{{{ fix
th2 = res$th
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

#{{{ read in again
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = th %>%
    mutate(lab = str_c(Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Replicate==1, Genotype, ''))
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Genotype',
    pal.col='viridis_d', expand.x=.1)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=8, height=12)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Genotype',
    pal.col='viridis_d', expand.x=.1)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=8, height=12)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=T, pca.scale=F,
    var.shape='', var.col='', var.lab='clab', var.ellipse='Genotype',
    legend.pos='top.left', legend.dir='v', pal.col='')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=10,iter=1500, seed=42,
    var.shape='', var.col='', var.lab='clab', var.ellipse='Genotype',
    legend.pos='top.right', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), width=8, height=8)
#}}}

gts = c("B73","Mo17")
gts1 = c("M0021","M0317","M0014","M0016")
gts2 = c('M0275','M0352','M0326','M0323','M0035','M0012','M0017','M0001')
sids_red = th %>% filter(Genotype %in% c(gts, gts1, gts2)) %>% pull(SampleID)
#{{{ RIL haplotype blocks
cp = res$ril$cp
p = plot_ril_genotype(cp, th, sids_red=sids_red)
fo = file.path(dirw, '41.ril.genotype.pdf')
ggsave(fo, p, width=8, height=15)
#}}}

#{{{ write RIL haplotype block coordinates
fw = '~/projects/genome/data/Zmays_B73/15_intervals/20.win11.tsv'
tw = read_tsv(fw)
offs = c(0, cumsum(tw$size)[-nrow(tw)]) + 0:10 * 10e6
tx = tw %>% mutate(off = offs) %>%
    mutate(gstart=start+off, gend=end+off, gpos=(gstart+gend)/2) %>%
    select(rid,chrom,gstart,gend,gpos,off) %>% filter(chrom!='B99')
#
tps = th %>% select(sid=SampleID,Genotype,Replicate) %>% arrange(Genotype) %>%
    mutate(y = 1:n())
tp = res$ril$cp %>% inner_join(tx, by='rid') %>%
    mutate(gstart=start+off, gend=end+off) %>%
    inner_join(tps, by='sid')
#
to = tp %>%
    mutate(sid=str_c(Genotype, Replicate, sep='_')) %>%
    select(sid, chrom,start,end, gt)

fo = file.path(dirw, '42.ril.genotype.tsv')
write_tsv(to, fo)
#}}}


#{{{ convert RIL genotype to phylip
fi = file.path(dird, 'raw', yid, 'ril.rds')
res = readRDS(fi)
write_phylip <- function(tp, fo) {
    #{{{
    nind = nrow(tp); nsite = nchar(tp$gt[1])
    write_tsv(tibble(x=nind, y=nsite), fo, append=T, col_names=F)
    write_tsv(tp, fo, append=T, col_names=F)
    #}}}
}

tp = res$bin %>% filter(rid != 'r11') %>%
    mutate(gt2=ifelse(gt=='a','0',ifelse(gt=='b','2','1'))) %>%
    arrange(SampleID, rid, start, end) %>%
    group_by(SampleID) %>%
    summarise(gt = str_c(gt2, collapse='')) %>%
    ungroup()

fo = file.path(dirw, '41.phy')
write_phylip(tp, fo)
#}}}

#{{{ plot tree
require(ape)
require(treeio)
require(ggtree)
require(tidytree)
fi = file.path(dirw, '41.phy.treefile')
tree = read.newick(fi)
tp = fortify(tree) %>% filter(isTip) %>%
    arrange(y) %>%
    left_join(th, by=c('label'='SampleID')) %>%
    #filter(!is.na(label)) %>%
    mutate(gt = ifelse(Genotype %in% gts, 'ref', ifelse(Genotype %in% gts1, 'select1', ifelse(Genotype %in% gts2,'select2','all'))))

g1 = tp %>% filter(row_number() %% 2 != 0) %>% pull(Genotype)
g2 = tp %>% filter(row_number() %% 2 == 0) %>% pull(Genotype)
identical(g1, g2)
sids = tp %>% filter(row_number() %% 2 != 0) %>% pull(label)

cols10 = c('black', pal_aaas()(10)[2:10])
p1 = ggtree(tree) %<+% tp +
    geom_tiplab(aes(label=lab,col=gt), size=2, hjust=0, align=T, linesize=.5) +
    scale_color_manual(values = cols10) +
    scale_x_continuous(expand=expand_scale(mult=c(.02,.1))) +
    scale_y_continuous(expand=expand_scale(mult=c(.01,.01))) +
    guides(color=F)
fo = file.path(dirw, '43.phylogeny.pdf')
ggsave(p1, file=fo, width=6, height=15)
#}}}


