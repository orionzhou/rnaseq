source("functions.R")

yid = 'rn99d'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ reformat meta-table
fi = file.path(dirw, '00.meta.tsv')
ti = read_tsv(fi)
hr_map = c("11am"=11,'12pm'=12,'2pm'=14,'4pm'=16)
th = ti %>%
    separate(treatment,c('rep','day','hr','Treatment'), sep='_') %>%
    mutate(day = as.integer(str_sub(day,4))) %>%
    mutate(hr = hr_map[hr]) %>%
    mutate(hac = (day-9)*24 + hr - 11) %>%
    select(SampleID=sid,Tissue=tissue,Genotype=genotype,Treatment,hac) %>%
    mutate(gt_treat=sprintf("%s_%s",Genotype,Treatment)) %>%
    mutate(tis_txt=ifelse(Tissue=='leaf_2', '', Tissue))
#}}}

#{{{ read in, filter/fix samples
ref = t_cfg %>% filter(yid == !!yid) %>% pull(ref)
#th = rnaseq_sample_meta(yid)
tt = rnaseq_mapping_stat(yid)
res = rnaseq_cpm_raw(yid)
tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
sum_stat_tibble(tt)

sids_keep = tt %>% filter(mapped>1.5) %>% pull(SampleID)
sum_stat_tibble(tt %>% filter(SampleID %in% sids_keep))

# fix th
th2 = th %>% filter(SampleID %in% sids_keep)

th = th2
tt = tt %>% filter(SampleID %in% th$SampleID)

fh = file.path(dirw, 'meta.tsv')
write_tsv(th, fh, na='')
# run snakemake again
#}}}
tm = tm %>% filter(SampleID %in% th$SampleID)

res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
th = th %>% mutate(lab = sprintf("%s", Genotype)) %>%
    replace_na(list(tis_txt=''))

#{{{ hclust
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
e = tw %>% filter(gid %in% gids) %>% select(-gid)
dim(e)

cor_opt = "spearman"
cor_opt = "pearson"
hc_opt = "ward.D"
hc_title = sprintf("dist: %s\nhclust: %s", cor_opt, hc_opt)
edist <- as.dist(1-cor(e, method = cor_opt))
ehc <- hclust(edist, method = hc_opt)
tree = as.phylo(ehc)
lnames = ehc$labels[ehc$order]
#
tp = th %>% mutate(taxa = SampleID) %>%
    mutate(lab = sprintf("%s_%s_%d", gt_treat, tis_txt, hac)) %>%
    select(taxa, everything())
p1 = ggtree(tree, layout = 'rectangular') +
    scale_x_continuous(expand = expand_scale(0,.8)) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+%
    tp + geom_tiplab(aes(label=lab, color=gt_treat), size=2.5) +
    scale_color_aaas()
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width=8, height=10)
#}}}


th = th %>%
    #filter(Tissue=='leaf_2', hac %in% c(5,29,48,53,77,149)) %>%
    filter(Tissue=='leaf_2', Genotype=='B73') %>%
    mutate(hac=sprintf("hac%03d",hac)) %>%
    mutate(grp = sprintf("%s_%s_%s", hac, Genotype, Treatment))

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .05) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=10,
              pca = T, max_iter = 1500)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp) +
    geom_text_repel(aes(x=V1,y=V2,label=hac), size=2.5) +
    geom_point(aes(x=V1, y=V2, color=hac, shape=gt_treat), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(name='Genotype',values = c(15,0,16,1)) +
    scale_color_viridis_d(name = 'hours after cold',direction=-1) +
    #scale_color_ucscgb(name = 'hours after cold') +
    otheme(legend.pos='bottom.left', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(color = F)
fp = file.path(dirw, "25.tsne.1.pdf")
ggsave(p_tsne, filename = fp, width=6, height=6)
#}}}


#{{{ DESeq2
require(DESeq2)
require(edgeR)
th1 = th
tm1 = tm %>% filter(SampleID %in% th1$SampleID)
#{{{ prepare data
vh = th1 %>% mutate(Genotype = factor(Genotype)) %>%
    mutate(Treatment=factor(Treatment)) %>%
    mutate(hac = factor(hac)) %>% arrange(SampleID)
vh.d = column_to_rownames(as.data.frame(vh), var = 'SampleID')
gids = tm1 %>% group_by(gid) %>% summarise(n.sam = sum(ReadCount >= 10)) %>%
    filter(n.sam > .05 * nrow(vh)) %>% pull(gid)
vm = tm1 %>% filter(gid %in% gids) %>%
    select(SampleID, gid, ReadCount)
x = readcount_norm(vm)
mean.lib.size = mean(x$tl$libSize)
vm = x$tm
vm.w = vm %>% select(SampleID, gid, ReadCount) %>% spread(SampleID, ReadCount)
vm.d = column_to_rownames(as.data.frame(vm.w), var = 'gid')
stopifnot(identical(rownames(vh.d), colnames(vm.d)))
#}}}
# DESeq2
dds = DESeqDataSetFromMatrix(countData=vm.d, colData=vh.d, design= ~ grp)
dds = estimateSizeFactors(dds)
dds = estimateDispersions(dds, fitType = 'parametric')
disp = dispersions(dds)
#dds = nbinomLRT(dds, reduced = ~ 1)
dds = nbinomWaldTest(dds)
resultsNames(dds)

call_deg <- function(t1, t2, dds, gids) {
    #{{{
    res1 = results(dds, contrast=c("grp",t1,t2), pAdjustMethod="fdr")
    stopifnot(rownames(res1) == gids)
    tibble(gid = gids, padj = res1$padj, lfc = res1$log2FoldChange) %>%
        replace_na(list(padj = 1))
    #}}}
}

t_ds = th1 %>% distinct(hac) %>% arrange(hac)
t_ds = crossing(t_ds, gt=c('B73','Mo17')) %>%
    mutate(t1 = sprintf("%s_%s_control", hac, gt)) %>%
    mutate(t2 = sprintf("%s_%s_cold", hac, gt)) %>%
    mutate(res = map2(t1,t2, call_deg, dds=dds, gids=gids)) %>%
    select(-t1,-t2) %>% unnest()

fo = file.path(dirw, '31.deg.rds')
saveRDS(t_ds, file=fo)

to = t_ds %>% filter(padj < .01) %>% dplyr::count(hac, gt) %>%
    spread(hac, n)
fo = file.path(dirw, '32.deg.tsv')
write_tsv(to, fo)
#}}}

#{{{ set overlap
require(UpSetR)
fi = file.path(dirw, '31.deg.rds')
t_ds = readRDS(fi)
hacs = t_ds %>% distinct(hac) %>% pull(hac)

gt = 'Mo17'
ts1 = t_ds %>% filter(gt==!!gt, padj < .01) %>%
    select(hac,gid)
x = list()
for (hac in hacs) {
    x[[hac]] = ts1 %>% filter(hac==!!hac) %>% pull(gid)
}

fo = sprintf('%s/42.%s.pdf', dirw, gt)
pdf(fo, width=8, height=6, onefile=F)
upset(fromList(x), sets=hacs,
      number.angles = 0, point.size = 2, line.size = 1,
      mainbar.y.label = "Genre Intersections",
      sets.x.label = "Num. DEGs at each time point",
      text.scale = c(1.1, 1.1, 1, 1, 1.3, 1),
      order.by='degree', keep.order=T)
dev.off()

require(ggalluvial)

gt = 'Mo17'
degs = c("up-DE",'down-DE','non-DE')
tp = t_ds %>% filter(gt == !!gt) %>%
    mutate(deg = ifelse(padj>.01, 'non-DE', ifelse(lfc<0, 'up-DE', 'down-DE'))) %>%
    select(gid,hac,deg) %>% spread(hac, deg) %>%
	count(hac005,hac029,hac048,hac053,hac077,hac149) %>%
	mutate(subject=1:length(n)) %>%
	gather(hac, deg, -subject, -n) %>%
	mutate(deg = factor(deg, levels=degs))
#
cols3 = c(pal_aaas()(2), 'gray')
p = ggplot(tp, aes(x=hac, stratum=deg, alluvium=subject, y=n, fill=deg, label=deg)) +
    scale_fill_manual(values=cols3) +
    #geom_flow(stat='alluvium', lode.guidance='rightleft', color='darkgray') +
	geom_flow() +
	geom_text(stat = 'stratum', size=2) +
    geom_stratum(alpha = .6) +
	otheme(xtext=T,ytext=T,xtick=T,ytick=T)
fo = sprintf('%s/43.%s.pdf', dirw, gt)
ggsave(p, file=fo, width=8, height=6)
#}}}



