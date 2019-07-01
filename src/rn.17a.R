source("functions.R")

yid = 'rn17a'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in, filter/fix samples
ref = t_cfg %>% filter(yid == !!yid) %>% pull(ref)
th = rnaseq_sample_meta(yid)
tt = rnaseq_mapping_stat(yid)
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
sum_stat_tibble(tt)

sids_keep = tt %>% filter(mapped>5) %>% pull(SampleID)
sum_stat_tibble(tt %>% filter(SampleID %in% sids_keep))

# fix th
th2 = th %>% filter(SampleID %in% sids_keep) %>%
    mutate(Tissue = ifelse(SampleID == 'SRR445601', 'tassel', Tissue)) %>%
    mutate(Tissue = ifelse(SampleID == 'SRR445416', 'tassel', Tissue)) %>%
    mutate(Genotype = ifelse(SampleID == 'SRR426798', 'Mo17', Genotype)) %>%
    mutate(Genotype = ifelse(SampleID == 'SRR426814', 'M37W', Genotype))
th2 = complete_sample_list(th2)

th = th2
tt = tt %>% filter(SampleID %in% th$SampleID)

fh = file.path(dirw, 'meta.tsv')
write_tsv(th, fh, na='')
# run snakemake again
#}}}

res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
th = th %>% mutate(lab = sprintf("%s_%s", Tissue, Genotype))

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
    select(taxa, everything())
p1 = ggtree(tree, layout = 'rectangular') +
    scale_x_continuous(expand = expand_scale(0,2)) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+%
    tp + geom_tiplab(aes(label=lab, color=Tissue), size=2.5) +
    scale_color_aaas()
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width=8, height=18)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=5,
              pca = T, max_iter = 1500)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp) +
    geom_text_repel(aes(x=V1,y=V2,label=Genotype), size=2.5) +
    geom_point(aes(x=V1, y=V2, color=Tissue, shape=Tissue), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(0:6)) +
    scale_color_aaas() +
    otheme(legend.pos='top.left', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines'))
fp = file.path(dirw, "25.tsne.pdf")
ggsave(p_tsne, filename = fp, width=8, height=8)
#}}}


