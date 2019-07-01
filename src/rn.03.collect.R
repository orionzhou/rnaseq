source("functions.R")
require(ape)
require(ggtree)
t_cfg %>% print(n=40)

yid = 'rn17d'
#{{{ config
Sys.setenv(R_CONFIG_ACTIVE = yid)
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))
#
cfg = t_cfg %>% filter(yid == !!yid)
stopifnot(nrow(cfg) == 1)
study = cfg %>% pull(study)
genome = cfg %>% pull(ref)
gcfg = read_genome_conf(genome)
#}}}

#{{{ [meta=F] mapping stats
th = rnaseq_sample_meta(yid)
tt = read_mapping_stat(yid)
sum_stat_tibble(tt)

fo = file.path(dirw, '10.mapping.stat.tsv')
write_tsv(tt, fo)
#}}}

#{{{ [meta=T] merge datasets
yids = str_split(t_cfg$study[t_cfg$yid==yid], "[\\+]")[[1]]
yids
#
th = tibble(); t_rc = tibble()
for (yid1 in yids) {
    diri = file.path(dird, 'raw', yid1)
    th1 = rnaseq_sample_meta(yid1)
    if(yid1 == 'rn12a') {
        th1 = th1 %>% filter(Treatment == 'WT')
    } else if(yid1 == 'rn13b') {
        th1 = th1 #%>% filter(!str_detect(Treatment, "ET"))
    } else if(yid1 == 'rn17c') {
        th1 = th1 %>% filter(Treatment == 'con')
    } else if(yid1 == 'rn18g') {
        th1 = th1 %>% filter(Genotype == 'B73')
    }
    oyids = th1$SampleID
    nyids = sprintf("%s_%s", yid1, oyids)
    names(oyids) = nyids
    th1 = th1 %>% mutate(SampleID = nyids) %>%
        replace_na(list(Treatment='?')) %>%
        mutate(Treatment=sprintf("%s|%s", yid1, Treatment)) %>%
        select(SampleID, Tissue, Genotype, Treatment, everything())
    th = rbind(th, th1)
    fi = file.path(diri, 'featurecounts.tsv')
    t_rc1 = read_tsv(fi) %>% select(one_of(c('gid', oyids))) %>%
        rename(!!oyids)
    stopifnot(ncol(t_rc1) == nrow(th1) + 1)
    if(nrow(t_rc) == 0) {
        t_rc = t_rc1
    } else {
        t_rc = t_rc %>% inner_join(t_rc1, by = 'gid')
    }
}
dim(th); dim(t_rc)
th %>% count(Tissue) %>% print(n=40)
th %>% count(Tissue,Genotype,Treatment) %>% print(n=40)
th %>% count(Tissue,Genotype,Treatment) %>% count(n) %>% print(n=40)

fo = sprintf("%s/05_read_list/%s.tsv", dird, yid)
write_tsv(th, fo)
diro = file.path(dird, '08_raw_output', yid)
system(sprintf("mkdir -p %s", diro))
map_int(sprintf("touch %s/%s", diro, c("multiqc.html",'trimming.tsv','bamstats.tsv')), system)
fo = sprintf("%s/featurecounts.tsv", diro)
write_tsv(t_rc, fo)
# run Snakemake rc2cpm
#}}}

th = rnaseq_sample_meta(yid)
fi = file.path(dird, 'raw', yid, 'cpm.rds')
res = readRDS(fi)
tl = res$tl; tm = res$tm
tiss = unique(th$Tissue); genos = unique(th$Genotype);
treas = unique(th$Treatment); reps = unique(th$Replicate)

#{{{ [optional] fix/remove mis-labelled samples, save to mexx.c.tsv
#{{{ cut hc tree
#cls = cutree(ehc, h = .1)
#tcl = tibble(SampleID = names(cls), grp = as.integer(cls))
#th2 = th %>% inner_join(tcl, by = 'SampleID')
#th3 = th2 %>% group_by(Tissue, Genotype) %>%
    #summarise(nrep = length(Replicate), ngrp = length(unique(grp))) %>%
    #ungroup() %>%
    #filter(nrep > 1, ngrp > 1)
#th2 %>% inner_join(th3, by = c("Tissue", "Genotype")) %>% print(n=40)
#}}}
ft = file.path(dirw, '10.mapping.stat.tsv')
tt = read_tsv(ft)
fh1 = sprintf("%s/05_read_list/%s.tsv", dird, yid)
th = read_tsv(fh1)
fh2 = sprintf("%s/05_read_list/%s.c.tsv", dird, yid)
if(yid == 'me13c') {
    #{{{
    th = th %>%
        mutate(Genotype = ifelse(SampleID=='SRR767691','Oh43',Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='SRR651079','Oh7b',Genotype))
    #}}}
} else if(yid == 'me14c') {
    th = th %>% filter(! SampleID %in% c("SRR254169"))
} else if(yid == 'me14d') {
    th = th %>% filter(! SampleID %in% c("SRR1573518", 'SRR1573513'))
} else if(yid == 'me16b') {
    #{{{ old
    th = th %>% filter(!SampleID %in% c("SRR1620930","SRR1620929","SRR1620927",
                                        "SRR1620908","SRR1620913"))
    #}}}
} else if(yid == 'me17a') {
    #{{{
    th = th %>%
        mutate(Tissue = ifelse(SampleID == 'SRR445601', 'tassel', Tissue)) %>%
        mutate(Tissue = ifelse(SampleID == 'SRR445416', 'tassel', Tissue)) %>%
        mutate(Genotype = ifelse(SampleID == 'SRR426798', 'Mo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID == 'SRR426814', 'M37W', Genotype))
    #}}}
} else if(yid == 'me18b') {
    #{{{
    th = th %>%
        filter(!SampleID %in% c("SRR5786263", "SRR5786217")) %>%
        mutate(Treatment=ifelse(SampleID=='SRR5786370','II',Treatment)) %>%
        mutate(Treatment=ifelse(Genotype=='B73' & Treatment=='I' &
                                Replicate %in% 3:8, 'Z1', Treatment)) %>%
        mutate(Treatment=ifelse(Genotype=='B73' & Treatment=='II' &
                                Replicate %in% c(3,4,10:12), 'Z2', Treatment)) %>%
        mutate(Treatment=ifelse(Genotype=='B73' & Treatment=='III' &
                                Replicate %in% c(4,7:10), 'Z3', Treatment))
    #}}}
} else if(yid == 'me99a') {
    #{{{
    samples_low = tt %>% filter(passed < 1e7) %>% pull(SampleID)
    samples2 = c("SRR8043600","SRR8043169","SRR8043151","SRR8043196","SRR8043556")
    th = th %>%
        filter(!SampleID %in% c(samples_low,samples2)) %>%
        mutate(Tissue=ifelse(SampleID=='SRR8043188','leaf',Tissue))
    i = which(th$SampleID == 'SRR5691477')
    th$Genotype[i] = 'PHN11x?'; th$inbred[i] = F
    #}}}
} else if(yid == 'me99b') {
    #{{{
    gts = c("B73", "Mo17", "B73xMo17")
    tissues = sort(unique(th$Tissue))
    th = th %>%
        #filter(Genotype %in% gts) %>%
        filter(! SampleID %in% c('BR207', 'BR230', "BR235")) %>%
        mutate(Genotype = ifelse(SampleID=='BR003', 'Mo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR004', 'B73', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR006', 'B73xMo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR007', 'Mo17', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR029', 'B73', Genotype)) %>%
        mutate(Genotype = ifelse(SampleID=='BR032', 'Mo17', Genotype))
    th = th %>% mutate(Replicate = '')
    th = sra_fill_replicate(th)
    #}}}
} else if(yid == 'me99c') {
    #{{{
    th = th %>%
        mutate(Tissue = ifelse(SampleID == 'bm318', 'Leaf', Tissue))
    #}}}
}
write_tsv(th, fh2, na = '')
# re-normalize everything [if no sample is removed then no need to do this]
#}}}

#{{{ hclust
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
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
tp = th %>% mutate(taxa = SampleID, lab = SampleID)
if(length(tiss)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Tissue), lab)
if(length(genos)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Genotype), lab)
if(length(treas)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Treatment), lab)
if(length(reps)>1) tp = tp %>% mutate(lab = sprintf("%s %s", lab, Replicate), lab)
tp = tp %>% select(taxa, everything())
cols1 = c('gray80','black','red','seagreen3', pal_d3()(5))
if(length(unique(tp$Tissue)) > 15) tp = tp %>% mutate(Tissue='')
p1 = ggtree(tree, layout = 'rectangular') +
    #geom_tiplab(size = labsize, color = 'black') +
    scale_x_continuous(expand = expand_scale(.001,1)) +
        #mult=c(.02, config::get("hc.x.expand")))) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+% tp + geom_tiplab(
        #aes(label=lab, color = Tissue),
        aes(label=lab, color = Genotype), size=2.5
        #aes(label=lab, color = Treatment),
        #size = config::get("hc.labsize"),
        #offset = config::get("hc.x.off"),
        #family='mono'
        ) +
    #scale_color_npg()
    scale_color_aaas()
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
#ggsave(p1, filename = fo, width=config::get("hc.wd"), height=config::get("hc.ht"))
ggsave(p1, filename = fo, width=8, height=12)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .8) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=10,
              pca = T, max_iter = 1500)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
    #replace_na(list(Treatment='')) %>%
    #mutate(txt = sprintf("%s_%s", Treatment, Replicate))
    #mutate(txt = sprintf("%s_%s", Genotype, Replicate))
if(length(unique(tp$Tissue)) > 15) tp = tp %>% mutate(Tissue='')
p_tsne = ggplot(tp) +
    geom_text_repel(data=tp, aes(x=V1,y=V2,label=Tissue), size=2.5, alpha=.99) +
    #geom_point(aes(x=V1, y=V2, color=Tissue), size=2) +
    #geom_point(aes(x=V1, y=V2, color=Tissue, shape=inbred), size=2) +
    geom_point(aes(x=V1, y=V2, color=hr, shape=Genotype), size=2) +
    #geom_point(aes(x=V1, y=V2, color=Treatment), size=2) +
    #stat_ellipse(aes(x=V1, y=V2, fill=Tissue), linetype=1, alpha=.4) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(15,0,16,1)) +
    #scale_color_simpsons() +
    scale_color_viridis() +
    #scale_color_gradientn(name='hours after day9', colors=pal_gradient()) +
    #scale_color_manual(values = pal_aaas()(10)) +
    otheme(legend.pos = 'top.left', legend.dir = 'v',
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(color = guide_legend(ncol=1, byrow = T))
fp = file.path(dirw, "25.tsne.pdf")
ggsave(p_tsne, filename = fp, width = 8, height = 8)
#}}}

#{{{ # ggtree + heatmap
is_tip <- tree$edge[,2] <= length(tree$tip.label)
ordered_tips = tree$edge[is_tip,2]
lnames = tree$tip.label[ordered_tips]
df_mat = data.frame(as.matrix(1-edist))[,rev(lnames)]
#
labsize = config::get("hc.labsize")
x.expand = config::get("hc.x.expand")
x.off = config::get("hc.x.off")
cols1 = c('gray80','black','red','seagreen3', pal_d3()(5))
p1 = ggtree(tree) +
    #geom_tiplab(size = labsize, color = 'black') +
    scale_x_continuous(expand = expand_scale(mult=c(.02,x.expand))) +
    scale_y_discrete(expand = c(.01,0)) +
    theme_tree2()
p1 = p1 %<+% tp +
    geom_tiplab(aes(label = lab), size = labsize, offset = x.off, family = 'mono')
    #geom_text(aes(label = Rep), size = 2, nudge_x = .022, hjust = 0)
p2 = gheatmap(p1, df_mat, offset = 10, width = 10, colnames = F)
fo = file.path(dirw, '22.heatmap.pdf')
ggsave(p2, filename = fo, width = 30, height = 20)
#}}}
#{{{ # pheatmap
mat = 1 - as.matrix(edist)
mat = mat[lnames, lnames]
t_hm = t_hc %>% mutate(SampleID = factor(SampleID, levels = lnames)) %>%
    arrange(SampleID)
fo = sprintf("%s/22.heatmap.pdf", dirw)
pheatmap(
    mat               = mat,
    #color             = inferno(length(mat_breaks) - 1),
    #breaks            = mat_breaks,
    border_color      = NA,
    cluster_cols      = F,
    cluster_rows      = F,
    show_colnames     = F,
    show_rownames     = T,
    labels_row        = t_hm$lab,
    #annotation_row    = th[,c('SampleID','yid')],
    #annotation_colors = pal_d3()(10),
    drop_levels       = T,
    fontsize          = 6,
    main              = study,
    filename          = fo,
    width             = config::get("hm.wd"),
    height            = config::get("hm.ht")
)

  cellwidth = 30, cellheight = 30, scale = "none",
  treeheight_row = 200,
  kmeans_k = NA,
  show_rownames = T, show_colnames = F,
  clustering_method = "complete",
  cluster_rows = T, cluster_cols = T,
  #clustering_distance_rows = drows1, 
  #clustering_distance_cols = dcols1,
  #annotation_col = ta,
  #annotation_colors = ann_colors,
#}}}
#{{{ # PCA
pca <- prcomp(asinh(e), center = F, scale. = F)
x = pca['rotation'][[1]]
y = summary(pca)$importance
y[,1:5]
xlab = sprintf("PC1 (%.01f%%)", y[2,1]*100)
ylab = sprintf("PC2 (%.01f%%)", y[2,2]*100)
#
tp = as_tibble(x[,1:5]) %>%
    add_column(SampleID = rownames(x)) %>%
    left_join(th, by = 'SampleID') %>%
    mutate(Treatment = factor(Treatment), Replicate = factor(Replicate))
fo = file.path(dirw, '25.pca.pdf')
plot_pca(tp, fo, opt = config::get("pca.opt"), labsize = config::get("pca.labsize"),
         wd = config::get("pca.wd"), ht = config::get("pca.ht"))
#}}}




