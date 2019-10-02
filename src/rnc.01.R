source("functions.R")
gcfg = read_genome_conf()
tsyn = read_syn(gcfg)
size.gene = gcfg %>% select(gid, size=size.exon)
yid = 'rnc01'
yids = c('rn10a','rn11a','rn13b','rn14b','rn14c','rn14e',"rn16b","rn16c","rn18g")
t_cfg %>% filter(yid %in% yids)
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ merge datasets
r = tibble(yid = yids) %>% mutate(data = map(yid, read_rnaseq))
th = r %>% mutate(d2 = map(data, 'th')) %>% select(d2) %>% unnest()
t_rc = r %>% mutate(d2 = map(data, 'tm')) %>% select(d2) %>% unnest() %>%
    select(gid, SampleID, ReadCount)
dim(th); dim(t_rc)
th %>% dplyr::count(Tissue) %>% print(n=40)
th %>% dplyr::count(Tissue,Genotype,Treatment) %>% print(n=40)
th %>% distinct(study,Tissue,Genotype,Treatment) %>% dplyr::count(study)
#
res = readcount_norm(t_rc, size.gene)
#
ths = th %>% distinct(Tissue, Genotype, Treatment) %>%
    mutate(nSampleID = sprintf("%s_%d", yid, 1:length(Tissue)))
t_map = th %>% inner_join(ths, by = c("Tissue", "Genotype", "Treatment")) %>%
    select(SampleID, nSampleID)
th_m = ths %>% select(SampleID=nSampleID, Tissue, Genotype, Treatment)
#
t_rc_m = t_rc %>% inner_join(t_map, by = 'SampleID') %>%
    mutate(SampleID = nSampleID) %>%
    group_by(gid, SampleID) %>%
    summarise(ReadCount = sum(ReadCount)) %>%
    ungroup()
resm = readcount_norm(t_rc_m, size.gene)
#
res$th = th
res$th_m = th_m
res$tm_m = resm$tm
#
fo = file.path(dirw, 'cpm.rds')
saveRDS(res, file = fo)
#}}}

res = rnaseq_cpm(yid)
th = res$th_m; tm = res$tm_m
th = th %>% rename(yid = Treatment) %>%
    inner_join(t_cfg, by='yid') %>%
    mutate(lgd = ifelse(author=='zhou2018', 'Zhou2018 B&M atlas [23]', lgd)) %>%
    mutate(lab = str_c(author,Tissue,sep="_"))

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
    scale_x_continuous(expand = expand_scale(0,2.5)) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+%
    tp + geom_tiplab(aes(label=lab, color=lgd), size=2.5) +
    scale_color_aaas()
fo = sprintf("%s/21.cpm.hclust.pdf", dirw)
ggsave(p1, filename = fo, width=8, height=20)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=9,
              pca = T, max_iter = 2000)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp) +
    geom_text_repel(aes(x=V1,y=V2,label=Tissue), size=2) +
    geom_point(aes(x=V1, y=V2, color=lgd), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_color_npg(name = 'Study') +
    otheme(legend.pos='top.center.out', legend.title=F,
           xtitle=T, ytitle=T,
           margin = c(2,.2,.2,.2)) +
    guides(col=guide_legend(nrow = 3)) +
    theme(axis.ticks.length = unit(0, 'lines'))
fp = file.path(dirw, "25.tsne.pdf")
ggsave(p_tsne, filename = fp, width=9, height=9)
#}}}


#{{{ #genes expressed in 0-23 tissues
tsh_e = tm %>% filter(gid %in% gcfg$gene$gid[gcfg$gene$ttype=='mRNA']) %>%
    mutate(silent = CPM < 1) %>%
    group_by(gid, silent) %>% summarise(n.tis = n()) %>% ungroup()
tsh_es = tsh_e %>%
    group_by(gid) %>% summarise(n.tis.tot = sum(n.tis)) %>% ungroup()
tsh_es %>% count(n.tis.tot)
etags = c('Silent', 'Tissue specific', 'Intermediate frequency', 'Constitutive')
tsh_e = tsh_e %>% filter(!silent) %>%
    right_join(tsh_es, by = 'gid') %>%
    replace_na(list(n.tis = 0)) %>%
    mutate(prop.tis = n.tis / n.tis.tot,
           etag = ifelse(prop.tis == 0, etags[1],
                  ifelse(prop.tis <= 0.2, etags[2],
                  ifelse(prop.tis < 0.8, etags[3], etags[4])))) %>%
    mutate(etag = factor(etag, levels = etags)) %>%
    select(gid, n.tis, prop.tis, etag)
tsh_e %>% count(etag)

tp = tsh_e %>% count(n.tis, etag) %>% rename(num_genes = n)
cat("genes expressed in >=1 tissues:\n")
sum(tp %>% filter(n.tis > 0) %>% pull(num_genes))
cat("prop. genes silent, constitutive. etc:\n")
tp %>% group_by(etag) %>% summarise(n = sum(num_genes)) %>% mutate(p=n/sum(n))
p = ggplot(tp) +
    geom_bar(aes(x = n.tis, y = num_genes, fill = etag), stat = 'identity', width = .8) +
    scale_x_continuous(name = 'Number Tissues with Expression', expand=expand_scale(mult=c(.03,.03))) +
    scale_y_continuous(name = "Number Genes", expand=expand_scale(mult=c(0,.03))) +
    scale_fill_npg() +
    otheme(xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T,
           margin = c(.3,1.3,.3,.3), "lines", legend.title=F) +
    theme(legend.pos = c(.5,1), legend.justification = c(.5,1))
#
fo = file.path(dirw, '31.tis.expression.pdf')
ggsave(fo, p, width=6, height=6)
#}}}

#{{{ synteny proportion
tp  = tsh_e %>% inner_join(tsyn, by = 'gid') %>% group_by(ftype,etag) %>%
    summarise(n = n()) %>%
    mutate(ntot = sum(n), prop = n/ntot) %>%
    ungroup() %>%
    select(ftype, ntot, etag, n, prop)

#}}}


