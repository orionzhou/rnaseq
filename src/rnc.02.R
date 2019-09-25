source("functions.R")
size.gene = read_genome_conf()$gene %>% select(gid, size=size.exon)

yid = 'rnc02'
yids = c('rn13a','rn15c','rn15d','rn18c')
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
th = th %>% mutate(lab = str_c(Treatment,Tissue,sep="_"))

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .7) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=10,
              pca = T, max_iter = 1500, check_duplicates = F)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp) +
    #geom_text_repel(aes(x=V1,y=V2,label=Tissue), size=2) +
    geom_point(aes(x=V1, y=V2, color=Treatment, shape=Treatment), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_color_aaas(name = 'Study') +
    scale_shape(name = 'Study', solid=F) +
    otheme(legend.pos='top.left', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines'))
fp = file.path(dirw, "25.tsne.pdf")
ggsave(p_tsne, filename = fp, width=9, height=9)
#}}}

