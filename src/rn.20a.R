source("functions.R")
require(fastqcr)
dird = '~/projects/rnaseq/data'
dirw = file.path(dird, '11_qc', 'rn20a')
diri = '/home/springer/zhoux379/projects/rnaseq/data/cache/rn20a'

sid = 'b01'

#{{{ read1 UMI
#{{{ read umi_cnt and save
ti = crossing(x=1:12, y=1:8) %>%
    mutate(yl = LETTERS[y]) %>%
    mutate(fn = str_c(yl,x,sep='')) %>%
    rbind(tibble(x=13, y=9, yl='I',fn='undetermined')) %>%
    mutate(fi = sprintf("%s/04_umi_cnt/%s/%s.tsv", diri, sid, fn)) %>%
    mutate(res = map(fi, read_tsv, col_names=c('umi','cnt'))) %>%
    select(x, y, yl, fn, res) %>%
    unnest()

fo = sprintf("%s/00.umi_cnt.%s.rds", dirw, sid)
saveRDS(ti, file=fo)
#}}}

fi = sprintf("%s/00.umi_cnt.%s.rds", dirw, sid)
ti = readRDS(fi)
n_read = sum(ti$cnt); n_umi = nrow(ti); pct_dedup = percent(n_umi/n_read)
tit = sprintf("%s UMIs / %s total reads = %s", number(n_umi,big.mark=','), number(n_read,big.mark=','), pct_dedup)

#{{{ plot reads + umi
tp = ti %>% group_by(x,y,fn) %>%
    summarise(cnt = sum(cnt), umi=n()) %>% ungroup() %>%
    mutate(lab = str_c(number(cnt), number(umi), sep="\n"))
mid = (min(tp$cnt) + max(tp$cnt)) / 2
tp = tp %>% mutate(color = ifelse(cnt < mid, 'white','black'))
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=cnt), color='black') +
    geom_text(aes(x,y, label=lab, color=cnt<mid), hjust=1, size=2.5, nudge_x=.4) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = sprintf("%s/01.nseq.%s.pdf", dirw, sid)
ggsave(fo, p, width=8, height=5)
#}}}

#{{{ fastqc
f_qc = sprintf("%s/02_trim/%s_2_fastqc.zip", diri, sid)
qc = qc_read(f_qc, modules=c("Basic Statistics", "Sequence Duplication Levels"))
n_read2 = as.integer(qc$basic_statistics$Value[[4]])
pct_dedup2 = qc$total_deduplicated_percentage
cat(str_c(number(n_read2), pct_dedup2, '\n', sep=', '))
#}}}

#{{{ umi - distri [obsolete]
tp = ti %>% filter(x<13, cnt<=2000) %>%
    mutate(cnti = cut_interval(cnt, n=10)) %>%
    count(x,yl,cnti)
p = ggplot(tp) +
    geom_bar(aes(cnti, n),stat='identity') +
    scale_x_discrete(name='# UMI occurences', expand=expand_scale(mult=c(.05,.05))) +
    scale_y_continuous(expand=expand_scale(mult=c(.05,.05))) +
    facet_grid(yl ~ x) +
    otheme(xtitle=T, xtext=T, ytext=T, legend.pos='none') +
    theme(axis.text.x=element_text(angle=60,hjust=1,vjust=1)) +
    #theme(panel.border = element_blank()) +
    ggtitle('UMI distribution') +
    theme(plot.title=element_text(hjust=.5))
fo = file.path(dirw, 'b05.umi.distr.pdf')
ggsave(fo, p, width=12, height=8)
#}}}
#}}}

#{{{ read sample meta, DGE, normalize and Save
read_samplelist_stress <- function() {
    #{{{
    dir_st = '/home/springer/zhoux379/projects/stress/data/01_exp_design'
    f1 = file.path(dir_st, "00.set1.tsv")
    f2 = file.path(dir_st, "00.set2.tsv")
    t1 = read_tsv(f1)
    t2 = read_tsv(f2)
    t1 %>% bind_rows(t2) %>% replace_na(list(Rep=1)) %>% rename(SampleID=Code)
    #}}}
}
th = read_samplelist_stress()

read_dge <- function(sid = 'b03') {
    #{{{
    diri = '/home/springer/zhoux379/projects/rnaseq/data/cache/rn20a'
    fi1 = sprintf("%s/21_dge/%s/output.dge.reads.txt", diri, sid)
    fi2 = sprintf("%s/21_dge/%s/output.dge.umis.txt", diri, sid)
    ti1 = read_tsv(fi1)
    ti2 = read_tsv(fi2)
    #
    tu1 = ti1 %>% rename(gid=1) %>% gather(sid, n_read, -gid)
    tu2 = ti2 %>% rename(gid=1) %>% gather(sid, n_umi, -gid)
    tu = tu1 %>% inner_join(tu2, by=c('gid','sid'))
    tu
    #}}}
}

tu1 = read_dge('b01') %>% mutate(batch='b01')
tu2 = read_dge('b02') %>% mutate(batch='b02')
tu3 = read_dge('b03') %>% mutate(batch='b03')
tu = tu1 %>% bind_rows(tu2) %>% bind_rows(tu3)

t_rc = tu %>% dplyr::rename(coord=sid) %>%
    inner_join(th, by=c('batch','coord')) %>%
    select(SampleID, gid, ReadCount = n_umi)
res = readcount_norm(t_rc)
tm_u = res$tm
t_rc = tu %>% dplyr::rename(coord=sid) %>%
    inner_join(th, by=c('batch','coord')) %>%
    select(SampleID, gid, ReadCount = n_read)
res = readcount_norm(t_rc)
tm_r = res$tm

res = list(tu=tu, th=th, tm_u=tm_u, tm_r = tm_r)
fo = file.path(dirw, '00.dge.rds')
saveRDS(res, fo)

#{{{ fcnt
fi = file.path(diri, '22.fcnt.rds')
tc = readRDS(fi)

tc %>% group_by(sid) %>%
    summarise(total=sum(cnt), ng=sum(cnt>=10)) %>% ungroup()
#}}}
#}}}

fi = file.path(dirw, '00.dge.rds')
res = readRDS(fi)

#{{{ ngene umi
l2n = 1:8; names(l2n) = LETTERS[l2n]
min_nr = 2; min_nu = 2 #min_nr = 10; min_nu = 10
tp = res$tu %>% dplyr::rename(SampleID = sid) %>% filter(batch == sid) %>%
    group_by(SampleID) %>%
    summarise(nr = sum(n_read >= min_nr), nu = sum(n_umi >= min_nu)) %>%
    ungroup() %>%
    mutate(yl = str_sub(SampleID, 1, 1), x = str_sub(SampleID, 2)) %>%
    mutate(y = as.integer(l2n[yl]), x = as.integer(x))

mid = (min(tp$nu) + max(tp$nu)) / 2
tit=sprintf('# genes with >= %d UMIs', min_nu)
p = ggplot(tp) +
    geom_tile(aes(x,y,fill=nu), color='black') +
    geom_text(aes(x,y, label=number(nu), color=nu<mid), hjust=1, size=2.5, nudge_x=.3) +
    scale_x_continuous(breaks=1:12, position='top', expand=expand_scale(mult=c(.01,.01))) +
    scale_y_reverse(breaks=1:8, labels=LETTERS[1:8], expand=expand_scale(mult=c(.01,.01))) +
    scale_fill_viridis(option='viridis', direction=-1) +
    scale_color_manual(values=c('white','black')) +
    otheme(xtext=T, ytext=T, legend.pos='none') +
    theme(panel.border = element_blank()) +
    ggtitle(tit) +
    theme(plot.title=element_text(hjust=.5))
fo = sprintf("%s/05.ngene.umi.%s.pdf", dirw, sid)
ggsave(fo, p, width=6, height=4)
#}}}

tid = 'set1'; sids = 'b01'
th = res$th %>% filter(batch %in% sids) %>%
    mutate(lab = str_c(Genotype, Treatment, Timepoint, sep='_')) %>%
    mutate(grp = str_c(Genotype, Treatment, sep='_'))
tm = res$tm_r %>% filter(SampleID %in% th$SampleID)

# set1
#{{{ hclust
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .5) %>% pull(gid)
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
    scale_x_continuous(expand = expand_scale(0,.2)) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+%
    tp + geom_tiplab(aes(label=lab, color=Genotype), size=2.5) +
    scale_color_aaas()
fo = sprintf("%s/21.hclust.set1.pdf", dirw)
ggsave(p1, filename = fo, width=6, height=8)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .6) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=6,
              pca = T, max_iter = 1200)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp, aes(x=V1,y=V2)) +
    #geom_mark_ellipse(aes(fill=grp,label=grp),
        #expand=unit(3,'mm'), alpha=0, size = .2,
        #con.type='none',label.fontsize=8,label.minwidth=unit(0,'mm'),
        #label.buffer=unit(0,'mm'),label.margin = margin(0,0,0,0,"mm")) +
    geom_point(aes(color=Genotype,shape=Treatment), size=2) +
    geom_text_repel(aes(label=sprintf("%g", Timepoint)), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(0:5)) +
    scale_color_aaas() +
    scale_fill_viridis_d() +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(fill=F)
fp = sprintf("%s/25.tsne.set1.pdf", dirw)
ggsave(p_tsne, filename = fp, width=8, height=8)
#}}}

# set2
tid = 'set2'; sids = c('b02','b03')
th = res$th %>% filter(batch %in% sids) %>%
    mutate(lab = str_c(Genotype, Treatment, Timepoint, Rep, sep='_')) %>%
    mutate(grp = str_c(Genotype, Treatment, Timepoint, sep='_'))
tm = res$tm_r %>% filter(SampleID %in% th$SampleID)

#{{{ hclust
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .5) %>% pull(gid)
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
    scale_x_continuous(expand = expand_scale(0,.2)) +
    scale_y_discrete(expand = c(.01,0))
p1 = p1 %<+%
    tp + geom_tiplab(aes(label=lab, color=Genotype), size=2.5) +
    scale_color_aaas()
fo = sprintf("%s/21.hclust.set2.pdf", dirw)
ggsave(p1, filename = fo, width=9, height=15)
#}}}

#{{{ tSNE
require(Rtsne)
tw = tm %>% select(SampleID, gid, CPM) %>% mutate(CPM=asinh(CPM)) %>% spread(SampleID, CPM)
t_exp = tm %>% group_by(gid) %>% summarise(n.exp = sum(CPM>=1))
gids = t_exp %>% filter(n.exp >= (ncol(tw)-1) * .6) %>% pull(gid)
tt = tw %>% filter(gid %in% gids)
dim(tt)
tsne <- Rtsne(t(as.matrix(tt[-1])), dims=2, verbose=T, perplexity=10,
              pca = T, max_iter = 1500)

tp = as_tibble(tsne$Y) %>%
    add_column(SampleID = colnames(tt)[-1]) %>%
    inner_join(th, by = 'SampleID')
x.max=max(tp$V1)
p_tsne = ggplot(tp, aes(x=V1,y=V2)) +
    #geom_mark_ellipse(aes(fill=grp,label=grp),
        #expand=unit(3,'mm'), alpha=0, size = .2,
        #con.type='none',label.fontsize=8,label.minwidth=unit(0,'mm'),
        #label.buffer=unit(0,'mm'),label.margin = margin(0,0,0,0,"mm")) +
    geom_point(aes(color=Genotype,shape=Treatment), size=2) +
    geom_text_repel(aes(label=sprintf("%dh_%d", Timepoint, Rep)), size=2) +
    scale_x_continuous(name = 'tSNE-1') +
    scale_y_continuous(name = 'tSNE-2') +
    scale_shape_manual(values = c(0:5)) +
    scale_color_aaas() +
    scale_fill_viridis_d() +
    otheme(legend.pos='top.right', legend.dir='v', legend.title=T,
           xtitle=T, ytitle=T,
           margin = c(.2,.2,.2,.2)) +
    theme(axis.ticks.length = unit(0, 'lines')) +
    guides(fill=F)
fp = sprintf("%s/25.tsne.set2.pdf", dirw)
ggsave(p_tsne, filename = fp, width=8, height=8)
#}}}



