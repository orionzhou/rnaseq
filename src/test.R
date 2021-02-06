source("functions.R")
require(lubridate)
tgt = 'B73'
qry = 'Mo17'
gconf1 = read_genome_conf(qry)
gconf2 = read_genome_conf(tgt)
#
fv = sprintf('~/projects/wgc/data/05_stats/10.%s_%s.tsv', tgt,qry)
tv = read_tsv(fv)
gids1 = gconf1$t_gs$gid
gids2 = gconf2$loc.gene %>% filter(ttype=='mRNA') %>% distinct(gid) %>% pull(gid)
gids1[gids1 %in% gids2]
t_pa = tibble(pa = rep(c('m','b'), c(length(gids1),length(gids2))),
    gid = c(gids1, gids2))

#{{{ july
ta = tibble(date1 = ymd('2018-11-3') + 0:6, Freq = c(5,5,0,1,0,1,1)) %>%
    mutate(Date = sprintf("%s %s", month(date1, label=T), mday(date1)))
p = ggplot(ta) +
    geom_bar(aes(x=Date, y=Freq, fill=as.character(Date)), stat='identity', width=.8) +
    scale_x_discrete(breaks = ta$Date, labels = ta$Date) +
    scale_y_continuous(name = 'Frequency of Angry') +
    scale_fill_npg() +
    otheme(legend.pos = 'none',
    xtext=T, xtitle=T, ytext=T, ytitle=T, xgrid=F, ygrid=T,
    xtick = T, ytick = T)
fo = file.path(dird, 'july.pdf')
ggsave(p, file = fo, width = 5, height = 4)
#}}}

concat_gid <- function(gid1, gid2)
    ifelse(gid1=='.', gid2,
    ifelse(gid2=='.', gid1, str_c(str_sort(c(gid1, gid2)), collapse="--")))
dird = '~/projects/wgc/data'
qry = 'Mo17'; tgt = 'B73'
diri = sprintf('%s/raw_output/%s_%s/20_synteny', dird, qry, tgt)
fi = file.path(diri, '07.t.ortholog')
ti = read_tsv(fi, col_names=c('tid2','tid1')) %>%
    mutate(rbh=str_detect(tid1, "\\'$")) %>%
    mutate(tid1 = str_remove(tid1, "\\'$")) %>%
    mutate(type2 = ifelse(tid1 == '.', 'lost1', ifelse(rbh, 'rbh', 'synteny')))
tp = ti %>% mutate(gid1=str_remove(tid1, "_\\w+$"),
    gid2=str_remove(tid2, "_\\w+$")) %>%
    select(gid1,gid2,type2) %>%
    mutate(cid = map2_chr(gid1, gid2, concat_gid)) %>%
    rename(type=type2)
tp %>% count(type)

fj = file.path(diri, '07.q.ortholog')
tj = read_tsv(fj, col_names=c('tid1','tid2')) %>%
    mutate(rbh=str_detect(tid2, "\\'$")) %>%
    mutate(tid2 = str_remove(tid2, "\\'$")) %>%
    mutate(type1 = ifelse(tid2 == '.', 'lost2', ifelse(rbh, 'rbh', 'synteny')))
tp = tj %>% mutate(gid1=str_remove(tid1, "_\\w+$"),
    gid2=str_remove(tid2, "_\\w+$")) %>%
    select(gid1,gid2,type1) %>%
    mutate(cid = map2_chr(gid1, gid2, concat_gid)) %>%
    rename(type=type1)
tp %>% count(type)

tg = tv %>% filter(syn=='syntenic') %>%
    inner_join(tp, by=c('gid'='gid2'))
tg %>% count(type,impact)

tp %>% inner_join(tv, by=c('gid2'='gid')) %>% count(type,impact)

yid = 'me17c'
sid = 'SRR2043331'
fi = sprintf('~/projects/rnaseq/data/cache/%s/31_mmquant/%s.tsv', yid, sid)
ti = read_tsv(fi) %>% select(gid=Gene,rc=`sid`) %>%
    mutate(ng = str_count(gid, '--')+1)
rct = sum(ti$rc)
sum(ti$rc[ti$ng==1])/rct
sum(ti$rc[ti$ng==2])/rct
sum(ti$rc[ti$ng>2])/rct
#
ti2 = ti %>% filter(ng == 1) %>% select(-ng) %>%
    inner_join(t_pa, by = 'gid')
ti2 %>%
    group_by(pa) %>% summarise(rc = sum(rc)) %>% ungroup() %>%
    spread(pa, rc) %>% mutate(prop = b/(b+m)) %>%
    print(n=10)

gids = tg %>% filter(po >= .95, syn == 'syntenic') %>% pull(gid)
ti2 %>% filter(gid %in% gids) %>% mutate(prop.b = b/(b+m)) %>% summary


ti %>% filter(gid %in% tp$cid[tp$type=='synteny']) %>%
    pull(rc) %>% sum()/rct
ti %>% filter(gid %in% tp$gid1[tp$type=='synteny']) %>%
    pull(rc) %>% sum()/rct
ti %>% filter(gid %in% tp$gid2[tp$type=='synteny']) %>%
    pull(rc) %>% sum()/rct

impacts = c("low",'modifier','no_change','moderate')
impacts = c('modifier','moderate')
to = tg %>% filter(type=='synteny',impact %in% impacts) %>%
    rename(gid2=gid) %>% select(gid1,gid2) %>%
    inner_join(ti, by=c("gid1"="gid")) %>%
    select(-ng) %>% rename(rc1=rc) %>%
    inner_join(ti, by=c("gid2"="gid")) %>%
    select(-ng) %>% rename(rc2=rc) %>%
    filter(rc1+rc2 >= 50) %>%
    inner_join(gconf1$size.gene, by=c("gid1"='gid')) %>% rename(size1=size) %>%
    inner_join(gconf2$size.gene, by=c("gid2"='gid')) %>% rename(size2=size) %>%
    mutate(rc1=rc1/size1, rc2=rc2/size2, rc = rc1+rc2, prop.b = rc2/rc)
to %>%
    summary()


fi = '~/projects/rnaseq/data/05_read_list/me99b.c.tsv'
diri = '~/projects/rnaseq/data/cache/me99b'

ti = read_tsv(fi) %>%
    select(sid=SampleID,tis=Tissue,gt=Genotype,cond=Treatment,rep=Replicate) %>%
    mutate(fs=sprintf("%s/22_bam/%s.tsv", diri, sid)) %>%
    mutate(data=map(fs, read_tsv, col_names=c('type','rc'))) %>%
    select(-fs) %>%
    unnest()

ti %>% filter(gt %in% c("B73",'Mo17')) %>%
    mutate(gt=sprintf("%s_%d",gt,rep)) %>%
    group_by(tis, cond, gt) %>%
    summarise(prop.map = rc[type=='pair_map']/rc[type=='pair']) %>%
    ungroup() %>% spread(gt, prop.map) %>%
    filter(tis %in% c("embryo",'leaf')) %>%
    print(n=24)


ti = read_tsv(fi) %>%
    select(sid=SampleID,tis=Tissue,gt=Genotype,cond=Treatment,rep=Replicate) %>%
    mutate(fs=sprintf("%s/test/%s.tsv", diri, sid)) %>%
    mutate(data=map(fs, read_tsv, col_names=c('chrom','rc','rcu'))) %>%
    select(-fs) %>%
    unnest()

ti2 = ti %>% select(-rcu) %>%
    filter(! tis %in% c('endosperm','kernel'), !cond %in% c('blade_v12')) %>%
    filter(! gt %in% c("B73",'Mo17')) %>%
    mutate(chrom.gt=str_sub(chrom,1,1)) %>%
    mutate(chrom2=str_sub(chrom,2,3)) %>%
    select(-chrom)# %>% filter(gt=='Mo17xB73')

ti2 %>% group_by(sid,tis,gt,cond,rep,chrom.gt) %>%
    summarise(rc=sum(rc)) %>%
    spread(chrom.gt, rc) %>%
    mutate(b_m = B/M) %>% summary()

ti2 %>%
    spread(chrom.gt, rc) %>%
    mutate(b_m = B/M) %>% select(-B, -M) %>%
    spread(chrom2, b_m) %>% summary()


