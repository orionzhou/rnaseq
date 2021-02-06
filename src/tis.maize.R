source("functions.R")
gcfg = read_genome_conf()
tsyn = read_syn(gcfg)
dirw = glue("{dird}/15_tissue_spec")

fi = file.path(dirw, 'cpm_combined.rds')
r = readRDS(fi)
th = r$th; tm = r$tm

#{{{ expression broadth
calc_tau <- function(x) {
    #{{{
    x[x <= 1] = 1
    x1 = log(x)
    xmax = max(x1)
    if (xmax == 0) {
        NA
    } else {
        sum(1 - x1/xmax) / (length(x) - 1)
    }
    #}}}
}

th0 = th %>% select(sid=SampleID, tis=Tissue2)
to0 = tm %>% select(gid,sid=SampleID,cpm=CPM) %>%
    mutate(silent = cpm < 1) %>%
    inner_join(th0, by='sid') %>%
    group_by(gid) %>%
    summarise(n.tis.tot = n(), n.tis = sum(!silent),
              tiss = str_c(unique(tis), collapse=','), cpms = list(cpm)) %>%
    ungroup() %>%
    mutate(tau = map_dbl(cpms, calc_tau)) %>%
    mutate(prop.tis = n.tis / n.tis.tot)

#{{{ prop.tis
tags = c('Silent', 'Tissue specific', 'Intermediate frequency', 'Constitutive')
to1 = to0 %>% select(gid,prop.tis) %>%
    mutate(tag = ifelse(prop.tis == 0, tags[1],
           ifelse(prop.tis <= 0.1, tags[2],
           ifelse(prop.tis < 0.9, tags[3], tags[4])))) %>%
    mutate(tag = factor(tag, levels = tags)) %>%
    select(gid, etag=tag)
to1 %>% count(etag)
#}}}
#{{{ tau
itv = .02
itvs = seq(0,1, by=itv)
tags = c('constitutive', 'intermediate', 'specific', 'silent')
to2 = to0 %>% select(gid,tau) %>%
    mutate(bin.tau = cut(tau,breaks=itvs, right=T, include.lowest=T, ordered_result=T)) %>%
    mutate(tau2 = as.integer(bin.tau) * itv) %>%
    mutate(tag = ifelse(tau>.9, tags[3], ifelse(tau<=.5, tags[1], tags[2]))) %>%
    replace_na(list(tag=tags[4])) %>%
    mutate(tag = factor(tag, levels=tags)) %>%
    select(gid, bin.tau, tau2, tag)
to2 %>% count(tag)
#}}}

to = to0 %>% inner_join(to1, by='gid') %>% inner_join(to2,by='gid') %>%
    select(gid,tau,bin.tau,tau2,tag,n.tis.tot,n.tis,prop.tis,etag,tiss)
to %>% count(tag,etag) %>% spread(etag,n)

fo = glue("{dirw}/10.tau.rds")
saveRDS(to, fo)
fo = glue("{dirw}/10.tau.tsv.gz")
write_tsv(to, fo, na='')

#{{{ plot
tp = to3 %>% count(bin.tau,tau2,tag)
tpx = tp %>% filter(tau2 %in% c(.2,.4,.6,.8)) %>% mutate(lab=tau2)
    #mutate(lab = c("<- narrow", "broad ->"))
p = ggplot(tp) +
    geom_col(aes(x=bin.tau, y=n, fill=tag), width=.8) +
    scale_x_discrete(name='Tau', breaks=tpx$bin.tau, labels=tpx$lab, expand=expansion(mult=c(.08,.02))) +
    scale_y_continuous(name="Number Genes", expand=expansion(mult=c(0,.03))) +
    scale_fill_npg() +
    otheme(legend.pos='top.left', 
           xtitle=T, ytitle=T, xtext=T, ytext=T, xtick=T, ytick=T,
           margin = c(.3,.3,.3,.3), "lines", legend.title=F)
#
fo = file.path(dirw, '11.tau.pdf')
ggsave(fo, p, width=4, height=4)
#}}}
#}}}

ft = glue("{dirw}/10.tau.rds")
tt = readRDS(ft)
tf = read_tf()

#{{{ synteny
tp = tt %>% inner_join(tsyn, by='gid') %>%
    select(gid, tag1=ftype, tag2=tag) %>% count(tag1, tag2)
p = cmp_proportion1(tp, xangle=10, ytext=T, legend.title='',
                    lab.size=2, alph=.6, fills = pal_lancet()(5)) +
    o_margin(.1,.3,.1,.3) +
    theme(legend.position='none')
fo = file.path(dirw, '13.syn.pdf')
ggsave(p, file=fo, width=4, height=4)

tp = tsh_e %>% filter(gid %in% tf$gid) %>%
    inner_join(tsyn, by='gid') %>%
    select(gid, tag1=ftype, tag2=etag) %>% count(tag1, tag2)
p1 = cmp_proportion(tp, xangle=10, legend.pos='left', legend.dir='v')
fo = file.path(dirw, 'syn.tis.tf.pdf')
ggsave(p1, file=fo, width=6, height=6)


tp = tsh_e %>% inner_join(tsyn, by='gid')
p1 = ggplot(tp) +
    geom_density(aes(tau, color=ftype), alpha=1) +
    scale_fill_npg() +
    scale_color_npg() +
    otheme(legend.pos='top.left')
fo = file.path(dirw, 'syn.density.pdf')
ggsave(p1, file=fo, width=6, height=6)
#}}}


