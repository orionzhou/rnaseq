source("functions.R")
require(lubridate)

#{{{ me18a
study = 'me18a'
dirw = file.path(dird, '08_raw_output', study)
diri = file.path(dirw, 'multiqc_data')
#
fi = file.path(dirw, 'featurecounts.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
t_rc = tcw
#
fo = file.path(dirw, 'featurecounts.tsv')
write_tsv(t_rc, fo)
#}}}

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


