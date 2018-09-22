require(tidyverse)
source("me.fun.r")

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
