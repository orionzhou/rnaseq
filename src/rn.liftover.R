source('functions.R')
dirw = '~/projects/rnaseq/data/raw'

yid = 'rn13c'
fi = sprintf("%s/%s/featurecounts.tsv", dirw, yid)
ti = read_tsv(fi)
to = ti %>% gather(sid, cnt, -gid) %>% select(sid, gid, cnt)
fo = file.path(dirw, yid, 'featurecounts.rds')
saveRDS(to, file=fo)
