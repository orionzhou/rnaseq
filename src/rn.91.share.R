source("functions.R")
t_cfg

gid = 'Zm00001d042315'

yid = 'me13a'
yid = 'me99b'
fc = sprintf("%s/08_raw_output/%s/cpm.rds", dird, yid)
x = readRDS(fc)

to = x$tm_m %>% filter(gid == !!gid) %>%
    inner_join(x$th_m, by='SampleID') %>%
    mutate(Tissue = sprintf("%s_%s", Tissue, Treatment)) %>%
    select(gid, Tissue, Genotype, ReadCount, CPM, FPKM) %>%
    arrange(gid, Tissue, Genotype)
fo = sprintf("%s/91_share/%s.tsv", dird, yid)
write_tsv(to, fo)
