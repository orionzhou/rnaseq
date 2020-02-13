source("functions.R")
t_cfg %>% print(n=50)

yid = 'rn17c'
yid = 'rn18g'
diri = file.path(dird, 'raw', yid)
fi = file.path(diri, 'cpm.rds')
x = readRDS(fi)
fi = file.path(diri, 'ase.rds')
ase = readRDS(fi)
th = x$th %>% select(SampleID,Tissue,Genotype,Treatment)
fi = file.path(diri, 'ase0.tsv')
ase0 = read_tsv(fi)
fi = file.path(diri, 'bamstats.tsv')
bs = read_tsv(fi)

ase1 = ase %>%
    separate(sid, c('sid','pa'), sep='[\\.]') %>%
    spread(pa, cnt) %>%
    rename(pa1 = `1`, pa2 = `2`) %>%
    replace_na(list(pa1=0, pa2=0))

ase1 %>% group_by(sid) %>%
    summarise(prop.B = sum(pa1) / sum(pa1+pa2)) %>%
    ungroup() %>%
    inner_join(th, by = c('sid'='SampleID')) %>%
    print(n=120)
