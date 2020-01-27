source("functions.R")
t_cfg

#{{{ one gene
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
#}}}

#{{{ rn18g - TF expression w. Erika
yid = 'rn18g'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

ft = '~/projects/grn/data/09.tf.txt'
tids = read_tsv(ft, col_names = F) %>% pull(X1)

tiss23 = c("radicle_root", "seedling_root", "seedling_leaf", "seedling_meristem", "coleoptile", "auricle", "blade_leaf", "internode", "tassel", "sheath", "ear", "flag_leaf", "floret", "husk", "root", "silk", "spikelet", "tassel_stem", "endosperm14D", "kernel", "embryo", "seed_imbibed", "endosperm27D")
ths = th_m %>% filter(Genotype=='B73') %>% select(SampleID, Tissue) %>%
    mutate(Tissue=factor(Tissue, levels=tiss23))
tot = tm_m %>% filter(gid %in% tids) %>%
    inner_join(ths, by='SampleID') %>%
    select(gid, Tissue, CPM) %>%
    spread(Tissue, CPM)
toa = tm_m %>%
    inner_join(ths, by='SampleID') %>%
    select(gid, Tissue, CPM) %>%
    spread(Tissue, CPM)

res = list(tf = tot, gene = toa)
dirw = file.path(dird, "91_share")
fo = file.path(dirw, "rn18g_cpm.rds")
saveRDS(res, fo)

fs = '~/projects/genome/data/Zmays_B73/gene_mapping/synteny_maize.tsv'
ts = read_tsv(fs) %>% filter(maize1 %in% tids | maize2 %in% tids)

tf = data.frame(tot[,-1])
rownames(tf) = tot$gid
#}}}

#{{{ check Jackie's genotypes
dirw = file.path(dird, '91_share')
fi = file.path(dirw, 'genotypes.xlsx')
ti = read_xlsx(fi)

tc = t_cfg %>% filter(yid %in% c('rn13c','rn13e','rn14a','rn15a','rn17a','rn18a','rn18d','rn19a','rn99a')) %>%
    select(-ref, -ase, -stress, -lgd)

tc2 = tc %>% select(yid) %>%
    mutate(fm = sprintf("%s/11_qc/%s/meta.tsv", dird, yid)) %>%
    mutate(tm = map(fm, read_tsv)) %>%
    select(-fm) %>% unnest() %>%
    distinct(yid, Genotype) %>%
    mutate(gt = str_to_upper(Genotype))

tc3 = tc2 %>% filter(gt %in% str_to_upper(ti$Genotype)) %>%
    count(yid)

to = tc2 %>% count(yid) %>% rename(n_total=n) %>%
    left_join(tc3, by='yid') %>% rename(n_overlap=n) %>%
    inner_join(tc, by='yid')
#}}}


