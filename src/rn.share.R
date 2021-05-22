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

#{{{ jackie ASE
dirw = file.path(dird, '91_share')
yid = 'rn18g'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
gids = c("Zm00001d051166", "Zm00001d046632", "Zm00001d023985", "Zm00001d051525", "Zm00001d002079", "Zm00001d015268", "Zm00001d035283", "Zm00001d052530")

th1 = th %>% select(SampleID, Tissue, Genotype, Replicate)
to = res$ase_gene %>% filter(gid %in% gids) %>%
    inner_join(th1, by='SampleID') %>%
    filter(Genotype == 'B73xMo17') %>%
    select(Tissue, Genotype, Replicate, gid, allele1, allele2)

fo = file.path(dirw, 'ase_8genes.tsv')
write_tsv(to, fo)

th1 = th %>% select(SampleID, Tissue, Genotype, Replicate)
to = res$ase_gene %>%
    inner_join(th1, by='SampleID') %>%
    filter(Genotype == 'B73xMo17') %>%
    select(Tissue, Replicate, gid, allele1, allele2)
fo = file.path(dirw, 'ase_all.tsv')
write_tsv(to, fo)
#}}}

#{{{ rn18g,rnc01 - maize1-maize2 correlation w. jackie & Erika
ts = read_syn(opt=3) %>% filter(!is.na(maize1), !is.na(maize2)) %>%
    select(maize1,maize2)
yid = 'rn18g'
#yid = 'rnc01'
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

ths = th_m %>% filter(Genotype=='B73') %>% select(SampleID, Tissue)
to = tm_m %>%
    inner_join(ths, by='SampleID') %>%
    select(gid, Tissue, CPM) %>%
    arrange(gid, Tissue) %>%
    group_by(gid) %>% summarise(cpm = list(CPM)) %>% ungroup()

calc_cor <- function(v1, v2) {
    #{{{
    r = cor.test(v1, v2)
    list(pcc=as.double(r$estimate), pcc.pval=r$p.value)
    #}}}
}
to2 = ts %>% inner_join(to, by=c('maize1'='gid')) %>% rename(cpm1=cpm) %>%
    inner_join(to, by=c('maize2'='gid')) %>% rename(cpm2=cpm) %>%
    mutate(r = map2(cpm1, cpm2, calc_cor)) %>%
    mutate(pcc = map_dbl(r, 'pcc')) %>%
    mutate(pcc.pval = map_dbl(r, 'pcc.pval')) %>%
    select(-cpm1, -cpm2, -r)
skim(to2$pcc)

to3 = ts %>% inner_join(to, by=c('maize1'='gid')) %>% rename(cpm1=cpm) %>%
    inner_join(to, by=c('maize2'='gid')) %>% rename(cpm2=cpm) %>%
    mutate(cpm1 = sample(cpm1)) %>%
    mutate(r = map2(cpm1, cpm2, calc_cor)) %>%
    mutate(pcc = map_dbl(r, 'pcc')) %>%
    mutate(pcc.pval = map_dbl(r, 'pcc.pval')) %>%
    select(-cpm1, -cpm2, -r)
skim(to3$pcc)

dirw = file.path(dird, "91_share")
fo = file.path(dirw, "cor_maize1_maize2.tsv")
write_tsv(to2, fo)
#}}}

#{{{ zm.rn19i - widiv
yid = 'zm.rn19i'
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm

to1 = tm %>% select(gid, SampleID, rTPM) %>% spread(SampleID, rTPM)
to2 = tm %>% select(gid, SampleID, TPM) %>% spread(SampleID, TPM)
to3 = tm %>% select(gid, SampleID, ReadCount) %>% spread(SampleID, ReadCount)

diro = "~/projects/s3/zhoup-share/rnaseq/rn19i"
fo = file.path(diro, "widiv942_rTPM.tsv.gz")
write_tsv(to1, fo)
fo = file.path(diro, "widiv942_TPM.tsv.gz")
write_tsv(to2, fo)
fo = glue("{diro}/widiv942_ReadCount.tsv.gz")
write_tsv(to3, fo)
#}}}

#{{{ zm.rn19i2 - widiv
yid = 'zm.rn19i2'
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm

to1 = tm %>% select(gid, SampleID, rTPM) %>% spread(SampleID, rTPM)
to2 = tm %>% select(gid, SampleID, TPM) %>% spread(SampleID, TPM)

diro = "/home/springer/zhoux379/projects/s3/zhoup-share/rnaseq/rn19i2"
fo = file.path(diro, "widiv304_rTPM.tsv.gz")
write_tsv(to1, fo)
fo = file.path(diro, "widiv304_TPM.tsv.gz")
write_tsv(to2, fo)
#}}}



