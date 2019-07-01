source("functions.R")

yid = 'rn18e'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ read in, filter/fix samples
ref = t_cfg %>% filter(yid == !!yid) %>% pull(ref)
th = rnaseq_sample_meta(yid)
tt = rnaseq_mapping_stat(yid)
#res = rnaseq_cpm_raw(yid)
#th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m
sum_stat_tibble(tt)

sids_keep = tt %>% filter(mapped>0) %>% pull(SampleID)
sum_stat_tibble(tt %>% filter(SampleID %in% sids_keep))

# fix th
th2 = th %>% filter(SampleID %in% sids_keep)


th = th2
tt = tt %>% filter(SampleID %in% th$SampleID)

fh = file.path(dirw, 'meta.tsv')
write_tsv(th, fh, na='')
# run snakemake again
#}}}





