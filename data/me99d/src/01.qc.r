require(tidyverse)
dirp = '~/projects/3rnaseq'
dirw = file.path(dirp, 'data')

#{{{ read sample meta and create 01.reads.tsv
fi = file.path(dirw, "00.samplelist.tsv")
ti = read_tsv(fi)

th = ti %>% 
    separate(treatment, c('rep','day','time','trt'), sep = "_") %>%
    transmute(SampleID = sid,
              Tissue = tissue,
              Genotype = genotype,
              Treatment = sprintf("%s_%s_%s", day, time, trt),
              Replicate = rep,
              paired = F) %>%
    arrange(SampleID)
th %>% count(Replicate)
th %>% count(Genotype)
th %>% count(Tissue)

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
#}}}

# read feature-counts table
fi = file.path(dirp, 'cache/24.featurecounts/01.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
t_3rs = ti[,-c(2:6)]
colnames(t_3rs) = c('gid', sids)

fo = file.path(dirw, "10.RawReadCount.tsv")
write_tsv(t_3rs, fo)

# read briggs RNA-Seq data
fi = '~/projects/briggs/data/41.qc/10.byrep.RData'
x = load(fi)

# look at B73 leaf sample first
ttw = t_3rs
ttl = ttw %>% gather(sid, RawReadCount, -gid)
tr1 = ttl %>% filter(sid == 'A01') %>%
    mutate(CPM = RawReadCount/sum(RawReadCount) * 1000000)
tr2 = t_byrep %>% 
    filter(Genotype == "B73", Tissue == 'seedlingleaf_11DAS', Treatment == 1) %>%
    transmute(gid, RawReadCount2 = ReadCount, FPM = FPM, FPKM = FPKM)

tp = tr1 %>% inner_join(tr2, by = 'gid') %>%
    filter(RawReadCount>=5, RawReadCount2>=5)

cor(tp$RawReadCount, tp$RawReadCount2)
cor(tp$RawReadCount, tp$FPM)
cor(tp$RawReadCount, tp$FPKM)
cor(tp$CPM, tp$FPM)
cor(tp$CPM, tp$FPKM)
