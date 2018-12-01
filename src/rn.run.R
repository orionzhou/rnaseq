source("functions.R")
source(file.path(dirr, "sra.R"))
t_cfg

sid = 'mec03'
#{{{ config
Sys.setenv(R_CONFIG_ACTIVE = sid)
dirw = file.path(dird, '11_qc', sid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))
#
cfg = t_cfg %>% filter(sid == !!sid)
stopifnot(nrow(cfg) == 1)
study = cfg %>% pull(study)
readtype = cfg %>% pull(readtype)
mapper = cfg %>% pull(mapper)
genome = cfg %>% pull(reference)
meta = cfg %>% pull(meta)
x = load(file.path(dirg, genome, '55.rda'))
#}}}

#{{{ [meta=T] collect featurecounts data & normalize
sids = str_split(config::get("sids"), "[\\+]")[[1]] 
sids
#
th = tibble(); t_rc = tibble()
for (sid1 in sids) {
    diri = file.path(dird, '08_raw_output', sid1, 'multiqc_data')
    th1 = get_read_list(dird, sid1)
    if(sid1 == 'me12a') {
        th1 = th1 %>% filter(Treatment == 'WT')
    } else if(sid1 == 'me13b') {
        th1 = th1 %>% filter(!str_detect(Treatment, "ET"))
    } else if(sid1 == 'me18b') {
        th1 = th1 %>% filter(Genotype == 'B73')
    } else if(sid1 == 'me99b') {
        th1 = th1 %>% filter(Genotype == 'B73')
    }
    th1 = th1 %>% mutate(sid = sid1) %>% 
        select(sid, SampleID, Tissue, Genotype, Treatment, everything())
    th = rbind(th, th1)
    fi = file.path(diri, '../featurecounts.tsv')
    t_rc1 = read_tsv(fi) %>% select(one_of(c('gid', th1$SampleID)))
    stopifnot(ncol(t_rc1) == nrow(th1) + 1)
    if(nrow(t_rc) == 0) {
        t_rc = t_rc1
    } else {
        t_rc = t_rc %>% inner_join(t_rc1, by = 'gid')
    }
}
dim(th); dim(t_rc)
th %>% dplyr::count(sid)
th %>% dplyr::count(Tissue) %>% print(n=40)

tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm

res = merge_reps(th, tm, sid)
th = res$th; tm = res$tm
th = th %>% mutate(Replicate = NA)
tiss = unique(th$Tissue); genos = unique(th$Genotype); treas = unique(th$Treatment)
reps = unique(th$Replicate)
#}}}

#{{{ save to 20.rc.norm.rda
fo = file.path(dirw, '20.rc.norm.rda')
save(th, tl, tm, file = fo)
#}}}

