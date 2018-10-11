#{{{ header
source("me.fun.r")
t_cfg
fi = file.path('~/data/genome/B73', "v37/t2.tsv")
t_gs = read_tsv(fi, col_types = 'ccccciic') %>% 
    filter(etype == 'exon') %>% 
    group_by(gid, tid) %>% 
    summarise(size = sum(end - beg + 1)) %>%
    group_by(gid) %>%
    summarise(size = max(size))
#}}}

sid = 'me99c'

#{{{ mapping stats
Sys.setenv(R_CONFIG_ACTIVE = sid)
diri = file.path(dird, '08_raw_output', sid, 'multiqc_data')
dirw = file.path(dird, '11_qc', sid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))
th = get_read_list(dird, sid)
tiss = unique(th$Tissue); genos = unique(th$Genotype); treas = unique(th$Treatment)
reps = unique(th$Replicate)
#
tt = read_multiqc(diri, th)
fo = file.path(dirw, '10.mapping.stat.tsv')
write_tsv(tt, fo)
#}}}

#{{{ obtain raw read counts, normalize and save
fi = file.path(diri, '../featurecounts.tsv')
t_rc = read_tsv(fi)
#
tm = t_rc %>% gather(SampleID, ReadCount, -gid)
res = readcount_norm(tm, t_gs)
tl = res$tl; tm = res$tm
#
fo = file.path(dirw, '20.rc.norm.rda')
save(tl, tm, file = fo)
#}}}


