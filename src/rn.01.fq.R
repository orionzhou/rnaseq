source("functions.R")
t_cfg %>% select(yid, study, author) %>% print(n=80)

yid = 'mem13'
fi = sprintf("%s/03_sra_list/%s.csv", dird, yid)
fi2 = sprintf("%s/03_sra_list/%s_exp.csv", dird, yid)
ti = read_sra_run(fi, fi2)
ti %>% print(width=Inf)
t_cfg %>% filter(yid == !!yid) %>% print(width=Inf)

th = fix_read_list(ti, yid)
th %>% count(Tissue, Genotype, Treatment) %>% print(n=100)
th %>% count(Replicate)
fo = sprintf("%s/05_read_list/%s.tsv", dird, yid)
write_tsv(th, fo)

#{{{ briggs
sid = 'me99b'
fi = '~/projects/briggs/data/03_collect/01.reads.tsv'
ti = read_tsv(fi)

tismap = c(
'spikelets_0DAP' = 'spikelet_0DAP',
'blade_v12' = 'leaf_blade_v12',
'auricle_v12' = 'stem_auricle_v12',
'sheath_v12' = 'stem_sheath_v12',
'husk_0DAP' = 'stem_husk_0DAP',
'tasselstem_0DAP' = 'stem_tasselstem_0DAP',
'flagleaf_0DAP' = 'leaf_flag_0DAP',
'coleoptile_tip' = 'coleoptile_tip',
'radicle_root' = 'root_radicle',
'seedlingleaf_11DAS' = 'leaf_seedling_11DAS',
'seedlingmeristem_11DAS' = 'meristem_seedling_11DAS',
'seedlingroot_11DAS' = 'root_seedling_11DAS'
)
th = ti %>% select(SampleID, Tissue, Genotype, r1=R1,r2=R2, spots, avgLength) %>%
    mutate(Tissue = ifelse(Tissue %in% names(tismap), tismap[Tissue], Tissue)) %>%
    separate(Tissue, c("Tissue", "Treatment"), fill='right', extra='merge') %>%
    mutate(Replicate = '', paired = T) %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired, everything())
th %>% count(Tissue, Treatment) %>% print(n=23)
th = sra_fill_replicate(th)

fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)
#}}}

#{{{ biomap
sid = 'me99c'
diri = '~/projects/biomap/data/01_exp_design'

seqdirs = sprintf('/home/hirschc1/shared/Archive/reads/biomap/RNAseq/%s',
                  c('JGI.170615.ZeamayEProfiling_3','JGI.230715.ZeamayEProfiling_4'))
fis = sprintf("%s/jgi.%d.tsv", diri, 1:2)
tq = tibble(seqdir=seqdirs, fn=fis) %>%
    mutate(data = map(fn, read_tsv, skip=1,
               col_names=c('lib','sid2','spots','sname','cond','sequencer','runType','fname'))) %>% unnest() %>%
    separate(sname,c('gt','tis','suf'),sep='_',fill='right',extra='merge',remove=F)

tismap = c('L'='Leaf','I'='Internode','E'='Endosperm','R'='Root','S'='Shoot')
tq2 = tq %>% mutate(r0=file.path(seqdir,fname)) %>%
    mutate(r0 = str_replace(r0, '\\.anqrpt', '')) %>%
    select(sid2,sname,gt,tis,suf,cond,r0,spots) %>%
    mutate(tis = str_replace(tis, '\\?+$', '')) %>%
    mutate(tis = tismap[tis]) %>%
    mutate(gt = str_replace(gt, 'xSelf$', '')) %>%
    mutate(gt = str_replace(gt, 'LH123ht', 'LH123HT')) %>%
    mutate(gt = str_replace(gt, '^OH43$', 'Oh43')) %>%
    mutate(inbred = !str_detect(gt, 'x')) %>%
    separate(gt, c('pa1','pa2'), sep="x", fill='right',extra='merge', remove=F) %>%
    mutate(pa2 = ifelse(!is.na(pa2) & pa2=='Mo17XB73','B73',pa2))
tq2 %>% count(tis)
tq2 %>% count(inbred)
tq2 %>% count(pa1) %>% print(n=40)
tq2 %>% count(pa2)

to = tq2 %>% arrange(sid2) %>%
    mutate(sid = sprintf("bm%03d", 1:length(sid2))) %>%
    mutate(Treatment=NA, Replicate=NA, paired = T, avgLength=300) %>%
    select(SampleID=sid, Tissue=tis, Genotype=gt, Treatment, Replicate, paired,
           inbred, pa1, pa2, r0, spots, avgLength) %>%
    mutate(fv = file.exists(r0))
to %>% count(fv)
to = to %>% select(-fv)

gts = to %>% filter(inbred) %>% distinct(Genotype) %>% pull(Genotype)
pas = unique(to$pa1, to$pa2)
pas = pas[!is.na(pas)]
pas[! pas %in% gts]

to = sra_fill_replicate(to)

fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(to, fo, na = '')
#}}}

#{{{ me17b, me99f
sid2pid = c('me17b'='sp033', 'me99f'='sp064')
sid = 'me99f'
sid = 'me17b' # lack fastqc report
fi = file.path(dird, '03.raw.xlsx')
ti = read_xlsx(fi, sheet=sid2pid[sid])

ndig = floor(log10(nrow(ti))) + 1
if(is.na(ti$SampleID[1])) {
    fmt_sid = sprintf("s%%0%dd", ndig)
    nsids = sprintf(fmt_sid, 1:nrow(ti))
    ti$SampleID = nsids
}
dir0 = '/home/springer/data_release/umgc/hiseq'
diro1 = sprintf("%s/%s/09_fastq_raw", dirc, sid)
diro2 = sprintf("%s/%s/10_fastq", dirc, sid)
ti = ti %>% fill(Tissue, Genotype, directory) %>%
    mutate(f1 = sprintf("%s/%s/%s_R1_001.fastq", dir0, directory, file),
           f2 = sprintf("%s/%s/%s_R2_001.fastq", dir0, directory, file)) %>%
    mutate(gz = ifelse(file.exists(f1), F, T)) %>%
    mutate(f1 = ifelse(gz, sprintf("%s.gz", f1), f1),
           f2 = ifelse(gz, sprintf("%s.gz", f2), f2)) %>%
    mutate(nf1 = ifelse(gz,
                        sprintf("%s/%s_1.fq.gz", diro2, SampleID),
                        sprintf("%s/%s_1.fq", diro1, SampleID)),
           nf2 = ifelse(gz,
                        sprintf("%s/%s_2.fq.gz", diro2, SampleID),
                        sprintf("%s/%s_2.fq", diro1, SampleID))) %>%
    mutate(cmd1 = sprintf("ln -sf %s %s", f1, nf1),
           cmd2 = sprintf("ln -sf %s %s", f2, nf2)) %>%
    mutate(tag = file.exists(f1) & file.exists(f2))
sum(!ti$tag)

if(sid == 'me17b') {
    to = ti %>% mutate(spots = 5e7, avgLength=51)
} else {
tt = ti %>% distinct(directory) %>%
    mutate(rdir = file.path(dir0, directory)) %>%
    mutate(data = map(rdir, read_msi_fastqc)) %>%
    unnest() %>% select(-rdir)
ti2 = ti %>% mutate(sampleName = str_replace_all(file,'_','.')) %>%
    mutate(sampleName = str_replace(sampleName, '\\.S[0-9]+$',''))
to = ti2 %>% inner_join(tt, by=c('directory','sampleName'))
}

th = to %>% select(SampleID, Tissue, Genotype, r1=f1, r2=f2, spots,avgLength) %>%
    mutate(Treatment='', Replicate = '', paired = T, avgLength=avgLength*2) %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired, everything())
th = sra_fill_replicate(th)
th %>% count(Genotype, Tissue, Treatment) %>% print(n=50)

fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)
#}}}

#{{{ tenders_3rnaseq: me99d
sid = 'me99d'
fi = file.path(dird, '03.raw.xlsx')
ti = read_xlsx(fi, sheet=sid)
#
diro2 = sprintf("%s/%s/10_fastq", dirc, sid)
ti = ti %>% fill(Tissue, Genotype, directory) %>%
    mutate(f1 = sprintf("%s/%s.fq.gz", directory, file)) %>%
    mutate(nf1 = sprintf("%s/%s.fq.gz", diro2, SampleID)) %>%
    mutate(cmd1 = sprintf("ln -sf %s %s", f1, nf1)) %>%
    mutate(tag = file.exists(f1))
sum(!ti$tag)

#{{{ [obsolte] create sym-links, write read list
if(!dir.exists(diro2)) system(sprintf("mkdir -p %s", diro2))
map_int(ti$cmd1, system)
#}}}

th = ti %>% select(SampleID,Tissue,Genotype,Treatment,r0=f1) %>%
    mutate(Replicate = '', paired = F) %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired, everything())
th = sra_fill_replicate(th)
th %>% count(Genotype, Tissue, Treatment) %>% print(n=50)

fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)
#}}}


