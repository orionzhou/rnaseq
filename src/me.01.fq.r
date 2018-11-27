source("functions.R")
source(file.path(dirr, "sra.R"))
t_cfg %>% select(sid, study, author) %>% print(n=40)

get_read_list <- function(ti, sid) {
#{{{
if(sid == 'me10a') {
#{{{ Li2010
    th = ti %>%
        transmute(SampleID = Run,
                  Tissue = 'leaf',
                  Genotype = 'B73',
                  Treatment = SampleName,
                  Replicate = '',
                  paired = paired) %>%
        arrange(SampleID)
#}}}
} else if (sid == 'me11a') {
#{{{ Davidson2011
    tismap = c(
"Embryo 25 days after pollination" = 'embryo_25DAP',
"Endosperm 25 days after pollination" = 'endosperm_25DAP',
"Leaves 20-day old seedling - field" = 'leaf_20d_field',
"Leaves 20-day old seedling - growth chamber" = 'leaf_20d_gc',
"Mature silk" = 'silk',
"Ovule" = 'ovule',
"Pollen" = 'pollen',
"Post-emergence cob" = 'cob_post-emerg',
"Postemergence tassel" = 'tassel_post-emerg',
"Pre-emergence cob" = 'cob_pre-emerg',
"Preemergence tassel" = 'tassel_pre-emerg',
"Seed 10 days after pollination" = 'seed_10DAP',
"Seed 5 days after pollination" = 'seed_5DAP',
"Whole anthers" = 'anther')
    th = ti %>%
        separate(Title, c('pre','tis'), sep = ' B73 ') %>%
        separate(tis, c('tis', 'suf'), sep = ' RNA-Seq ') %>%
        mutate(Tissue = tismap[tis]) %>%
        separate(Tissue, c('Tissue','Treatment'), extra='merge', fill='right') %>%
        transmute(SampleID = Run,
                  Tissue = Tissue,
                  Genotype = 'B73',
                  Treatment = Treatment,
                  Replicate = '',
                  paired = paired) %>%
        arrange(SampleID)
#}}}
} else if (sid == 'me12a') {
#{{{ Bolduc2012
    th = ti %>%
        separate(LibraryName, c("gsm1", "gsm"), sep = ": ") %>%
        separate(gsm, c("gentisrep", 'zm', 'rna'), sep = "; ") %>%
        separate(gentisrep, c("gentis", 'rep'), sep = " #") %>%
        mutate(gentis = str_replace(gentis, " leaf (homo|het)", "_\\1 leaf")) %>%
        separate(gentis, c("gen", "tis"), sep = " ") %>%
        mutate(tis = str_replace(tis, 's$', '')) %>%
        transmute(SampleID = Run,
                  Tissue = tis,
                  Genotype = 'B73',
                  Treatment = gen,
                  Replicate = rep,
                  paired = paired) %>%
        arrange(SampleID)
#}}}
} else if (sid == 'me13a') {
#{{{ Li2013
    th = ti %>% separate("SampleName", c('org', 'ibm', 'parent', 'tis1', 'tis2', 'Genotype'), sep = "_", fill = 'left')
    th %>% count(parent)
    th %>% count(tis1)
    th %>% count(tis2)
    th %>% count(Genotype)
    th = th %>% transmute(SampleID = Run,
                          Tissue = 'SAM',
                          Genotype = Genotype,
                          Treatment = '',
                          Replicate = '',
                          paired = paired) %>%
        arrange(SampleID)
#}}}
} else if (sid == 'me13b') {
#{{{ Liu2013
    th = ti %>% 
        separate("SampleName", c("pre", "Treatment"), sep = "_", fill = "left") %>%
        mutate(Treatment = ifelse(is.na(pre), Treatment, sprintf("E%s", Treatment))) 
    th %>% count(paired)
    th %>% count(Treatment)
    th = th %>% transmute(SampleID = Run,
                          Tissue = 'leaf',
                          Genotype = 'B73',
                          Treatment = Treatment,
                          Replicate = '',
                          paired = paired) %>%
        arrange(SampleID)
#}}}
} else if (sid == 'me13c') {
#{{{ Eitchen2013
    th = ti %>%
        mutate(gt = str_replace(SampleName, "[ _](rep|R) ?[0-9]+", '')) %>%
        mutate(gt = ifelse(gt=='M37W','M37w',gt))
    th %>% count(paired)
    th %>% count(gt) 
    th = th %>% transmute(SampleID = Run,
                          Tissue = 'unknown',
                          Genotype = gt,
                          Treatment = '',
                          Replicate = '',
                          paired = paired) %>%
        arrange(SampleID)
#}}}
} else if (sid == 'me13d') {
#{{{ Waters2013
    th = ti %>% transmute(SampleID = Run,
                          Tissue = 'endosperm',
                          Genotype = SampleName,
                          Treatment = '',
                          Replicate = '',
                          paired = paired) %>%
        arrange(SampleID)
    th %>% count(Genotype)
    th %>% count(paired)
#}}}
} else if (sid == 'me14a') {
#{{{ Hirsch2014
th = ti %>% separate(Title, c("org", "Genotype"), sep = ", ") %>%
    separate(Genotype, c("Genotype", "suf"), sep = " RNAseq") %>%
    select(-org, -suf) %>%
    transmute(SampleID = Run,
              Tissue = "seedling",
              Genotype = Genotype,
              Treatment = '',
              Replicate = 1,
              paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me14b') {
#{{{ Li2014 endosperm
th = ti %>% 
    separate(Title, c("pre", "Treatment"), sep = " B73 ") %>%
    select(-pre) %>%
    separate(Treatment, c("Treatment", "suf"), sep = "DAP; ") %>%
    select(-suf) %>%
    mutate(Treatment = ifelse(Treatment %in% c("0a", "0b"), "0", Treatment)) %>%
    mutate(Treatment = as.integer(Treatment)) %>%
    transmute(SampleID = Run,
              Tissue = "endosperm",
              Genotype = 'B73',
              Treatment = Treatment,
              Replicate = '',
              paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me14c') {
#{{{ Chettoor gamete 7
th = ti %>%
    mutate(tisrep = str_replace(LibraryName, 'Maize', '')) %>%
    mutate(tis = str_to_lower(str_sub(tisrep, 1, -2)),
           rep = as.integer(str_sub(tisrep, -1, -1))) %>%
    mutate(tis = ifelse(tis == 'embryosac', 'embryo_sac', tis)) %>%
    separate(tis, c("tis","treat"), sep="_", fill='right',extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = tis,
              Genotype = 'B73',
              Treatment = treat,
              Replicate = rep,
              paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me14d') {
#{{{ johnston ligule
th = ti %>% separate(Title, c("gsm1", "gsm"), sep = ": ") %>%
    separate(gsm, c("tisrep", 'zm', 'rna'), sep = "; ") %>%
    mutate(tisrep = str_replace(tisrep, '[\\.]', '-')) %>%
    separate(tisrep, c('tis', 'rep'), sep = '-') %>%
    transmute(SampleID = Run,
              Tissue = 'ligule',
              Genotype = 'B73',
              Treatment = tis,
              Replicate = '',
              paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me15a') {
#{{{ Leiboff2015
th = ti %>% transmute(SampleID = Run,
                      Tissue = 'SAM',
                      Genotype = SampleName,
                      Treatment = '',
                      Replicate = '',
                      paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me15b') {
#{{{ Yu2015
    th = ti %>% separate(Title, c("pre", "Treatment"), sep = " at ") %>%
        transmute(SampleID = Run,
                  Tissue = "leaf",
                  Genotype = 'B73',
                  Treatment = Treatment,
                  Replicate = 1,
                  paired = paired) %>%
        filter(paired) %>%
        arrange(SampleID)
#}}}
} else if (sid == 'me16a') {
#{{{ Jin2016
th = ti %>% transmute(SampleID = Run,
                      Tissue = 'kernel',
                      Genotype = SampleName,
                      Treatment = '',
                      Replicate = '',
                      paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me16b') {
#{{{ Stelpflug2016
th = ti %>%
    mutate(Title = str_replace(Title, 'RNA-seq', 'RNA-Seq')) %>%
    separate(Title, c('pre', 'str1'), sep = ", B73 ", fill='left') %>%
    mutate(str1 = str_replace(str1, '^Zea mays ', '')) %>%
    mutate(str1 = str_replace(str1, ' \\(.*\\)$', '')) %>%
    separate(str1, c('tis_str', 'suf'), sep = " RNA-Seq", fill='right') %>%
    mutate(tis_str = str_replace(tis_str, '_R[1-3]$', '')) %>%
    mutate(tis_str = str_replace(tis_str, ' Rep[1-3]$', ''))
tmap = th %>% select(SampleID=Run, Tissue=tis_str) %>% 
    mutate(nTissue = '') %>% count(Tissue, nTissue)
fo = file.path(dird, '05_read_list/me16b_map_raw.tsv')
write_tsv(tmap, fo)
fo = file.path(dird, '05_read_list/me16b_map.tsv')
tmap = read_tsv(fo)
th = th %>% 
    inner_join(tmap, by = c('tis_str' = 'Tissue')) %>%
    separate(nTissue, c("Tissue",'Treatment'), sep='_', fill='right', extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = Tissue,
              Genotype = 'B73',
              Treatment = Treatment,
              Replicate = '',
              paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me16c') {
#{{{ Walley2016
tismap = c(
"2-4 mm from tip of ear primordium" = 'ear_2-4',
"6-8 mm from tip of ear primordium" = 'ear_6-8',
"Cortex" = 'root_cortex',
"EMBRYOS" = 'embryo',
"embryos_20DAP" = 'embryo_20DAP',
"endosperm_12DAP" = 'endosperm_12DAP',
"endosperm_crown" = 'endosperm_crown',
"EZ" = 'root_ez',
"Germinating Kernels" = 'kernel_germinating',
"GROWTH ZONE" = 'leaf_growth',
"Internode 6-7" = 'internode_6-7',
"Internode 7-8" = 'internode_7-8',
"mature female spikelets" = 'spikelet',
"MATURE LEAF TISSUE (leaf 8)" = 'leaf_mature_8',
"Mature pollen" = 'pollen',
"MZ" = 'root_mz',
"pericarp_aleurone" = 'seed_pericarp',
"PR" = 'root_primary',
"silks" = 'silk',
"SR" = 'root_secondary',
"STOMATAL_DIVISION_ZONE" = 'leaf_stomatal',
"SYMMETRICAL_DIVISION_ZONE" = 'leaf_symmetrical',
"Vegetative Meristem & Surrounding Tissue" = 'meristem')
th = ti %>% separate(Title, c('gsm1', 'gsm'), sep = ': ') %>%
    select(-gsm1) %>%
    separate(gsm, c('tisrep', 'suf1', 'suf2'), sep = '; ') %>%
    select(-suf1, -suf2) %>%
    mutate(tisrep = str_replace(tisrep, "_r([1-3])$", "=\\1")) %>%
    mutate(tisrep = str_replace(tisrep, "_rep([1-3])$", "=\\1")) %>%
    mutate(tisrep = str_replace(tisrep, " ([1-3])$", "=\\1")) %>%
    separate(tisrep, c('Tissue', 'Replicate'), sep = "=") %>%
    mutate(Tissue = tismap[Tissue]) %>%
    separate(Tissue, c("Tissue","Treatment"), sep='_', fill='right') %>%
    transmute(SampleID = Run,
              Tissue = Tissue,
              Genotype = 'B73',
              Treatment = Treatment,
              Replicate = Replicate,
              paired = paired) %>%
    arrange(SampleID)
th %>% count(Tissue) %>% print(n=23)
#}}}
} else if (sid == 'me17a') {
#{{{ Lin2017
th1 = ti %>% filter(paired) %>%
    mutate(LibraryName = ifelse(LibraryName=='Mo18W', 'Mo18W-root', LibraryName)) %>%
    separate(LibraryName, c('gt','tissue'), sep = "-") %>%
    mutate(tissue = ifelse(tissue == 'fieldear', 'ear', tissue)) %>%
    transmute(SampleID = Run,
              Tissue = tissue,
              Genotype = gt,
              Treatment = '',
              Replicate = '',
              paired = paired) %>%
    arrange(SampleID)
th1 %>% count(Tissue)
#
th2 = ti %>% filter(!paired) %>%
    transmute(SampleID = Run,
              Tissue = 'SAM',
              Genotype = LibraryName,
              Treatment = '',
              Replicate = '',
              paired = paired) %>% 
    arrange(SampleID)
th = rbind(th1,th2)
th %>% count(Tissue)
th %>% count(Genotype)
#}}}
} else if (sid == 'me18a') {
#{{{ Kremling2018
th1 = ti %>% 
    separate("LibraryName", c('lib1', 'lib2', 'tissue', 'genotype', 'suf'), 
             sep = "_", fill = 'left', extra = 'merge')
th1 %>% count(tissue)
tissues = "GRoot GShoot Kern L3Base L3Mid L3Tip LMAD26 LMAD8 LMAN26 LMAN8 LMid"
tissues = strsplit(tissues, " ")[[1]]
idxs = (!th1$tissue %in% tissues)
#
th2 = th1 %>% 
    mutate(
           genotype = ifelse(idxs, sprintf("%s_%s", tissue, genotype), genotype),
           tissue = ifelse(idxs, lib2, tissue),
           lib2 = ifelse(idxs, '', lib2))
th2 %>% count(tissue)
#
th3 = th2 %>% 
    mutate(tissue = ifelse(tissue %in% c("LMAD26", "LMAD8"), "LMAD", tissue)) %>%
    mutate(tissue = ifelse(tissue %in% c("LMAN26", "LMAN8"), "LMAN", tissue))
th3 %>% count(tissue)
#
tissues = "GRoot GShoot Kern L3Base L3Tip LMAD LMAN" 
tissues = strsplit(tissues, " ")[[1]]
th4 = th3 %>% filter(tissue %in% tissues) %>% 
    transmute(SampleID = Run,
              Tissue = tissue,
              Genotype = genotype,
              paired = paired) %>%
    arrange(SampleID)
#
th = th4 %>% mutate(Treatment = '', Replicate = '') %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired)
th %>% count(Tissue) %>% print(n=100)
#}}}
} else if (sid == 'me18b') {
#{{{ Baldauf2018
th = ti %>% 
    mutate(SampleName = str_replace(SampleName, '-', '_0_')) %>%
    separate(SampleName, c("gt", "stage", "rep"), sep = "_") %>%
    transmute(SampleID = Run,
              Tissue = 'root',
              Genotype = gt,
              Treatment = stage,
              Replicate = rep,
              paired = paired) %>% 
    arrange(SampleID)
#}}}
} else if (sid == 'me18d') {
#{{{ Schaefer2018
th = ti %>%
    separate(SampleName, c("gt", "suf1", "suf2"), sep="-", fill='right', extra='merge') %>%
    transmute(SampleID = Run,
              Tissue = 'root',
              Genotype = gt,
              Treatment = '',
              Replicate = '',
              paired = paired) %>%
    arrange(SampleID)
} else {
    cat("unknown study: ", sid, "\n")
}
th = sra_fill_replicate(th)
th
#}}}
#}}}
}

sid = 'me13d'
fi = sprintf("%s/03_sra_list/%s.csv", dird, sid)
fi2 = sprintf("%s/03_sra_list/%s_exp.csv", dird, sid)
ti = read_sra_run(fi, fi2)

th = get_read_list(ti, sid)
th %>% count(Tissue, Genotype, Treatment) %>% print(n=100)
th %>% count(Replicate)
fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)

#{{{ briggs
sid = 'me99b'
fi = '~/projects/briggs/data/03_collect/01.reads.tsv'
ti = read_tsv(fi)

diro1 = sprintf("%s/cache/%s/09_fastq_raw", dird, sid)
diro2 = sprintf("%s/cache/%s/10_fastq", dird, sid)
if(!dir.exists(diro1)) system(sprintf("mkdir -p %s", diro1))
if(!dir.exists(diro2)) system(sprintf("mkdir -p %s", diro2))

# create sym-links, write read list
ti = ti %>%
    mutate(nR1 = ifelse(gz, sprintf("%s/%s_1.fq.gz", diro2, SampleID), 
                        sprintf("%s/%s_1.fq", diro1, SampleID)),
           nR2 = ifelse(gz, sprintf("%s/%s_2.fq.gz", diro2, SampleID), 
                        sprintf("%s/%s_2.fq", diro1, SampleID))) %>%
    mutate(cmd1 = sprintf("ln -sf %s %s", R1, nR1),
           cmd2 = sprintf("ln -sf %s %s", R2, nR2))

map_int(ti$cmd1, system)
map_int(ti$cmd2, system)

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
th = ti %>% select(SampleID, Tissue, Genotype) %>%
    mutate(Tissue = ifelse(Tissue %in% names(tismap), tismap[Tissue], Tissue)) %>%
    separate(Tissue, c("Tissue", "Treatment"), fill='right', extra='merge') %>%
    mutate(Replicate = '', paired = T) %>%
    select(SampleID, Tissue, Genotype, Treatment, everything())
th %>% count(Tissue, Treatment) %>% print(n=23)
th = sra_fill_replicate(th)
fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)
#}}}

#{{{ biomap
sid = 'me99c'
diri = '~/projects/biomap/analysis/01_exp_design'
fi = file.path(diri, "04.typo.corrected.tsv")
ti = read_tsv(fi) %>%
    mutate(raw_idx = sid,
           sid = sprintf("bm%03d", idx)) %>%
    select(sid, raw_idx, sample, everything()) %>%
    select(-idx)
fo = file.path(diri, "06.sample.renamed.tsv")
#write_tsv(ti, fo)

# create sym-links, write read list
th = ti %>%
    transmute(SampleID = sid, Tissue = tissue, Genotype = genotype, 
              Treatment = NA, Replicate = NA,
              paired = T,
              inbred = ifelse(type == 'Inbred', T, F), 
              fi = Readfile,
              fv = file.exists(Readfile))
table(th$fv)

gts = th %>% filter(inbred) %>% distinct(Genotype) %>% pull(Genotype)
pas = strsplit(th$Genotype[!th$inbred], split = 'x')
pas = unique(unlist(pas))
pas[! pas %in% gts]

diro = sprintf("%s/cache/%s/09_fq_interleaved", dird, sid)
if(!dir.exists(diro)) system(sprintf("mkdir -p %s", diro))
cmds = th %>% 
    mutate(fl = sprintf("%s/%s.fq.gz", diro, SampleID)) %>%
    mutate(cmd = sprintf("ln -sf %s %s", fi, fl)) %>%
    pull(cmd)
map_int(cmds, system)

tp = th %>% select(-fi, -fv) %>%
    separate(Tissue, c('Tissue', 'rep'), by = '_')
tp %>% filter(!is.na(rep))
tp %>% filter(!is.na(rep)) %>% distinct(Tissue, Genotype)
tp = tp %>% mutate(Replicate = as.integer(ifelse(is.na(rep), 1, 2))) %>%
    select(-rep)

fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(tp, fo, na = '')
#}}}

#{{{ sna
#}}}

#{{{ Enders
sid = 'me99d'
fi = sprintf("%s/05_read_list/%s.raw.tsv", dird, sid)
ti = read_tsv(fi)
#}}}


