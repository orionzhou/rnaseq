#{{{
source("me.fun.r")
source("sra.R")
t_cfg
#}}}

get_read_list <- function(ti) {
#{{{
if(sid == 'me10a') {
#{{{ me10a - li2010
    th = ti %>%
        mutate(Tissue = sprintf("leaf_%s", SampleName)) %>%
        transmute(SampleID = Run,
                  Tissue = Tissue,
                  Genotype = 'B73',
                  Treatment = '',
                  Replicate = '',
                  paired = paired) %>% 
        arrange(SampleID)
#}}}
} else if(sid == 'me11a') {
#{{{ me11a - Davidson2011
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
        transmute(SampleID = Run,
                  Tissue = Tissue,
                  Genotype = 'B73',
                  Treatment = '',
                  Replicate = '',
                  paired = paired) %>% 
        arrange(SampleID)
#}}}
} else if(sid == 'me12a') {
#{{{ Bolduc2012
    th = ti %>%
        separate(LibraryName, c("gsm1", "gsm"), sep = ": ") %>%
        separate(gsm, c("gentisrep", 'zm', 'rna'), sep = "; ") %>% 
        separate(gentisrep, c("gentis", 'rep'), sep = " #") %>%
        mutate(gentis = str_replace(gentis, " leaf (homo|het)", "_\\1 leaf")) %>%
        separate(gentis, c("gen", "tis"), sep = " ") %>%
        transmute(SampleID = Run,
                  Tissue = tis,
                  Genotype = 'B73',
                  Treatment = gen,
                  Replicate = rep,
                  paired = paired) %>% 
        arrange(SampleID)
#}}}
} else if(sid == 'me13a') {
#{{{ me13a - li2013
    th = ti %>% separate("SampleName", c('org', 'ibm', 'parent', 'tis1', 'tis2', 'Genotype'), sep = "_", fill = 'left')
    th %>% count(parent)
    th %>% count(tis1)
    th %>% count(tis2)
    th %>% count(Genotype)
    th = th %>% transmute(SampleID = Run,
                          Tissue = sprintf("%s_%s", tis1, tis2),
                          Genotype = Genotype,
                          Treatment = '',
                          Replicate = '',
                          paired = paired) %>% 
        arrange(SampleID)
#}}}
} else if (sid == 'me13b') {
#{{{ me13b - liu2013
    th = ti %>% 
        separate("SampleName", c("pre", "Treatment"), sep = "_", fill = "left") %>%
        mutate(Treatment = ifelse(is.na(pre), Treatment, sprintf("E%s", Treatment))) 
    th %>% count(paired)
    th %>% count(Treatment)
    th = th %>% transmute(SampleID = Run,
                          Tissue = 'Leaf',
                          Genotype = 'B73',
                          Treatment = Treatment,
                          Replicate = '',
                          paired = paired) %>% 
        arrange(SampleID)
#}}}
} else if (sid == 'me14a') {
#{{{ me14a - hirsch2014
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
#{{{ me14b - li2014 endosperm
th = ti %>% 
    separate(Title, c("pre", "Treatment"), sep = " B73 ") %>%
    select(-pre) %>%
    separate(Treatment, c("Treatment", "suf"), sep = "DAP; ") %>%
    select(-suf) %>%
    mutate(Treatment = ifelse(Treatment %in% c("0a", "0b"), "0", Treatment)) %>%
    mutate(Treatment = as.integer(Treatment)) %>%
    transmute(SampleID = Run, 
              Tissue = "Endosperm",
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
    transmute(SampleID = Run, 
              Tissue = tis,
              Genotype = 'B73',
              Treatment = '', 
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
    mutate(tis = sprintf("ligule_%s", tis)) %>%
    transmute(SampleID = Run, 
              Tissue = tis,
              Genotype = 'B73',
              Treatment = '', 
              Replicate = '',
              paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me15a') {
#{{{ me15a - leiboff2015
th = ti %>% transmute(SampleID = Run,
                      Tissue = 'SAM',
                      Genotype = SampleName,
                      Treatment = '',
                      Replicate = '',
                      paired = paired) %>% 
    arrange(SampleID)
#}}}
} else if (sid == 'me15b') {
#{{{ me15b - yu2015
    th = ti %>% separate(Title, c("pre", "Treatment"), sep = " at ") %>%
        transmute(SampleID = Run, 
                  Tissue = "Leaf",
                  Genotype = 'B73',
                  Treatment = Treatment, 
                  Replicate = 1,
                  paired = paired) %>%
        filter(paired) %>%
        arrange(SampleID)
#}}}
} else if (sid == 'me16a') {
#{{{ me16a - jin2016
th = ti %>% transmute(SampleID = Run,
                      Tissue = 'kernel',
                      Genotype = SampleName,
                      Treatment = '',
                      Replicate = '',
                      paired = paired) %>% 
    arrange(SampleID)
#}}}
} else if (sid == 'me16b') {
#{{{ me16b - stelpflug2016
tismap = c(
    "Eighth Leaf" = 'leaf_8',
    "Eleventh Leaf" = 'leaf_11',
    "Embryo" = 'embryo',
    "Endopsperm" = 'endosperm',
    "Endosperm" = 'endosperm',
    "Germinating Seed" = 'seed_germ',
    "Immature Leaves" = 'leaf_immature',
    "Primary Root" = 'root_primary',
    "Stem and SAM" = 'stem_SAM',
    "Thirteenth Leaf" = 'leaf_13',
    "Tip of stage-2 Leaf" = 'leaf_tip',
    "Whole seed" = 'seed')
th = ti %>% 
    mutate(Title = str_replace(Title, '_GH_', '_')) %>%
    separate(Title, c('pre', 'str1'), sep = ", B73 ") %>%
    separate(str1, c('tis.age.rep', 'suf'), sep = " RNA-Seq") %>%
    separate(tis.age.rep, c('age','tis','rep'), sep = "_") %>%
    mutate(tis = tismap[tis])
th %>% count(tis)
th %>% count(age)
th %>% count(rep)
th = th %>% transmute(SampleID = Run,
                      Tissue = sprintf("%s_%s", tis, age),
                      Genotype = 'B73',
                      Treatment = '',
                      Replicate = rep,
                      paired = paired) %>%
    arrange(SampleID)
#}}}
} else if (sid == 'me16c') {
#{{{ me16c - walley2016
tismap = c(
"2-4 mm from tip of ear primordium" = 'ear_2-4',
"6-8 mm from tip of ear primordium" = 'ear_6-8', 
"Cortex" = 'root_cortex',
"EMBRYOS" = 'embryo', 
"embryos_20DAP" = 'embryo_20DAP',
"endosperm_12DAP" = 'endosperm_12DAP',
"endosperm_crown" = 'endosperm_crown',
"EZ" = 'root_ez',
"Germinating Kernels" = 'kernel_germ',
"GROWTH ZONE" = 'leaf_growth',
"Internode 6-7" = 'internode_6-7',
"Internode 7-8" = 'internode_7-8',
"mature female spikelets" = 'spikelet',
"MATURE LEAF TISSUE (leaf 8)" = 'leaf_mature_8',
"Mature pollen" = 'pollen',
"MZ" = 'root_mz',
"pericarp_aleurone" = 'pericarp',
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
    transmute(SampleID = Run,
              Tissue = Tissue,
              Genotype = 'B73',
              Treatment = '',
              Replicate = Replicate,
              paired = paired) %>% 
    arrange(SampleID)
th %>% count(Tissue) %>% print(n=23)
#}}}
} else if (sid == 'me17a') {
#{{{ me17a - lin2017
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
#{{{ me18a - kremling2018
th1 = ti %>% separate("LibraryName", c('lib1', 'lib2', 'tissue', 'genotype', 'suf'), sep = "_", fill = 'left')
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
th %>% count(Tissue)
#}}}
} else if (sid == 'me18b') {
#{{{ me18b - baldauf2018
th = ti %>% 
    mutate(SampleName = str_replace(SampleName, '-', '_0_')) %>%
    separate(SampleName, c("gt", "stage", "rep"), sep = "_") %>%
    transmute(SampleID = Run,
              Tissue = sprintf("root %s", stage),
              Genotype = gt,
              Treatment = '',
              Replicate = rep,
              paired = paired) %>% 
    arrange(SampleID)
#}}}
} else {
    cat("unknown study: ", sid, "\n")
}
th = sra_fill_replicate(th)
th
#}}}
}

sid = 'me14d'
fi = sprintf("%s/03_sra_list/%s.csv", dird, sid)
fi2 = sprintf("%s/03_sra_list/%s_exp.csv", dird, sid)
ti = read_sra_run(fi, fi2)

th = get_read_list(ti)
th %>% count(Tissue); th %>% count(Genotype); th %>% count(Replicate)
fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)
create_cache_dir(sid, dird, dirc)


