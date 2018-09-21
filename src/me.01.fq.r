#{{{
source("me.fun.r")
source("sra.R")
t_cfg
#}}}

sid = 'me13b'
fi = sprintf("%s/03_sra_list/%s.csv", dird, sid)
fi2 = sprintf("%s/03_sra_list/%s_exp.csv", dird, sid)
ti = read_sra_run(fi, fi2)

if(sid == 'me13a') {
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
} else if(sid == 'me13b') {
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
} else if(sid == 'me15b') {
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
} else if(sid == 'me14a') {
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
    select(Experiment, Treatment) %>%
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
th = ti %>% 
    mutate(Title = str_replace(Title, '_GH_', '_')) %>%
    separate(Title, c('pre', 'str1'), sep = ", B73 ") %>%
    separate(str1, c('tis.age.rep', 'suf'), sep = " RNA-Seq") %>%
    separate(tis.age.rep, c('age','tis','rep'), sep = "_") 
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
th = ti %>% 
    transmute(SampleID = Run,
              Tissue = Title,
              Genotype = Title,
              Treatment = '',
              Replicate = '',
              paired = paired) %>% 
    arrange(SampleID)
#}}}
} else if (sid == 'me17a') {
#{{{ me17a - lin2017
th1 = ti1 %>% 
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
th2 = ti2 %>% 
    transmute(SampleID = Run,
              Tissue = 'SAM',
              Genotype = LibraryName,
              Treatment = '',
              Replicate = '',
              paired = paired) %>% 
    arrange(SampleID)
th = rbind(th1,th2)
th %>% count(Replicate)
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
th3 %>% tissue
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

fo = sprintf("%s/05_read_list/%s.tsv", dird, sid)
write_tsv(th, fo)
create_cache_dir(sid, dird, dirc)
