require(tidyverse)
dirw = '~/projects/kremling2018'

#{{{ read sra table
fi = file.path(dirw, "data/SraRunInfo.csv")
ti = read_csv(fi) %>% select(Run, LibraryName)
colnames(ti) = c("sid", "libname")

th1 = ti %>% separate("libname", c('lib1', 'lib2', 'tissue', 'genotype', 'suf'), sep = "_", fill = 'left')
table(th1$tissue)
tissues = "GRoot GShoot Kern L3Base L3Mid L3Tip LMAD26 LMAD8 LMAN26 LMAN8 LMid"
tissues = strsplit(tissues, " ")[[1]]
idxs = (!th1$tissue %in% tissues)

th2 = th1 %>% 
    mutate(
           genotype = ifelse(idxs, sprintf("%s_%s", tissue, genotype), genotype),
           tissue = ifelse(idxs, lib2, tissue),
           lib2 = ifelse(idxs, '', lib2))
table(th2$tissue)

th3 = th2 %>% 
    mutate(tissue = ifelse(tissue %in% c("LMAD26", "LMAD8"), "LMAD", tissue)) %>%
    mutate(tissue = ifelse(tissue %in% c("LMAN26", "LMAN8"), "LMAN", tissue))
table(th3$tissue)

tissues = "GRoot GShoot Kern L3Base L3Tip LMAD LMAN" 
tissues = strsplit(tissues, " ")[[1]]
th4 = th3 %>% filter(tissue %in% tissues) %>% 
    transmute(sid = sid, Tissue = tissue, Genotype = genotype) %>%
    arrange(sid)
table(th4$Tissue)

to = th4 %>% transmute(SampleID = sid, Tissue = Tissue, Genotype = Genotype, 
                       Treatment = 1)
fo = file.path(dirw, "data/01.reads.tsv")
write_tsv(to, fo)
#}}}

