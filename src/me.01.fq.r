#{{{
source("me.fun.r")
source("sra.R")
source("snk.R")
tl
#}}}

#{{{ me13a - li2013
study = 'me13a'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

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
th = sra_fill_replicate(th)
th %>% count(Replicate)
#}}}

#{{{ me13b - liu2013
study = 'me13b'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

th = ti %>% separate("SampleName", c("pre", "Treatment"), sep = "_", fill = "left")
th %>% count(paired)
th %>% count(Treatment)
th = th %>% transmute(SampleID = Run,
                      Tissue = 'Leaf',
                      Genotype = 'B73',
                      Treatment = Treatment,
                      Replicate = '',
                      paired = paired) %>% 
    arrange(SampleID)
th = sra_fill_replicate(th)
th %>% count(Replicate)
#}}}

#{{{ me15b - yu2015
study = 'me15b'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

fm = file.path(dirw, "sra_result.csv")
tm = read_csv(fm)
tm = tm %>% separate(`Experiment Title`, c("pre", "Treatment"), sep = " at ") %>%
    transmute(Experiment = `Experiment Accession`, Treatment = Treatment)

th = ti %>% inner_join(tm, by = 'Experiment') %>% 
    transmute(SampleID = Run, 
              Tissue = "Leaf",
              Genotype = 'B73',
              Treatment = Treatment, 
              Replicate = 1,
              paired = paired) %>%
    filter(paired) %>%
    arrange(SampleID)
#}}}

#{{{ me14a - hirsch2014
study = 'me14a'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

fm = file.path(dirw, "biosample_result.tsv")
tm = read_tsv(fm)
tm = tm %>% separate(Title, c("org", "Genotype"), sep = ", ") %>%
    separate(Genotype, c("Genotype", "suf"), sep = " RNAseq") %>%
    select(-org, -suf)

th = ti %>% inner_join(tm, by = 'BioSample') %>% 
    transmute(SampleID = Run, 
              Tissue = "seedling",
              Genotype = Genotype,
              Treatment = '', 
              Replicate = 1,
              paired = paired) %>%
    arrange(SampleID)

th = sra_fill_replicate(th)
#}}}

#{{{ me14b - li2014 endosperm
study = 'me14b'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

fm = file.path(dirw, "sra_result.csv")
tm = read_csv(fm)
tm = tm %>% mutate(Experiment = `Experiment Accession`) %>%
    separate(`Experiment Title`, c("pre", "Treatment"), sep = " B73 ") %>%
    select(Experiment, Treatment) %>%
    separate(Treatment, c("Treatment", "suf"), sep = "DAP; ") %>%
    select(-suf) %>%
    mutate(Treatment = ifelse(Treatment %in% c("0a", "0b"), "0", Treatment)) %>%
    mutate(Treatment = as.integer(Treatment)) 

th = ti %>% inner_join(tm, by = 'Experiment') %>% 
    transmute(SampleID = Run, 
              Tissue = "Endosperm",
              Genotype = 'B73',
              Treatment = Treatment, 
              Replicate = '',
              paired = paired) %>%
    arrange(SampleID)
th = sra_fill_replicate(th)
th
#}}}

#{{{ me15a - leiboff2015
study = 'me15a'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

th = ti %>% transmute(SampleID = Run,
                      Tissue = 'SAM',
                      Genotype = SampleName,
                      Treatment = '',
                      Replicate = '',
                      paired = paired) %>% 
    arrange(SampleID)
th = sra_fill_replicate(th)
th %>% count(Replicate)
#}}}

#{{{ me16a - jin2016
study = 'me16a'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

th = ti %>% transmute(SampleID = Run,
                      Tissue = 'kernel',
                      Genotype = SampleName,
                      Treatment = '',
                      Replicate = '',
                      paired = paired) %>% 
    arrange(SampleID)
th = sra_fill_replicate(th)
th %>% count(Replicate)
#}}}

#{{{ me16b - stelpflug2016
study = 'me16b'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

fx = file.path(dirw, "sra_result.csv")
tx = read_csv(fx) %>% transmute(Experiment = `Experiment Accession`,
                                exp.title = `Experiment Title`)

th = ti %>% inner_join(tx, by = 'Experiment') %>%
    mutate(exp.title = str_replace(exp.title, '_GH_', '_')) %>%
    separate(exp.title, c('pre', 'str1'), sep = ", B73 ") %>%
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

#{{{ me16c - walley2016
study = 'me16c'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

fm = file.path(dirw, "biosample_result.tsv")
tm = read_tsv(fm)
tm = tm %>% select(BioSample, genotype, tissue)

th = ti %>% left_join(tm, by = 'BioSample') %>%
    transmute(SampleID = Run,
              Tissue = tissue,
              Genotype = genotype,
              Treatment = '',
              Replicate = '',
              paired = paired) %>% 
    arrange(SampleID)
th = sra_fill_replicate(th)
th %>% count(Replicate)
#}}}

#{{{ me17a - lin2017
study = 'me17a'
dirw = file.path(dirp, study, 'data') 

fi1 = file.path(dirw, "SraRunInfo1.csv")
ti1 = read_sra_run(fi1)
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

fi2 = file.path(dirw, "SraRunInfo2.csv")
ti2 = read_sra_run(fi2)
th2 = ti2 %>% 
    transmute(SampleID = Run,
              Tissue = 'SAM',
              Genotype = LibraryName,
              Treatment = '',
              Replicate = '',
              paired = paired) %>% 
    arrange(SampleID)
th2 %>% count(Genotype)

th = sra_fill_replicate(rbind(th1,th2))
th %>% count(Replicate)
#}}}

#{{{ me18a - kremling2018
study = 'me18a'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

th1 = ti %>% separate("LibraryName", c('lib1', 'lib2', 'tissue', 'genotype', 'suf'), sep = "_", fill = 'left')
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
    transmute(SampleID = Run,
              Tissue = tissue,
              Genotype = genotype,
              paired = paired) %>%
    arrange(SampleID)

th = th4 %>% mutate(Treatment = '', Replicate = '') %>%
    select(SampleID, Tissue, Genotype, Treatment, Replicate, paired)
th %>% count(Tissue)
th = sra_fill_replicate(th)
th %>% count(Replicate)
#}}}

#{{{ me18b - baldauf2018
study = 'me18b'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

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
th %>% count(Replicate)
#}}}

#{{{ me99a - kaeppler
study = 'me99a'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)
#}}}

#{{{
#}}}

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp, dirc)
