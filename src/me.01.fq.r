#{{{
source("me.fun.r")
studies
read_sra_run <- function(fi) {
    #{{{
    ti0 = read_csv(fi)
    ti = ti0 %>% select(Run, spots, spots_with_mates, avgLength, 
                        LibraryName, LibraryLayout, SampleName, 
                        BioSample, Sample, Experiment) %>%
        mutate(paired = ifelse(spots_with_mates/spots >= .5, T, F)) %>%
        print(width = Inf)
    ti %>% count(LibraryLayout, paired) %>% print(n=5)
    ti
    #}}}
}
sra_fill_replicate <- function(th) {
    #{{{
    cmap = c()
    for (i in 1:nrow(th)) {
        tis = th$Tissue[i]; gt = th$Genotype[i]; trt = th$Treatment[i]
        key = sprintf("%s-%s-%s", tis, gt, trt)
        if (key %in% names(cmap))
            cmap[key] = cmap[key] + 1
        else
            cmap[key] = 1
        th$Replicate[i] = cmap[key]
    }
    th %>% count(Replicate) %>% print(n=10)
    th
    #}}}
}
create_cache_dir <- function(study, dirp) {
    #{{{
    dircp = "/scratch.global/zhoux379/maize.expression"
    dirw = file.path(dirp, study)
    dirc = file.path(dircp, study)
    cmd = sprintf("mkdir -p %s", dirc)
    system(cmd)
    if(file.exists(file.path(dirw, 'cache'))) system(sprintf("rm %s/cache", dirw))
    cmd = sprintf("ln -sf %s/ %s/cache", dirc, dirw)
    system(cmd)
    fread = file.path(dirw, 'data/01.reads.tsv')
    stopifnot(file.exists(fread))
    cmd = sprintf("ln -sf %s %s/01.reads.tsv", fread, dirc)
    system(cmd)
    cmd = sprintf("mkdir -p %s/data/raw", dirw)
    system(cmd)
    cmd = sprintf("ln -sf %s/data/raw/ %s/data", dirw, dirc)
    system(cmd)
    #}}}
}
#}}}

#{{{ li2013
study = 'li2013'
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

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ hirsch2014
study = 'hirsch2014'
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
fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ walley2016
study = 'walley2016'
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

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ leiboff2015
study = 'leiboff2015'
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

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ kremling2018
study = 'kremling2018'
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

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ jin2016
study = 'jin2016'
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

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ lin2017
study = 'lin2017'
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

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ stelpflug2016
study = 'stelpflug2016'
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

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ baldauf2018
study = 'baldauf2018'
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

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
create_cache_dir(study, dirp)
#}}}

#{{{ kaeppler2018
study = 'kaeppler2018'
dirw = file.path(dirp, study, 'data') 
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)
#}}}

#{{{
#}}}

#{{{ hapmap3
dirw = '~/projects/hapmap3/data'
fi = file.path(dirw, "SraRunInfo.csv")
ti = read_sra_run(fi)

th = ti %>% 
    separate("SampleName", c('pre', 'Genotype'), sep = "_", fill = 'left') %>%
    transmute(SampleID = Run, 
              Tissue = '', 
              Genotype = Genotype, 
              Treatment = '',
              Replicate = '',
              paired = paired)
th = sra_fill_replicate(th)

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
#}}}

