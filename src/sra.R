read_sra_run <- function(fi, fi2) {
    #{{{
    ti0 = read_csv(fi)
    ti = ti0 %>% select(Run, spots, spots_with_mates, avgLength, 
                        LibraryName, LibraryLayout, SampleName, 
                        BioSample, Sample, Experiment) %>%
        mutate(paired = ifelse(spots_with_mates/spots >= .5, T, F)) %>%
        print(width = Inf)
    ti %>% count(LibraryLayout, paired) %>% print(n=5)
    if(file.exists(fi2)) {
        ti2 = read_csv(fi2) %>%
            transmute(Experiment = `Experiment Accession`, 
                      Title = `Experiment Title`)
        ti = ti %>% left_join(ti2, by = 'Experiment') %>% select(-Experiment)
    }
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
    th = th %>% mutate(Replicate = as.integer(Replicate))
    th %>% dplyr::count(Replicate) %>% print(n=10)
    th
    #}}}
}


