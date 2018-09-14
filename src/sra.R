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


