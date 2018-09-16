create_cache_dir <- function(study, dirp, dirc) {
    #{{{
    dirw = file.path(dirp, study)
    dirc1 = file.path(dirc, study)
    cmd = sprintf("mkdir -p %s", dirc1)
    system(cmd)
    if(file.exists(file.path(dirw, 'cache'))) system(sprintf("rm %s/cache", dirw))
    cmd = sprintf("ln -sf %s/ %s/cache", dirc1, dirw)
    system(cmd)
    fread = file.path(dirw, 'data/01.reads.tsv')
    stopifnot(file.exists(fread))
    cmd = sprintf("ln -sf %s %s/01.reads.tsv", fread, dirc1)
    system(cmd)
    cmd = sprintf("mkdir -p %s/data/raw", dirw)
    system(cmd)
    cmd = sprintf("ln -sf %s/data/raw/ %s/data", dirw, dirc1)
    system(cmd)
    #}}}
}

