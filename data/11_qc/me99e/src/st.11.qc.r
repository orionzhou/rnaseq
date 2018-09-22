require(tidyverse)
dirp = '~/projects/settles'
dirw = file.path(dirp, 'data')
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

#{{{ create sym-links and 01.reads.tsv
fi = file.path(dirw, "00.samplelist.tsv")
ti = read_tsv(fi)
th = ti %>% transmute(SampleID = sid,
                      Tissue = tissue,
                      Genotype = genotype,
                      Treatment = note,
                      Replicate = '',
                      paired = !is.na(Readfile2))
th = sra_fill_replicate(th)

for (i in 1:nrow(ti)) {
    sid = ti$sid[i]; r1 = ti$Readfile1[i]; r2 = ti$Readfile2[i]
    if(is.na(r2) | r2 == '') {
        ff1 = sprintf("%s/cache/09_raw_fastq/%s.fq.gz", dirp, sid)
        ft1 = sprintf("%s/cache/07_raw_fq/%s", dirp, r1)
        cmd1 = sprintf("ln -sf %s %s", ft1, ff1)
        system(cmd1)
    } else{
        ff1 = sprintf("%s/cache/09_raw_fastq/%s_1.fq.gz", dirp, sid)
        ff2 = sprintf("%s/cache/09_raw_fastq/%s_2.fq.gz", dirp, sid)
        ft1 = sprintf("%s/cache/07_raw_fq/%s", dirp, r1)
        ft2 = sprintf("%s/cache/07_raw_fq/%s", dirp, r2)
        cmd1 = sprintf("ln -sf %s %s", ft1, ff1)
        cmd2 = sprintf("ln -sf %s %s", ft2, ff2)
        system(cmd1)
        system(cmd2)
    }
}

fo = file.path(dirw, "01.reads.tsv")
write_tsv(th, fo)
#}}}

#{{{ collect featurecounts data
fi = file.path(dirp, 'cache/33.featurecounts/01.txt')
ti1 = read_tsv(fi, skip = 1)
ti = ti1 %>% bind_cols(ti2[,-c(1:6)])
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
tcl = tcw %>% gather(sid, RawReadCount, -gid)
fo = file.path(dirw, "10.RawReadCount.RData")
save(tcw, tcl, file = fo)
#}}}

#{{{ collect ASE data
t_ase = tibble()
for (i in 1:nrow(th)) {
    sid = th$sid[i]
    gt = th$genotype[i]
    fi1 = sprintf("%s/cache1/26.ase/%s_%s.tsv", dirw, sid, gt)
    fi2 = sprintf("%s/cache2/26.ase/%s_%s.tsv", dirw, sid, gt)
    fi = ifelse(file.exists(fi1), fi1, fi2)
    ti = read_tsv(fi)
    ti2 = ti %>% mutate(sid = sid) %>%
        select(sid, gid, nref, nalt, ncft)
    t_ase = rbind(t_ase, ti2)
    cat(i, sid, "\n")
}

fo = file.path(dirw, "11.ase.RData")
save(t_ase, file = fo)
#}}}

#{{{ first look at ASE
fi = file.path(dirw, "11.ase.RData")
x = load(fi)

tm = th %>% inner_join(t_ase, by = 'sid') %>%
    filter(nref + nalt >= 20) %>%
    mutate(prop.cft = ncft/(nref+nalt+ncft),
           prop.ref = nref/(nref+nalt))

tp = tm %>% 
    group_by(sid, type, tissue, genotype, note) %>%
    summarise(n.pass = sum(prop.cft < 0.05),
              p.fail = sum(prop.cft >= 0.05)/n(),
              pref.q25 = quantile(prop.ref, 0.25),
              pref.q50 = quantile(prop.ref, 0.5),
              pref.q75 = quantile(prop.ref, 0.75)) %>% 
    ungroup()

fo = file.path(dirw, "53.ase.tsv")
write_tsv(tp, fo, na = '')
#}}}
