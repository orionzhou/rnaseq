require(tidyverse)
dirw = '~/projects/biomap'

#{{{ sample list -> analysis/01.sample.RData
fh = file.path(dirw, "data/01.reads.tsv")
th = read_tsv(fh) %>%
    separate(Genotype, c("pa1", "pa2"), "x", remove = F) %>%
    replace_na(list(pa2=''))

thg = th %>% distinct(inbred, Genotype) %>%
    arrange(desc(inbred), Genotype)
gts = thg %>% pull(Genotype)
gts_i = thg %>% filter(inbred) %>% pull(Genotype)
gts_h = thg %>% filter(!inbred) %>% pull(Genotype)
tissues5 = c("Root", "Internode", "Leaf", "Seedling", "Endosperm")
tl = th %>% select(-pa1, -pa2) %>%
    mutate(Tissue = factor(Tissue, levels = tissues5),
           Genotype = factor(Genotype, levels = gts))
fo = file.path(dirw, 'analysis/01.sample.RData')
save(tl, tissues5, gts, gts_i, gts_h, file = fo)
#}}}

#{{{ collect featurecounts data
fi = file.path(dirw, 'cache/31_featurecounts/01.txt')
ti = read_tsv(fi, skip = 1)
cnames = colnames(ti)[-c(1:6)]
res = strsplit(cnames, split = ":")
sids = sapply(res, "[", 2)
tcw = ti[,-c(2:6)]
colnames(tcw) = c('gid', sids)
#t_rc = tcw %>% gather(sid, RawReadCount, -gid)
t_rc = tcw
fo = file.path(dirw, "analysis/10.RawReadCount.RData")
save(t_rc, file = fo)
#}}}

#{{{ collect ASE
t_ase = tibble()
for (i in 1:nrow(th)) {
    sid = th$sid[i]
    gt = th$genotype[i]
    fi = sprintf("%s/cache/33_ase/%s_%s.tsv", dirw, sid, gt)
    ti = read_tsv(fi)
    ti2 = ti %>% mutate(sid = sid) %>%
        select(sid, everything())
    t_ase = rbind(t_ase, ti2)
    cat(i, sid, "\n")
}
fo = file.path(dirw, "analysis/11.ase.RData")
save(t_ase, file = fo)
#}}}


