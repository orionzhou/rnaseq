source("functions.R")

yid = 'rn20d3'
dirw = file.path(dird, '11_qc', yid)
if(!dir.exists(dirw)) system(sprintf("mkdir -p %s", dirw))

#{{{ raw: read in
res = rnaseq_cpm_raw(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(SampleID, Tissue, Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Genotype=='B73' & Replicate==1, Tissue, '')) %>%
    arrange(Tissue, inbred, Genotype, Treatment)
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ raw: hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    expand.x=.2, lab.size=2)
ggsave(file.path(dirw, '11.hclust.p.pdf'), p1, width=8, height=30)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    expand.x=.2, lab.size=2)
ggsave(file.path(dirw, '11.hclust.s.pdf'), p1, width=8, height=30)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=F, pca.scale=F,
    var.shape='Tissue',var.col='Tissue',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=8,iter=1500, seed=42,
    var.shape='Tissue',var.col='Tissue',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '11.tsne.pdf'), p3, width=8, height=8)
#}}}

#{{{ raw: filter/fix samples
th2 = res$th %>%
    mutate(Tissue = ifelse(SampleID == 'bm318', 'Leaf', Tissue))
th2 = complete_sample_list(th2)

fh = file.path(dirw, '01.meta.tsv')
write_tsv(th2, fh, na='')
#}}}

#{{{ read in
res = rnaseq_cpm(yid)
th = res$th; tm = res$tm; tl = res$tl; th_m = res$th_m; tm_m = res$tm_m

th = res$th %>%
    mutate(lab = str_c(Tissue, Genotype, Replicate, sep='_')) %>%
    mutate(clab = ifelse(Genotype=='B73' & Replicate==1, Tissue, '')) %>%
    arrange(Tissue, inbred, Genotype, Treatment)
tm = res$tm %>% filter(SampleID %in% th$SampleID) %>%
    mutate(value=asinh(CPM))
#}}}

#{{{ hclust & tSNE
p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='pearson',var.col='Tissue',
    expand.x=.2, lab.size=2)
ggsave(file.path(dirw, '21.hclust.p.pdf'), p1, width=8, height=30)

p1 = plot_hclust(tm,th,pct.exp=.7,cor.opt='spearman',var.col='Tissue',
    expand.x=.2, lab.size=2)
ggsave(file.path(dirw, '21.hclust.s.pdf'), p1, width=8, height=30)

p2 = plot_pca(tm,th,pct.exp=.7, pca.center=F, pca.scale=F,
    var.shape='Tissue',var.col='Tissue',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.pca.pdf'), p2, width=8, height=8)

p3 = plot_tsne(tm,th,pct.exp=.7,perp=8,iter=1500, seed=42,
    var.shape='Tissue',var.col='Tissue',var.lab='clab',var.ellipse='Tissue',
    legend.pos='top.left', legend.dir='v', pal.col='aaas')
ggsave(file.path(dirw, '21.tsne.pdf'), p3, width=8, height=8)
#}}}

#{{{ ase
pa1 = plot_ase(res$ase_gene, th, val.col='Tissue', pal.col='aaas')
fo = file.path(dirw, '31.afs_gene.pdf')
ggsave(fo, pa1, width=7, height=35)

pa2 = plot_ase(res$ase_snp, th, val.col='Tissue', pal.col='aaas')
fo = file.path(dirw, '32.afs_site.pdf')
ggsave(fo, pa2, width=7, height=35)
#}}}



#{{{ # make sup-table for PeteC
th = rnaseq_sample_meta(yid)
tt = rnaseq_mapping_stat(yid)
th = th %>% inner_join(tt, by='SampleID') %>%
    select(SampleID, paired, avgLength,
        TotalReads=total, PassedTrimming=passed,
        MappedReads=mapped, UniquelyMappedReads=mappedu)

fi = '~/projects/barn/data/06_local_list/rn99c.tsv'
ti = read_tsv(fi) %>% mutate(fname=basename(r0)) %>%
    mutate(fname=str_replace(fname,'fastq','anqrpt.fastq')) %>%
    select(SampleID,fname)
fp = file.path(dirw, 'Supplemental_table_mRNA_libs_tmp.csv')
tp = read_csv(fp) %>% inner_join(ti, by=c('fastq_file_ID'='fname')) %>%
    inner_join(th, by='SampleID') %>% select(-SampleID)

fo = file.path(dirw, 'sup.csv')
write_csv(tp, fo)

### sanity check
fi = '~/projects/barn/data/06_local_list/rn99c.tsv'
ti = read_tsv(fi) %>% mutate(fname=basename(r0)) %>%
    mutate(fname=str_replace(fname,'fastq','anqrpt.fastq')) %>%
    select(SampleID,tis=Tissue,gt=Genotype,fname)
fp = file.path(dirw, 'Supplemental_table_mRNA_libs.csv')
tp = read_csv(fp) %>% inner_join(ti, by=c('fastq_file_ID'='fname')) %>%
    inner_join(th, by='SampleID')

tx1 = tp %>% distinct(tis, gt)
th0 = rnaseq_sample_meta(yid)
th0 %>%
    inner_join(tx1, by=c("Tissue"='tis','Genotype'='gt')) %>%
    count(Tissue, Genotype) %>% filter(n>1)
#}}}

