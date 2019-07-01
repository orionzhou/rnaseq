#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Normalize raw read count matrix to CPM and FPKM by DESeq2 and edgeR')
parser$add_argument("sam", nargs=1, help="sample list table (*.tsv)")
parser$add_argument("exp", nargs=1, help="raw read count matrix (*.rds)")
parser$add_argument("out", nargs=1, help="output file (*.rds)")
parser$add_argument("--opt", default='mmquant', help="data matrix format [default: %(default)s]")
parser$add_argument("--yid", default='mex', help="study ID [default: %(default)s]")
parser$add_argument("--config", default='none',
                    help="genome configuration file (*.rds) [default: %(default)s]")
args <- parser$parse_args()

f_sam = args$sam
f_exp = args$exp
f_out = args$out
opt = args$opt
yid = args$yid
f_cfg = args$config

if( file.access(f_exp) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_exp))
if( file.access(f_sam) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_sam))
size.gene = F

require(devtools)
load_all("~/git/rmaize")

if( f_cfg != 'none')
    size.gene = readRDS(f_cfg)$gene %>% select(gid, size=size.exon)

th = read_tsv(f_sam)

if(opt == 'mmquant') {
    t_rc = readRDS(f_exp) %>%
        rename(SampleID = 1, gid = 2, ReadCount = 3) %>%
        filter(SampleID %in% th$SampleID) %>%
        filter(!str_detect(gid, '--'))
} else if (opt == 'featurecounts') {
    t_rc = readRDS(f_exp) %>%
        rename(SampleID = 1, gid = 2, ReadCount = 3) %>%
        filter(SampleID %in% th$SampleID)
}

res = readcount_norm(t_rc, size.gene)

ths = th %>% distinct(Tissue, Genotype, Treatment) %>%
    mutate(nSampleID = sprintf("%s_%d", yid, 1:length(Tissue)))
t_map = th %>% inner_join(ths, by = c("Tissue", "Genotype", "Treatment")) %>%
    select(SampleID, nSampleID)
th_m = ths %>% select(SampleID=nSampleID, Tissue, Genotype, Treatment)

t_rc_m = t_rc %>% inner_join(t_map, by = 'SampleID') %>%
    mutate(SampleID = nSampleID) %>%
    group_by(gid, SampleID) %>%
    summarise(ReadCount = sum(ReadCount)) %>%
    ungroup()
resm = readcount_norm(t_rc_m, size.gene)

res$th = th
res$th_m = th_m
res$tm_m = resm$tm
saveRDS(res, file = f_out)

