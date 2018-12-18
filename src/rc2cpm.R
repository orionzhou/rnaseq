#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'Normalize raw read count matrix to CPM and FPKM by DESeq2 and edgeR')
parser$add_argument("exp", nargs=1, help="raw read count table (*.tsv)")
parser$add_argument("out", nargs=1, help="output file (*.rds)")
parser$add_argument("--sample", default='none',
                    help="sample list table (*.tsv) [default: %(default)s]")
parser$add_argument("--config", default='none',
                    help="genome configuration file (*.rds) [default: %(default)s]")
args <- parser$parse_args()

f_exp = args$exp
f_out = args$out
fh = args$sample
f_cfg = args$config

if( file.access(f_exp) == -1 )
    stop(sprintf("file ( %s ) cannot be accessed", f_exp))
size.gene = F
if( f_cfg != 'none')
    size.gene = readRDS(f_cfg)$size.gene

source("~/projects/rnaseq/src/functions.R")

t_rc = read_tsv(f_exp) %>% gather(SampleID, ReadCount, -gid)
if( fh != 'none') {
    th = read_tsv(fh)
    t_rc = t_rc %>% filter(SampleID %in% th$SampleID)
}
res = readcount_norm(t_rc, size.gene)

saveRDS(res, file = f_out)

