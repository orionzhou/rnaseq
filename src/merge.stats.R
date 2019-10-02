#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser(description = 'merge stats tables')
parser$add_argument("fi", nargs='+', help = "tabular stats file(s)")
parser$add_argument("-o", dest = 'fo', metavar = 'output',
                    nargs=1, default="stats.tsv",
                    help = "output file [default: %(default)s]")
parser$add_argument("--opt", metavar = 'option',
                    nargs=1, default="bam_stat",
                    help = "stats file format [default: %(default)s]")
args <- parser$parse_args()

fis = args$fi
fo = args$fo
opt = args$opt

require(dplyr)
require(readr)
require(stringr)
require(tidyr)
require(purrr)

nfile = length(fis)
ti = tibble(fi = fis) %>%
    mutate(fname = map_chr(fi, basename)) %>%
    mutate(sid = str_replace(fname, '[\\.]gz$', '')) %>%
    mutate(sid = str_replace(sid, '[\\.]\\w+$', ''))

read_featurecount <- function(fi)
    read_tsv(fi, skip = 1) %>%
        select(gid = 1, cnt = 7)

read_mmquant <- function(fi)
    read_tsv(fi, col_types = 'ci') %>%
        select(gid = 1, cnt = 2)

if (opt == 'bam_stat') {
    to = ti %>%
        mutate(data = map(fi, read_tsv, col_names=c('type','cnt'), col_types='cd')) %>%
        select(sid, data) %>%
        unnest() %>%
        spread(type, cnt)
    write_tsv(to, fo)
} else if (opt == 'featurecounts') {
    to = ti %>%
        mutate(data = map(fi, read_featurecount)) %>%
        select(sid, data) %>%
        unnest()
    saveRDS(to, file=fo)
} else if (opt == 'mmquant') {
    to = ti %>%
        mutate(data = map(fi, read_mmquant)) %>%
        select(sid, data) %>%
        unnest()
    saveRDS(to, file=fo)
} else if (opt == 'ase') {
    to = ti %>%
        mutate(data = map(fi, read_featurecount)) %>%
        select(sid, data) %>%
        unnest() %>%
        separate(sid, c('sid','allele'), sep='[\\.]', extra='merge') %>%
        mutate(allele=str_c("allele", allele, sep='')) %>%
        spread(allele, cnt)
    saveRDS(to, file=fo)
} else if (opt == 'ase_snp') {
    to = ti %>%
        mutate(data = map(fi, read_tsv)) %>%
        select(sid, data) %>%
        unnest() %>%
        mutate(allele1 = ifelse(gt=='1|0', altsupport, refsupport)) %>%
        mutate(allele2 = ifelse(gt=='1|0', refsupport, altsupport)) %>%
        select(sid,chr,pos,ref,alt,gt,allele1,allele2)
    saveRDS(to, file=fo)
} else {
    stop("unknown option", opt, "\n")
}
