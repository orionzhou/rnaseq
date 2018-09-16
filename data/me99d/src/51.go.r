require(tidyverse)
source("/home/springer/zhoux379/git/luffy/r/enrich.R")
dirp = '/home/springer/zhoux379/projects/3rnaseq'
dirw = file.path(dirp, "data/51.go.enrich")

# read gene IDs
fi = file.path(dirw, "gene_IDs_TAE_061118.csv")
ti = read_csv(fi)
gids_cold_d = ti %>% filter(!is.na(down_cold)) %>% pull(down_cold)
gids_cold_u = ti %>% filter(!is.na(up_cold)) %>% pull(up_cold)
gids_heat_d = ti %>% filter(!is.na(down_heat)) %>% pull(down_heat)
gids_heat_u = ti %>% filter(!is.na(up_heat)) %>% pull(up_heat)

# run GO enrichment test 
tgo.cd = go_enrich_all(gids_cold_d)
tgo.cu = go_enrich_all(gids_cold_u)
tgo.td = go_enrich_all(gids_heat_d)
tgo.tu = go_enrich_all(gids_heat_u)
fo = file.path(dirw, "10.cold.down.tsv")
write_tsv(tgo.cd, fo)
fo = file.path(dirw, "10.cold.up.tsv")
write_tsv(tgo.cu, fo)
fo = file.path(dirw, "10.heat.down.tsv")
write_tsv(tgo.td, fo)
fo = file.path(dirw, "10.heat.up.tsv")
write_tsv(tgo.tu, fo)


