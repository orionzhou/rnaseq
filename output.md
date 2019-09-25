### Sample list / meta table: `meta.tsv`
- `SampleID`, `Tissue`, `Gentoype`, `Treatment`, `Replicate`
- `paired`: paired end or not
- `spots`: number of reads (single-end) or pairs (paired-end)
- `avgLength`: average read length

### Trimming statitics: `trimming.tsv`
- `sid`: SampleID
- `passed_filter_reads`, `low_quality_reads`, `too_many_N_reads`, `too_short_reads`, `too_long_reads`

### Mapping statistics: `bamstats.tsv`
- `sid`: SampleID
- `pair`: read pairs
  - `pair_bad`, `pair_dup`: pairs that failed QC or duplicates
  - `pair_map`: mapped pairs (both ends)
  - `pair_orphan`: pairs with one end mapped
  - `pair_unmap`: unmapped pairs
- `unpair`: singletons (single-end reads or pairs with one end failing QC)
  - `unpair_bad`, `unpair_dup`: singletons that failed QC or duplicates
  - `unpair_map`: mapped reads
  - `unpair_unmap`: unmapped reads
- `pair_map_hq`, `pair_orphan_hq`, `unpair_map_hq`: pairs/reads mapped
  with high quality (i.e., unique)
- `pair_map0`, `pair_orphan0`, `unpair_map0`: pairs/reads mapped with 0 mismatch
- `pair_map_hq0`, `pair_orphan_hq0`, `unpair_map_hq0`: pairs/reads mapped 
  with high quality (i.e., unique) and with 0 mismatch

### Raw read count and normalized CPM / FPKM table: `cpm.rds`
- can be loaded into R using `x = readRDS("cpm.rds")`, contains the following data frames:
- `tl` - library stats
  - `SampleID`, `libSize`
  - `sizeFactor`: DESeq2 library size factor
  - `normFactor`: edgeR library normalization factor
- `th` - sample list / meta table, same with `meta.tsv`
- `tm` - expression table
  - `gid`: Gene ID (AGP_v4, Ensembl Plants v37, 46,117 in total)
  - `SampleID`
  - `ReadCount`: raw read count
  - `nRC`: normalized read count (`nRC = ReadCount / sizeFactor`)
  - `rCPM`: raw CPM (adds up to 1,000,000 for each sample/library)
  - `rFPKM`: raw FPKM calculated using rCPM and gene exon length
  - `CPM`: CPM calculated by edgeR (`CPM = rCPM / normFactor`)
  - `FPKM`: FPKM calculated using CPM and gene exon length
- `th_m`: replicate merged sample list / meta table
- `tm_m`: replicate merged expression table

### Gene-based allele specific read counts: `ase.rds`
- can be loaded into R using `x = readRDS("ase.rds")`, contains the following columns:
- `sid`: Sample ID, two for each sample ("Sample.1" and "Sample.2"). For example, in the case of B73xMo17, "Sample.1" represents the B73 (first) allele and "Sample.2" represents the Mo17 (second) allele
- `gid`: gene ID
- `cnt`: allele-specific read count

### SNP-based allele specific read counts: `ase2.rds`
- can be loaded into R using `x = readRDS("ase2.rds")`, contains the following columns:
- `sid`: Sample ID
- `chr`, `pos`, `ref`, `alt`: SNP information
- `cnt_ref`: ref allele-specific read count
- `cnt_alt`: alt allele-specific read count
- `genotype`: sample genotype at this site: in cases of "1|0", allele 1 (maternal allele) is in `alt` state while allele 2 (paternal allele) is the `ref` state
