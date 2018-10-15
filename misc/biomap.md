# biomAP

[This table](/data/05_read_list/me99c.tsv) contains raw sample meta-data, while [this table](/data/05_read_list/me99c.c.tsv) has corrected meta-data (sample `bm252` corrected from `Root` to `Leaf` tissue).
  - msi path: `/home/springer/zhoux379/projects/maize.expression/data/05_read_list/me99c.c.tsv`

Intermediate files (\*.fastq, \*.bam) are under msi scratch space:
  - `/scratch.global/zhoux379/maize.expression/me99c/`, with sub-directories:
  - `10_fastq`: deinterleaved fastq files in gzip format
  - `14_trim`: adapter trimmed fastq files by Trimmomatic
  - `21_star`: raw BAM files mapped by STAR
  - `22_bam`: sorted and indexed BAM files
  - `31_featurecounts`: Read count tables for 46,117 Ensembl Plants v37 genes
  - `33_ase`: intermediate files in the allele count pipeline

[Mapping stats table](/data/11_qc/me99c/10.mapping.stat.tsv) contains trimming, mapping and counting statistics for each sample with the following columns:
  - msi path: `/home/springer/zhoux379/projects/maize.expression/data/11_qc/me99c/10.mapping.stat.tsv`
  - Sample meta data: `SampleID`, `Tissue`, `Genotype`, `Treatment`, `Replicate`: 
  - Trimming statistics: `total`, `surviving`, `surviving_f`, `surviving_r`, `dropped`
  - Mapping statistics:
    - `pair`: read pairs
      - `pair_bad`, `pair_dup`: pairs that failed QC or duplicates
      - `pair_map`: mapped pairs (both ends)
      - `pair_orphan`: pairs with one end mapped
      - `pair_unmap`: unmapped pairs
    - `unpair`: singletons (single-end reads or pairs with one end failing QC)
      - `unpair_bad`, `unpair_dup`: singletons that failed QC or duplicates
      - `unpair_map`: mapped reads
      - `unpair_unmap`: unmapped reads
    - `pair_map_hq`, `pair_orphan_hq`, `unpair_map_hq`: pairs/reads mapped with high quality (i.e., unique)
    - `pair_map0`, `pair_orphan0`, `unpair_map0`: pairs/reads mapped with 0 mismatch
    - `pair_map_hq0`, `pair_orphan_hq0`, `unpair_map_hq0`: pairs/reads mapped with high quality (i.e., unique) and with 0 mismatch
  - Read counting statistics:
    - `Assigned`: reads assigned to exonic regions and thus counted
    - `Unassigned_MultiMapping`, `Unassigned_NoFeatures`, `Unassigned_Ambiguity`, `Unassigned_Unmapped`: reads not counted due to various reasons
    
R data file (msi path: `/home/springer/zhoux379/projects/maize.expression/data/11_qc/me99c/20.rc.norm.rda`) containing raw read count tables and normalized expression values, with the following objects (tibbles):
* tl - tibble for library (sample), with columns:
  * `SampleID`: bm001 - bm467
  * ` Tissue`: Leaf, Internode, Root, etc.
  * `Genotype`: B73, Mo17xPH207, etc
  * `Treatment`: replicate 1 or 2
  * `inbred`: whether this is the inbred parent (TRUE) or hybrid (FALSE)
  * `sizeFactor`, `libSize`: library size and normalization factor calculated using the median log ratio approach by DESeq2, accounts for library size
  * `normFactor`: library normalzation factor computed by edgeR using the TMM approach, does NOT account for library size
* tm - tibble for biomap expression data
  * `gid`: Gene ID (AGP_v4, Ensembl Plants v37, 46,117 in total)
  * `SampleID`: bm001 - bm467
  * `ReadCount`: raw read count
  * `nRC`: normalized read count (`nRC = ReadCount / sizeFactor`)
  * `rCPM`: raw CPM (adds up to 1,000,000 for each sample/library)
  * `rFPKM`: raw FPKM calculated using rCPM and gene exon length
  * `CPM`: CPM calculated by edgeR (`CPM = rCPM / normFactor`)
  * `FPKM`: FPKM calculated using CPM and gene exon length


