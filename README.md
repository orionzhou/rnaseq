# Collection of large-scale maize expression datasets

This repo provides a collection of recently published (after 2010) maize gene expression datasets. The focus is on large-scale (multiple tissues / developmental stages, multiple inbred / hybrid lines) Illumina RNA-Seq experiments, but also inlcude experiments done with other sequencing technologies (3' RNA-Seq, MNase-Seq, etc.).

Raw sequencing reads were downloaded from [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra), trimmed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and mapped to the [maize B73 AGP_v4 genome](http://plants.ensembl.org/Zea_mays/Info/Index) using [STAR](https://github.com/alexdobin/STAR).  Uniquely mapped reads were assigned to and counted for the 46,117 reference gene models ([Ensembl Plants v37](http://plants.ensembl.org/Zea_mays/Info/Index)) using [featureCounts](http://bioinf.wehi.edu.au/featureCounts/).  Raw read counts were then normalized using the [TMM normalization approach](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to give CPMs (Counts Per Million reads) and then further normalized by gene CDS lengths to give FPKM (Fragments Per Kilobase of exon per Million reads) values.  Hierarchical clustering and principal component analysis were used to explore sample clustering pattern.

See [this table](/data/01.cfg.tsv) for the list of collected datasets.

| Study | Line/Genotype | Tissue | Reads/Sample | note | results |
| ----- | ------------- | ------ | ------------------- | --- | ---- |
| [Li2013 eQTL](https://www.ncbi.nlm.nih.gov/pubmed/23341782) | 105 RILs | shoot apices | 20-30M | B73 and Mo17 IBM | [link](/li2013/data) |
| [Hirsch2014 pan-transcriptome](https://www.ncbi.nlm.nih.gov/pubmed/24488960 ) | 503 inbreds | seedling | 10-30M | | [link](/hirsch2014/data)|
| [Leiboff2015 SAM](https://www.ncbi.nlm.nih.gov/pubmed/26584889 ) | 380 | SAM | 20-50 | | [link](/leiboff2015/data) |
| [Jin2016 ePAV GWAS](https://www.ncbi.nlm.nih.gov/pubmed/26729541 ) | 368 inbred lines | 15 DAP kernels | 20-50M | phenotypic data | [link](/jin2016/data) |
| [Stelpflug2016 atlas](https://www.ncbi.nlm.nih.gov/pubmed/27898762) | B73 | 18 tissues | 10-20M |  | [link](/stelpflug2016/data) |
| [Walley2016 proteome](https://www.ncbi.nlm.nih.gov/pubmed/27540173 ) | B73 | 23 tissues | 20-30M | proteomics | [link](/walley2016/data) |
| [Lin2017 eRD-GWAS](https://www.ncbi.nlm.nih.gov/pubmed/29041960 ) | 27 inbreds | 5 tissues | 20-30M |  | [link](/lin2017/data) |
| [Kremling2018 dysregulation](https://www.ncbi.nlm.nih.gov/pubmed/29539638 ) | 255 lines | 7 tissues | 2-6M | 3' RNA-Seq | [link](/kremling2018/data) |
| [Baldauf2018 SPE](https://www.ncbi.nlm.nih.gov/pubmed/29358068 ) | 7 inbreds + 6 hybrids | 3 root stages | 20M * 3-4 reps |  | |
| [Wang2018 sorghum pacbio](https://www.ncbi.nlm.nih.gov/pubmed/29712755) | B73 and sorghum | 11 tissues | |  | |
| Kaeppler2018 biomass | 78 inbreds | ? | 25M | |
| Briggs B73&Mo17 atlas | B73, Mo17 and BxM | | 20-30M | | [link](/briggs/data) |
| biomAP | 35 inbreds + 96 hybrids | 5 | 20-40M | | |


