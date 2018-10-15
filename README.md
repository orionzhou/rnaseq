# Collection of large-scale maize expression datasets

This repo provides a collection of recently published (after 2010) maize gene expression datasets. The focus is on large-scale (multiple tissues / developmental stages, multiple inbred / hybrid lines) Illumina RNA-Seq experiments, but also inlcude experiments done with other sequencing technologies (3' RNA-Seq, MNase-Seq, etc.).

Raw sequencing reads were downloaded from [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra), trimmed using [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) and mapped to the [maize B73 AGP_v4 genome](http://plants.ensembl.org/Zea_mays/Info/Index) using [STAR](https://github.com/alexdobin/STAR).  Uniquely mapped reads were assigned to and counted for the 46,117 reference gene models ([Ensembl Plants v37](http://plants.ensembl.org/Zea_mays/Info/Index)) using [featureCounts](http://bioinf.wehi.edu.au/featureCounts/).  Raw read counts were then normalized using the [TMM normalization approach](https://bioconductor.org/packages/release/bioc/html/edgeR.html) to give CPMs (Counts Per Million reads) and then further normalized by gene CDS lengths to give FPKM (Fragments Per Kilobase of exon per Million reads) values.  Hierarchical clustering and principal component analysis were used to explore sample clustering pattern.

See [this table](/data/01.cfg.tsv) for a list of collected datasets.

