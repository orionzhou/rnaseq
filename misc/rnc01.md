# B73 developmental atlas

- [9 previous studies](data/11_qc/rnc01/00.studies.tsv):
  - a total of 247 distinct tissue / developmental stages
- [hierarchical clustering](data/11_qc/rnc01/21.cpm.hclust.pdf)
- [tSNE plot](data/11_qc/rnc01/25.tsne.pdf)
- [Classification of 39,591 protein-coding genes based on expression broadness](data/11_qc/rnc01/31.tis.expression.pdf):
  - `Silent`: never expressed
  - `Tissue specific`: expressed in less than 20% tissues (out of total 247)
  - `Intermediate`: expressed in 20-80% tissues
  - `Constitutive`: expressed in >80% tissues
- [Classification table of 39,391 protein coding genes based on expression broadness](data/11_qc/rnc01/30.tis.expression.tsv.gz) with columns:
  - `gid`: B73 v4 gene ID
  - `n.tis`: number of tissues the gene showed detectable expression (CPM >= 1)
  - `etag`: classification of expression broadness
  - `tiss`: tissues in which the gene is expressed