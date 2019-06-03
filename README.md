# CellSIUS: Cell Subtype Identification from Upregulated gene Sets
CellSIUS is an R package enabling the identification and characterization
of (rare) cell sub-populations from complex scRNA-seq datasets:
it takes as input expression values of N cells grouped into M(>1)
clusters. Within each cluster, genes with a bimodal distribution
are selected and only genes with cluster-specific expression are
retained. Among these candidate marker genes, sets with correlated
expression patterns are identified by graph-based clustering.
Finally, cells are assigned to subgroups based on their average
expression of each gene set. The CellSIUS algorithm output provides the
rare/ sub cell types by cell indices and their transcriptomic signatures.

[![DOI](https://zenodo.org/badge/189385961.svg)](https://zenodo.org/badge/latestdoi/189385961)

## Reference
Wegmann *et Al.*, **CellSIUS provides sensitive and specific detection of rare cell populations from complex single cell RNA-seq data**, Genome Biology, 2019 [accepted]

