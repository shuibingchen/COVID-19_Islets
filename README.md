# COVID-19_Islets

This repository contains the R scripts necessary to perform the analysis in our paper [Tang, Xuming et al. “SARS-CoV-2 Infection Induces Beta Cell Transdifferentiation.” *Cell Metabolism*, 19 May. 2021, doi:10.1016/j.cmet.2021.05.015](https://doi.org/10.1016/j.cmet.2021.05.015), as described in the supplementary methods and main text.

### Input data
The single cell RNA-seq data were generated with the 10X Chromium and pre-processed using 10X cellranger pipeline. The raw data are available in the GEO database with accession# [GSE159556](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159556).

### Requirements
The following R packages were used:
- Seurat
- scran
- scater
- slingshot
- batchelor
- gam
- future
- dplyr
- reshape2
- magrittr
- ggplot2
- pheatmap
- cowplot
- RColorBrewer

