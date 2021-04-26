# COVID-19_Islets

This repository contains the R scripts necessary to perform the analysis in paper [SARS-CoV-2 Infection Induces Beta Cell Transdifferentiation](https://github.com/shuibingchen/COVID-19_Islets).

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

