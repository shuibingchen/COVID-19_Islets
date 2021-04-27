# subcluster_beta.R
# perform subclustering on beta cells
# Author: Tuo Zhang
# Date: 07/19/2020
# Version: 1.0
# 

library(scran)
library(Seurat)
library(batchelor)
library(scater)
library(future)
library(dplyr)
library(magrittr)
library(ggrepel)

workdir <- "project.folder"
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure", "beta")
infodir <- file.path(workdir, "info", "beta")

# project name
project <- "islets"

# pattern for defining mitochondrial/ribosomal/SARS-CoV-2 genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"
cov2.pattern <- "^CoV2"

# load functions
setwd(workdir)
source("my_functions.R")

# set a random seed
set.seed(98)

# set parallelization for Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=10*1024^3)

# read in dissociation-related genes
dissofile <- file.path(sourcedir, "dissociation", "dissociation_related_genes.human.txt")
disso.genes <- as.vector(read.table(dissofile, header=F, check.names=F, sep="\t")$V1)

# -------------------------------------------------- extract beta cells ------------------------------------------------- #
# sample info
sample.info <- data.frame(Name=c("Healthy_7_Mock","Healthy_7_SARS-CoV-2","Healthy_8_Mock","Healthy_8_SARS-CoV-2","Healthy_8_SARS-CoV-2_ISIRB"), 
                          Condition=c("Mock","Infected","Mock","Infected","Treated"), 
                          Donor=c("Donor7","Donor7","Donor8","Donor8","Donor8"))
rownames(sample.info) <- c("D1S1","D1S2","D2S1","D2S2","D2S3")

# load Seurat object
panc <- readRDS(file.path(workdir, "info", "panc.rds"))

# subset to get beta cells
panc.beta <- subset(panc, idents=c(2))

# free memory
rm(panc)
# ----------------------------------------------------------------------------------------------------------------------- #

# ---------------------------------------------- run MNN-based correction ----------------------------------------------- #
# load rescaled SingleCellExperiment object
rescaled.sce.list <- readRDS(file=file.path(workdir, "info", "rescaled.sce.list.rds"))

# subset beta cells only
rescaled.sce.list %<>% lapply(function(x) {x[, intersect(colnames(x), colnames(panc.beta))]})

# Identification of highly variable features (feature selection)
panc.beta %<>% FindVariableFeatures(selection.method="vst", nfeatures=3500)

# remove dissociation-related genes and mitochondrial/ribosomal/SARS-CoV-2 genes from variable gene list
# and select the remaining top 3000 genes for MNN-based correction
variable.genes <- setdiff(VariableFeatures(panc.beta), c(disso.genes, 
                                                         grep(mito.pattern, rownames(panc.beta), value=T), 
                                                         grep(ribo.pattern, rownames(panc.beta), value=T), 
                                                         grep(cov2.pattern, rownames(panc.beta), value=T)))
variable.genes <- head(variable.genes, 3000)

# perform MMN-based correction
original <- lapply(rescaled.sce.list, function(x) {logcounts(x)[variable.genes,]})
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.merge=TRUE)))

# set column names
colnames(reducedDim(mnn.out)) = paste0("MNN_", 1:ncol(reducedDim(mnn.out)))

# add MNN correction results to Seurat object
panc.beta[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(mnn.out)[rownames(panc.beta@meta.data),], key="MNN_", assay=DefaultAssay(panc.beta))

# free memory
rm(mnn.out)
rm(original)
rm(rescaled.sce.list)
# ----------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------- perform dimension reduction and clustering -------------------------------------- #
# Run non-linear dimensional reduction (UMAP)
panc.beta %<>% RunUMAP(dims=1:50, reduction="mnn", n.components=3, seed.use=42, n.neighbors=35, n.epochs=1000)

# Cluster the cells
# FindNeighbors: Shared Nearest Neighbor(SNN) Graph Construction
panc.beta %<>% FindNeighbors(reduction="mnn", dims=1:50)
# FindClusters
panc.beta %<>% FindClusters(resolution=seq(0.05,1,by=0.05), verbose=T)

# set cell identity
panc.beta %<>% SetIdent(value="RNA_snn_res.0.1")

# reorder clusters
Idents(panc.beta) <- factor(Idents(panc.beta), levels=0:(length(unique(Idents(panc.beta)))-1))

# Finding differentially expressed features (cluster biomarkers)
for (k in c(0:(length(unique(Idents(panc.beta)))-1))){
  # wilcox test
  myDETest(panc.beta, k,"wilcox",FALSE,infodir,figdir,tassay="RNA")
}
# ----------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------------- merge clusters ---------------------------------------------------- #
# save current cluster
panc.beta[["base.clust"]] <- Idents(panc.beta)

# merge the following clusters
# C0 + C1 + C2 + C3           ==>  C0 (GCG-/INS+)
# C4                          ==>  C2 (GCG+/INS+)
merge.clust <- c(0,0,0,0,1)
names(merge.clust) <- 0:4

final.clust <- as.vector(merge.clust[as.vector(Idents(panc.beta))])
names(final.clust) <- names(Idents(panc.beta))
final.clust <- factor(final.clust, levels=0:1)

# add final clusters to meta data
panc.beta[["final.clust"]] <- final.clust

# set final clusters
Idents(panc.beta) <- final.clust

# save seurat object
saveRDS(panc.beta, file=file.path(infodir, "panc.beta.rds"))

# free memory
rm(merge.clust)
rm(final.clust)
# ----------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------------- generate figures ------------------------------------------------ #
# Figure S4A
# UMAP plot showing the two sub-clusters of beta cell population
# load all cell Seurat object
panc <- readRDS(file.path(workdir,"info","panc.rds"))
# prepare cell labeling
subcluster.labels <- FetchData(panc.beta, vars='ident')
non.beta.cells <- setdiff(colnames(panc), colnames(panc.beta))
non.beta.labels <- data.frame(ident=rep("None",length(non.beta.cells)))
rownames(non.beta.labels) <- non.beta.cells
all.labels <- rbind(subcluster.labels, non.beta.labels)
colnames(all.labels) <- c('subcluster')
# assign to seurat object
panc %<>% AddMetaData(metadata = all.labels[rownames(panc@meta.data),,drop=F])
# set color
my.color <- c("None"="#d9d9d9","0"="#a6cee3","1"="#fb9a99")
# plot UMAP
g <- myDimPlot(tobj=panc, treduct="umap", tcate="subcluster", tsuffix="subcluster", tcolor=my.color, tlabel=TRUE, tsplit=FALSE, tptsize=1)
ggsave(file.path(figdir, "UMAPPlot.highlight_refined_subcluster.png"), plot=g, width=11, height=8, dpi=300)

# Figure S4B
# Dot plot illustrating expression level of INS, GCG, PRSS1, PRSS2 in two sub-clusters of mock versus SARS-CoV-2 infected human islets at 24 hpi (MOI=1).
for (c in as.vector(unique(Idents(panc.beta)))){
  tcells <- rownames(subset(FetchData(panc.beta, vars=c('ident', 'Condition')), ident %in% c & Condition %in% c('Mock','Infected')))
  print(paste('cluster',c,'Mock vs Infected:',length(tcells),'cells.'))
  plot <- myDotPlot.2(tobj=panc.beta, tgenes=c('INS','PRSS1','PRSS2','GCG'), tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=tcells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
  ggsave(file.path(figdir, paste('Dot.exp.INS_PRSS1_2_GCG.SARS-CoV-2_infected.vs.Mock.beta_cells','refined','subcluster',c,'png', sep='.')), plot=plot, width=4.5, height=4, dpi=300)
}
# ----------------------------------------------------------------------------------------------------------------------- #

sessionInfo()
