# run_slingshot.R
# run trajectory analysis on beta cells on both donor (after integration using MNN)
# Author: Tuo Zhang
# Date: 07/19/2020
# Version: 1.0
# 

library(scran)
library(Seurat)
library(batchelor)
library(slingshot)
library(gam)
library(future)
library(dplyr)
library(magrittr)
library(RColorBrewer)
library(ggsci)

workdir <- "project.folder"
sourcedir <- file.path(workdir, "data")
figdir <- file.path(workdir, "figure", "slingshot")
infodir <- file.path(workdir, "info", "slingshot")

seurat.obj <- file.path(workdir, "info", "panc.rds")
rescaled.sce.obj <- file.path(workdir, "info", "rescaled.sce.list.rds")

# project name
project <- "islets"

# pattern for defining mitochondrial/ribosomal genes
mito.pattern <- "^MT-"
ribo.pattern <- "^RPL|^RPS"
cov2.pattern <- "^CoV2"

# load functions
setwd(workdir)
source("my_functions.R")

# set a random seed
set.seed(98)

# set parallelization in Seurat
plan("multiprocess", workers=6)
options(future.globals.maxSize=10*1024^3)

# read in dissociation-related genes
dissofile <- file.path(sourcedir, "dissociation", "dissociation_related_genes.human.txt")
disso.genes <- as.vector(read.table(dissofile, header=F, check.names=F, sep="\t")$V1)

# load seurat object
panc <- readRDS(seurat.obj)

# select beta cells in mock and infected samples from both donors
panc.combined <- subset(panc, subset=Condition %in% c('Mock','Infected'), ident=c(2))

# Identification of highly variable features (feature selection)
panc.combined %<>% FindVariableFeatures(selection.method="vst", nfeatures=3500)

# remove dissociation-related genes and ribosomal genes from variable gene list
# and select the remaining top 3000 genes for MNN-based correction
variable.genes <- setdiff(VariableFeatures(panc.combined), c(disso.genes, 
                                                             grep(mito.pattern, rownames(panc.combined), value=T), 
                                                             grep(ribo.pattern, rownames(panc.combined), value=T), 
                                                             grep(cov2.pattern, rownames(panc.combined), value=T)))
variable.genes <- head(variable.genes, 3000)

# load rescaled SingleCellExperiment object
rescaled.sce.list <- readRDS(rescaled.sce.obj)
# subsetting rescaled objects
target.cells <- rownames(subset(FetchData(panc, vars=c('ident','Condition')), Condition %in% c('Mock','Infected') & ident %in% c(2)))
target.rescaled.sce.list <- lapply(rescaled.sce.list[c('D1S1','D1S2','D2S1','D2S2')], function(x) { x[, intersect(colnames(x), target.cells)] })

# perform MMN-based correction
original <- lapply(target.rescaled.sce.list, function(x) {logcounts(x)[variable.genes,]})
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.merge=TRUE)))

# set column names
colnames(reducedDim(mnn.out)) = paste0("MNN_", 1:ncol(reducedDim(mnn.out)))

# add MNN correction results to Seurat object
panc.combined[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(mnn.out)[rownames(panc.combined@meta.data),], key="MNN_", assay=DefaultAssay(panc.combined))

# Run non-linear dimensional reduction (UMAP)
panc.combined %<>% RunUMAP(dims=1:50, reduction="mnn", n.components=3, seed.use=42, n.neighbors=30, n.epochs=500)

# Cluster the cells
# FindNeighbors: Shared Nearest Neighbor(SNN) Graph Construction
panc.combined %<>% FindNeighbors(reduction="mnn", dims=1:50)
# FindClusters
panc.combined %<>% FindClusters(resolution=seq(0.05,1,by=0.05), verbose=T)

# set cell identity
panc.combined %<>% SetIdent(value="RNA_snn_res.0.2")

# UMAP by Cluster
g <- myDimPlot(tobj=panc.combined, treduct="umap", tcate="ident", tsuffix="Cluster", tlabel=TRUE, tsplit=FALSE, tptsize=1) 
ggsave(file.path(figdir, "UMAPPlot.by_Cluster.labeled.png"), plot=g, width=11, height=8, dpi=300)

# save seurat object
saveRDS(panc.combined, paste(infodir,'panc.combined.rds',sep="/"))

# free memory
rm(panc)
rm(rescaled.sce.list)
rm(mnn.out)
rm(original)
rm(target.rescaled.sce.list)

# start to run slingshot
# remove cells from cluster 5 (likely alpha cells) + a couple of cells that are close to cluster 5
keep.cells <- rownames(subset(FetchData(panc.combined, vars=c('ident','UMAP_1','UMAP_2')), (! ident %in% c(5)) & UMAP_2 < 5))

# get raw counts
raw.counts <- as.matrix(panc.combined[["RNA"]]@counts[,keep.cells])

# get normalized counts
normalized.counts <- as.matrix(panc.combined[["RNA"]]@data[,keep.cells])

# prepare cell annotations
metadata <- FetchData(panc.combined, vars=c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo", "percent.virus", "Name", "Condition", "ident"), cells=keep.cells)
metadata$ident <- factor(metadata$ident, levels=0:4)

# create SingleCellExperiment object
sim <- SingleCellExperiment(assays=List(counts=raw.counts[,rownames(metadata)], norm=normalized.counts[,rownames(metadata)]), colData=metadata)

# add dimensionality reduction data
# UMAP
umap.rd <- FetchData(object=panc.combined, vars=c("UMAP_1","UMAP_2"), cells=keep.cells)
# MNN
mnn.rd <- FetchData(object=panc.combined, vars=c("MNN_1","MNN_2"), cells=keep.cells)

reducedDims(sim) <- SimpleList(UMAP=as.matrix(umap.rd[rownames(metadata),]), MNN=as.matrix(mnn.rd[rownames(metadata),]))

# run slingshot:
# - identifying global lineage structure
# - constructing smooth curves and ordering cells
sce <- slingshot(data=sim, clusterLabels="ident", reducedDim="UMAP", start.clus=0, end.clus=4)

# save slingshot object
saveRDS(sce, paste(infodir,"sce.slingshot.rds",sep='/'))

# Figure 5A
# Ordering beta cells along a transdifferentiation transition in mock and SARS-CoV-2 infected islets
# plot inferred lineage with cells colored by pseudo time
# curve coordinates
curve.coord <- as.data.frame(SlingshotDataSet(sce)@curves$curve1$s)
# reorder points (cells)
ordered.cells <- rownames(curve.coord)[SlingshotDataSet(sce)@curves$curve1$ord]
curve.coord <- curve.coord[ordered.cells,]
# point coordinates
point.coord <- merge(umap.rd, as.data.frame(colData(sce)[,c("slingPseudotime_1","Name","Condition","ident"),drop=F]), by=0)
rownames(point.coord) <- point.coord$Row.names
point.coord <- point.coord[ordered.cells,-1]
# plot
g <- ggplot(data=point.coord)
g <- g + geom_point(mapping=aes(x=UMAP_1, y=UMAP_2, color=slingPseudotime_1), shape=19, size=2, alpha=0.7)
g <- g + scale_colour_gradientn(colours=colorRampPalette(brewer.pal(11, "Spectral")[-6])(100), name="Pseudo-time")
g <- g + geom_path(data=curve.coord, mapping=aes(x=UMAP_1, y=UMAP_2), size=1, lineend="round")
g <- g + xlab("UMAP 1") + ylab("UMAP 2")
g <- g + theme_bw()
g <- g + theme(axis.text=element_blank(), axis.ticks=element_blank())
g <- g + theme(axis.title=element_text(size=18, face="bold"), legend.title=element_text(size=12, face="bold"))
ggsave(file.path(figdir, "trajectory.color_by_pseudo_time.v2.png"), width=8, height=6, dpi=300)

# identify temporally expressed genes
# find genes that change their expression over the course of development
# by selecting a set of variable genes and regressing each gene on the pseudo-time variable
# using a general additive model (GAM)

# pseudo-time
t <- sce$slingPseudotime_1

# normalized expression
texp <- assays(sce)$norm

# pre-filter genes (expressed in more than 1% of all cells)
ncells.per_gene <- apply(texp, MARGIN=1, FUN=function(x) { sum(x>0) })
filtered.genes <- rownames(texp)[ncells.per_gene > ncol(texp) * 0.01]
length(filtered.genes)
# 13001 genes

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(texp[filtered.genes,], 1, function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

# multiple testing by Benjamini & Hochberg (1995)
gam.pval.out <- data.frame(pval=gam.pval, padj=p.adjust(gam.pval, method="BH"))

# output gam results
write.table(gam.pval.out[with(gam.pval.out, order(padj)),], file=file.path(infodir, "gam.pval.txt"), quote=FALSE, na="NA", sep="\t", col.names=NA)

# manually generate heatmap plot
myPlotHeatmap <- function(texp, tgenes, tpseudotime, tpt.colors=colorRampPalette(brewer.pal(11, "Spectral")[-6])(100), tcellann=NULL, tanncolors=list(), toutfile, tcolor, tmax=2.5, tmin=-2.5, twidth=15, theight=12, tunits="in", tres=300, tfsize=7, tfsrow=6){
  # reorder cells by pseudo-time
  torder.cells <- names(tpseudotime)[order(tpseudotime, na.last=NA)]
  # get valid gene list
  tvalid.genes <- intersect(rownames(texp), tgenes)
  # prepare expression matrix
  texp.raw <- texp[tvalid.genes, torder.cells]
  # z-score transform expression matrix
  texp.zscore <- t(scale(t(texp.raw)))
  # set minimum/maximum display value (clip all values that are out of range)
  texp.zscore <- apply(texp.zscore, 2, function(x) ifelse(x > tmax, tmax, x))
  texp.zscore <- apply(texp.zscore, 2, function(x) ifelse(x < tmin, tmin, x))
  # add time to column annotations
  tanncols <- data.frame(PseudoTime=tpseudotime[torder.cells])
  rownames(tanncols) <- torder.cells
  if (! is.null(tcellann)){
    tanncols <- tcellann[torder.cells,,drop=F]
    tanncols$PseudoTime <- tpseudotime[torder.cells]
  }
  tanncolors[['PseudoTime']] <- tpt.colors
  # plot heatmap
  png(file=toutfile, width=twidth, height=theight, units=tunits, res=tres, pointsize=6)
  tout <- pheatmap(texp.zscore, color=tcolor, scale="none", cluster_rows=T, cluster_cols=F, annotation_col=tanncols, annotation_color=tanncolors, show_colnames=F, show_rownames=T, fontsize=tfsize, fontsize_row=tfsrow)
  dev.off()
}

# manually generate scatter plot (expression v.s. pseudo-time) for a given gene
myPlotExpTime <- function(texp, tgene, tpseudotime, toutfile, tcolor, twidth=10, theight=6, tunits="in", tres=300){
  # valid gene?
  if (! tgene %in% rownames(texp)){
    print(paste("gene",tgene,"is","unavailable.",sep=" "))
    return(NULL)
  }
  # reorder cells by pseudo-time
  torder.cells <- names(tpseudotime)[order(tpseudotime, na.last=NA)]
  # prepare data frame for plotting
  texp.forPlot <- data.frame(expression=texp[tgene, torder.cells], pseudotime=tpseudotime[torder.cells])
  # local polynomial regression fitting (loess)
  loessMod75 <- loess(expression ~ pseudotime, data=texp.forPlot, span=0.75)
  loessMod100 <- loess(expression ~ pseudotime, data=texp.forPlot, span=1)
  # get smoothed output
  texp.forPlot$smoothed75 <- predict(loessMod75) 
  texp.forPlot$smoothed100 <- predict(loessMod100) 
  # plot
  tg <- ggplot(texp.forPlot)
  tg <- tg + geom_point(aes(x=pseudotime, y=expression, color=pseudotime), shape=19, size=2, alpha=0.7)
  tg <- tg + scale_colour_gradientn(colours=tcolor)
  tg <- tg + geom_path(mapping=aes(x=pseudotime, y=smoothed100), color='gray20', size=1, lineend="round")
  tg <- tg + ggtitle(tgene)
  tg <- tg + theme_bw()
  tg <- tg + theme(plot.title=element_text(hjust=0.5, face="bold", size=16))
  ggsave(toutfile, width=twidth, height=theight, dpi=tres)
}

# color
my.color <- colorRampPalette(c("#2166ac","#f7f7f7","#b2182b"))(100)

# pseudo-time
t <- sce$slingPseudotime_1
names(t) <- colnames(assays(sce)$norm)

# Figure 5B
# Changes in expression of INS, GCG, CPA1, PRSS1, and PRSS2 during beta cell transdifferentiation.
for (gene in c('INS','GCG','CPA1','PRSS1','PRSS2')){
  print(paste("processing",gene,sep=" "))
  myPlotExpTime(texp=assays(sce)$norm, tgene=gene, tpseudotime=t,
                toutfile=file.path(figdir,paste("expression","vs","pseudotime",gene,"png",sep=".")),
                tcolor=colorRampPalette(brewer.pal(11, "Spectral")[-6])(100))
}

# Figure 5C
# Heatmap showing expression changes in INS, GCG, CPA1, CPA2, CPB1, PRSS1, PRSS2, SARS-CoV-2-N, and SARS-CoV-2-ORF1ab during beta cell transdifferentiation.
genes <- c('CPA2','CoV2-N','CPB1','CPA1','GCG','PRSS1','PRSS2','INS','CoV2-orf1ab')
myPlotHeatmap(texp=assays(sce)$norm, tgenes=genes, tpseudotime=t, 
              tcellann=NULL, tanncolors=list(),
              toutfile=file.path(figdir,"heatmap.selected_genes.png"),
              tcolor=my.color, twidth=10, theight=1.9, tres=300, tfsize=8, tfsrow=8)

# Figure 5E
# Heatmap showing expression changes in eIF2 pathway-associated genes during beta cell transdifferentiation.
# manually selected EIF pathway genes
eif.genes <- c('RPS29','INS','RPS21','HNRNPA1','EIF3K','RPL27','RPL27A','RPL31','RPL10A','RPLP2','RPL13A','PRS16','RPL36','RPS5',
               'FAU','RPL18','RPL6','RPL28','RPL12','RPS15','RPS7','RPL3','RPS18','RPL32','RPL11','RPS8','RPL7A','RPS3','RPL41',
               'RPL13','RPL8','RPLP1','RPS12','RPLP0','RPS4X','RPS24','RPS25','RPS26','RPS13','RPL10','RPS14','RPL30','RPL5','RPL14',
               'RPL9','RPL24','RPS23','RPS6','RPS3A','RPL15','RPL19','PRS27A','PRL18A','PRL29','RPL23A','PRSA','RPS19','RPL35',
               'RPL37A','RPL39','RPL37','RPS27','RPL21','RPS15A','RPS28','RPS9','RPL34','RPL35A','RPL7','RPL4','EIF1','RPL22',
               'RPS11','UBA52','RPL38','RPS20','EIF3E','EIF3L','EIF4A2','EIF3F','RPL17','RPL26','EIF3H','RPL36AL')
myPlotHeatmap(texp=assays(sce)$norm, tgenes=c(eif.genes, 'INS', 'CoV2-N'), tpseudotime=t, 
              tcellann=NULL, tanncolors=list(),
              toutfile=file.path(figdir,"heatmap.EIF.pathway.genes.png"),
              tcolor=my.color, twidth=10, theight=16, tres=300, tfsize=8, tfsrow=8)
# ----------------------------------------------------------------------------------------------------------------- #

sessionInfo()
