# run_cluster.mnn.R
# run MNN-based alignment with Scran, followed by clustering with Seurat3
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
figdir <- file.path(workdir, "figure")
infodir <- file.path(workdir, "info")

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

# ----------------------------------------------------- filter cells ---------------------------------------------------- #
# sample info
sample.info <- data.frame(Name=c("Healthy_7_Mock","Healthy_7_SARS-CoV-2","Healthy_8_Mock","Healthy_8_SARS-CoV-2","Healthy_8_SARS-CoV-2_ISIRB"), 
                          Condition=c("Mock","Infected","Mock","Infected","Treated"), 
                          Donor=c("Donor7","Donor7","Donor8","Donor8","Donor8"))
rownames(sample.info) <- c("D1S1","D1S2","D2S1","D2S2","D2S3")

# load raw UMI counts table per patient
raw.counts.list <- list()
for (k in 1:nrow(sample.info)){
  pid <- rownames(sample.info)[k]
  sid <- sample.info$Name[k]
  raw.counts.list[[k]] <- my.Read10X(file.path(sourcedir, sid, "filtered_feature_bc_matrix"), pid)
}
names(raw.counts.list) <- rownames(sample.info)

# merge raw UMI counts tables
raw.counts.all <- my.MergeMatrix(raw.counts.list)

# Initialize the Seurat object with the raw (non-normalized data).
panc.initial <- CreateSeuratObject(counts=raw.counts.all, project=project, assay="RNA", min.cells=0, min.features=0, 
                                   names.field=1, names.delim="_", meta.data=NULL)

# Calculates the content of mitochondrial/ribosomal/SARS-CoV-2 genes per cell
panc.initial[["percent.mito"]] <- PercentageFeatureSet(panc.initial, pattern=mito.pattern)
panc.initial[["percent.ribo"]] <- PercentageFeatureSet(panc.initial, pattern=ribo.pattern)
panc.initial[["percent.virus"]] <- PercentageFeatureSet(panc.initial, pattern=cov2.pattern)

# Add sample information to the Seurat object
tmeta <- data.frame(row.names=rownames(panc.initial@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc.initial@meta.data[,"orig.ident"])])
}
panc.initial %<>% AddMetaData(metadata=tmeta)

# perform cell filtering
# nGene > 500, nGene <= 6000, nUMI > 1000, nUMI <= 60000, percent.mito < 15%
panc.initial %<>% subset(subset=nFeature_RNA > 500 & nFeature_RNA <= 6000 & nCount_RNA > 1000 & nCount_RNA <= 60000 & percent.mito < 15)

# free memory
rm(raw.counts.list)
rm(raw.counts.all)
rm(tmeta)
rm(tdic)
# ----------------------------------------------------------------------------------------------------------------------- #

# ---------------------------------------------- run MNN-based correction ----------------------------------------------- #
# prepare raw UMI counts table from each donor
selected.donors <- c("D1S1", "D1S2", "D2S1", "D2S2", "D2S3")
sample.list <- list()
for (donor in selected.donors){
  sample.list[[donor]] <- panc.initial[["RNA"]]@counts[, rownames(subset(panc.initial@meta.data, orig.ident == donor))]
}

# create SingleCellExperiment object
sce.list <- list()
for (donor in names(sample.list)){
  sce.list[[donor]] <- SingleCellExperiment(list(counts=as.matrix(sample.list[[donor]])))
}

# run a pre-clustering to avoid pooling together very different cells
# normalization will be performed for cells within each cluster
preclust.list <- lapply(sce.list, function(x) quickCluster(x=x, min.size=200, assay.type="counts", method="hclust", min.mean=0.1))

# normalize data by deconvolving size factors from cell pools
sce.list <- mapply(FUN=function(x,y) {computeSumFactors(x=x, min.mean=0.1, cluster=y)}, x=sce.list, y=preclust.list)

# compute normalized log-expression values
sce.list %<>% lapply(FUN=function(x) {normalize(object=x)})

# rescale among donors
rescaled.sce.list <- do.call(multiBatchNorm, sce.list)

# create a seurat object with raw UMI counts
panc <- CreateSeuratObject(counts=as(do.call(cbind, lapply(rescaled.sce.list, function(x) counts(x))), "dgCMatrix"), 
                           project=project, assay="RNA", min.cells=0, min.features=0,
                           names.field=1, names.delim="_", meta.data=NULL)

# Calculates the content of mitochondrial/ribosomal/SARS-CoV-2 genes per cell
panc[["percent.mito"]] <- PercentageFeatureSet(panc, pattern=mito.pattern)
panc[["percent.ribo"]] <- PercentageFeatureSet(panc, pattern=ribo.pattern)
panc[["percent.virus"]] <- PercentageFeatureSet(panc, pattern=cov2.pattern)

# Add sample information
tmeta <- data.frame(row.names=rownames(panc@meta.data))
for (tx in colnames(sample.info)){
  tdic <- as.vector(sample.info[,tx])
  names(tdic) <- rownames(sample.info)
  tmeta[,tx] <- as.vector(tdic[as.vector(panc@meta.data[,"orig.ident"])])
}
panc %<>% AddMetaData(metadata=tmeta)


# replace normalized data with the scran normalized data
panc[["RNA"]]@data <- as(do.call(cbind, lapply(rescaled.sce.list, function(x) logcounts(x))) * log(2), "dgCMatrix")

# Identification of highly variable features (feature selection)
panc %<>% FindVariableFeatures(selection.method="vst", nfeatures=3500)

# remove dissociation-related genes and ribosomal genes from variable gene list
# and select the remaining top 3000 genes for MNN-based correction
variable.genes <- setdiff(VariableFeatures(panc), c(disso.genes, 
                                                    grep(mito.pattern, rownames(panc), value=T), 
                                                    grep(ribo.pattern, rownames(panc), value=T), 
                                                    grep(cov2.pattern, rownames(panc), value=T)))
variable.genes <- head(variable.genes, 3000)

# perform MMN-based correction
original <- lapply(rescaled.sce.list, function(x) {logcounts(x)[variable.genes,]})
mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, auto.merge=TRUE)))
# set column names
colnames(reducedDim(mnn.out)) = paste0("MNN_", 1:ncol(reducedDim(mnn.out)))

# add MNN correction results to Seurat object
panc[["mnn"]] <- CreateDimReducObject(embeddings=reducedDim(mnn.out)[rownames(panc@meta.data),], key="MNN_", assay=DefaultAssay(panc))

# free memory
rm(preclust.list)
rm(panc.initial)
rm(sce.list)
rm(tmeta)
rm(tdic)
rm(mnn.out)
rm(original)
rm(rescaled.sce.list)
# ----------------------------------------------------------------------------------------------------------------------- #

# ------------------------------------- perform dimension reduction and clustering -------------------------------------- #
# Run non-linear dimensional reduction (UMAP)
panc %<>% RunUMAP(dims=1:50, reduction="mnn", n.components=3, seed.use=42, n.neighbors=30, n.epochs=4000)

# Cluster the cells
# FindNeighbors: Shared Nearest Neighbor(SNN) Graph Construction
panc %<>% FindNeighbors(reduction="mnn", dims=1:50)
# FindClusters
panc %<>% FindClusters(resolution=seq(0.05,1,by=0.05), verbose=T)

# set cell identity
panc %<>% SetIdent(value="RNA_snn_res.0.4")

# reorder clusters
Idents(panc) <- factor(Idents(panc), levels=0:(length(unique(Idents(panc)))-1))

# Finding differentially expressed features (biomarkers for each cluster)
for (k in c(0:(length(unique(Idents(panc)))-1))){
  print(paste("C",k,sep=""))
  # wilcox test
  myDETest(panc, k,"wilcox",FALSE,infodir,figdir,tassay="RNA",tsuf=NULL)
}
# ----------------------------------------------------------------------------------------------------------------------- #

# --------------------------------------------------- merge clusters ---------------------------------------------------- #
# save current cluster
panc[["base.clust"]] <- Idents(panc)

# merge clusters after manual review
# C0 + C5 + C8                ==>  C0 (Acinar cells)
# C1 + C2 + C10               ==>  C1 (Alpha cells)
# C4 + C7 + C9 + C19          ==>  C2 (Beta cells)
# C3 + C12 + C13 + C16 + C18  ==>  C3 (Ductal cells)
# C6                          ==>  C4 (Fibroblast)
# C11                         ==>  C5 (Delta cells)
# C14                         ==>  C6 (PP cells)
# C15                         ==>  C7 (Endothelial cells)
# C17                         ==>  C8 (Immune cells)
merge.clust <- c(0,1,1,3,2,0,4,2,0,2,1,5,3,3,6,7,3,8,3,2)
names(merge.clust) <- 0:19

final.clust <- as.vector(merge.clust[as.vector(Idents(panc))])
names(final.clust) <- names(Idents(panc))
final.clust <- factor(final.clust, levels=0:8)

# add final clusters to meta data
panc[["final.clust"]] <- final.clust

# set final clusters
Idents(panc) <- final.clust

# Finding differentially expressed features (biomarkers for each cluster)
for (k in c(0:8)){
  print(paste("Cluster",k))
  # wilcox test
  myDETest(panc, k,"wilcox",FALSE,infodir,figdir,tassay="RNA",tsuf="merged_clusters")
}

# save seurat object
saveRDS(panc, file=file.path(infodir, "panc.rds"))

# free memory
rm(merge.clust)
rm(final.clust)
# ----------------------------------------------------------------------------------------------------------------------- #

# ----------------------------------------------------- generate figures ------------------------------------------------ #
# get x/y-axis boundaries in UMAP
apply(FetchData(panc, vars=c("UMAP_1","UMAP_2")), MARGIN=2, FUN=min)
#    UMAP_1    UMAP_2
# -7.756379 -7.701544 
apply(FetchData(panc, vars=c("UMAP_1","UMAP_2")), MARGIN=2, FUN=max)
#    UMAP_1   UMAP_2
#  7.072962 5.866345
# set x/y-axis boundaries 
x.upper <- 7.5
x.lower <- -8
y.upper <- 6
y.lower <- -8

# set colors for cell types
my.cluster.color.1 <- c('0'='#fb9a99','1'='#a6cee3','2'='#33a02c','3'='#b2df8a','4'='#1f78b4','5'='#fdbf6f','6'='#ff7f00','7'='#cab2d6','8'='#e31a1c')

# set color for MyFeaturePlot
myExpLowColor <- '#d9d9d9'
myExpHighColor <- '#b30000'

# Figure 2A
# UMAP plots presenting the nine cell types
g <- myDimPlot(tobj=panc, treduct="umap", tcate="ident", tsuffix="Cluster", tcolor=my.cluster.color.1, tlabel=FALSE, tsplit=FALSE, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tptsize=1.5, talpha=0.7, tltsize=20, tatlsize=22)
ggsave(file.path(figdir, "UMAPPlot.by_Cluster.png"), plot=g, width=10, height=8, dpi=300)

# Figure S2C
# UMAP plots presenting the nine cell types in each sample 
plot <- myDimPlot3(tobj=panc, treduct='umap', tgroup_by='ident', tgroup_order=0:8, thighlight=NULL, tsuffix='cluster', tcells=NULL, tcolor=my.cluster.color.1, 
                   tlabel=TRUE, tsplit_by='Name', tsplit_order=c('ICRH134_Mock','ICRH134_SARS-CoV-2','HPAP066_Mock','HPAP066_SARS-CoV-2','HPAP066_SARS-CoV-2_ISIRB'), 
                   txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=3, tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20)
ggsave(paste(revision.figdir, "UMAPPlot.by_cluster.split_by_sample.png", sep="/"), plot=plot, width=16, height=9, dpi=300)

# Find markers for every cluster compared to all remaining cells
# report only the positive ones (for plotting)
markers.wilcox.all <- FindAllMarkers(panc, only.pos=TRUE, logfc.threshold=0.25, test.use="wilcox", min.pct=0.25, assay="RNA")
# filter dissociation related genes, mitochondrial/ribosomal/SARS-CoV-2 genes
markers.wilcox.cleaned <- subset(markers.wilcox.all, ! gene %in% c(disso.genes, 
                                                                   grep(mito.pattern, rownames(panc), value=T), 
                                                                   grep(ribo.pattern, rownames(panc), value=T), 
                                                                   grep(cov2.pattern, rownames(panc))))
# select top ranked marker genes per cluster
markers.wilcox.top10 <- markers.wilcox.cleaned %>% group_by(cluster) %>% top_n(n=10, wt=avg_logFC)

# scaling the data on selected top ranked marker genes and SARS-CoV-2 genes
panc %<>% ScaleData(features=c(unique(markers.wilcox.top10$gene), grep(cov2.pattern, rownames(panc), value=T)))

# Figure S2A
# heatmap plot highlighting top ranked marker genes defining the nine cell types in human islets
myHeatmap(tobj=panc, tgenes=unique(markers.wilcox.top10$gene), toutfile=file.path(figdir,'myheatmap.pos.markers.wilcox.top10.png'), tassay="RNA", tgroup_by='ident', tgroup_order=0:8, tmax=2.5, tmin=-2.5, tcolor=colorRampPalette(c("#4575B4", "#f7f7f7", "#D73027"))(100), tanncolor=my.cluster.color.1, twidth=8400, theight=9600, tunits="px", tfontsize_row=10, tres=600)

# Visualize expressions of selected genes in UMAP/Violin plots
# Figure 2B
# expression levels of SARS-CoV-2 entry factors, including FURIN and CTSL
fig.2b.genes <- c('FURIN','CTSL')
for (gene in fig.2b.genes){
  print(gene)
  # UMAP
  plot <- MyFeaturePlot(tobj=panc, tgenes=gene, tcells=NULL, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=1, tlegend=NULL)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"png",sep='.')), plot=plot, width=10, height=8, dpi=300)
  # Violin plot
  plot <- MyExpViolin(tobj=panc, tgene=gene, tgroup_by='ident', tgroup_order=0:8, tcolor_by='ident', tcolor_order=0:8, tcolor=my.cluster.color.1, tcells=NULL, tassay="RNA", tncol=1)
  ggsave(file.path(figdir, paste("Violin","exp",gene,"png",sep=".")), plot=plot, height=4, width=6, dpi=300)
}

# Figure S2B
# expression of pancreatic cell markers
fig.s2b.genes <- c('KRT19','PRSS1','INS','GCG','COL1A1','LAPTM5','SST','PECAM1','PPY')
for (gene in fig.s2b.genes){
  print(gene)
  # UMAP
  plot <- MyFeaturePlot(tobj=panc, tgenes=gene, tcells=NULL, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=1, tlegend=NULL)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"png",sep='.')), plot=plot, width=10, height=8, dpi=300)
  # Violin plot
  plot <- MyExpViolin(tobj=panc, tgene=gene, tgroup_by='ident', tgroup_order=0:8, tcolor_by='ident', tcolor_order=0:8, tcolor=my.cluster.color.1, tcells=NULL, tassay="RNA", tncol=1)
  ggsave(file.path(figdir, paste("Violin","exp",gene,"png",sep=".")), plot=plot, height=4, width=6, dpi=300)
}

# Figure 2B and Figure S3A
# expression levels of SARS-CoV-2 genes in the mock and infected human islets.
cells.mock <- rownames(subset(FetchData(panc, vars=c('Condition')), Condition %in% c('Mock')))
cells.infected <- rownames(subset(FetchData(panc, vars=c('Condition')), Condition %in% c('Infected')))
cov2.genes <- c('CoV2-E','CoV2-M','CoV2-orf1ab','CoV2-ORF8','CoV2-ORF10','CoV2-S')
for (gene in cov2.genes){
  print(gene)
  # UMAP - mock
  plot <- MyFeaturePlot(tobj=panc, tgenes=c(gene), tcells=cells.mock, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=1, tlegend=NULL)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"Mock","png",sep='.')), plot=plot, width=10, height=8, dpi=300)
  # UMAP - infected
  plot <- MyFeaturePlot(tobj=panc, tgenes=c(gene), tcells=cells.infected, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=1, tlegend=NULL)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"Infected","png",sep='.')), plot=plot, width=10, height=8, dpi=300)
  # Violin plot
  plot <- MyExpViolin(tobj=panc, tgene=gene, tgroup_by='ident', tgroup_order=0:8, tsplit_by='Condition', tsplit_order=c('Mock','Infected'), tcolor_by='ident', tcolor_order=0:8, tcolor=my.cluster.color.1, tcells=c(cells.mock, cells.infected), tassay="RNA", tncol=1)
  ggsave(file.path(figdir, paste("Violin","exp",gene,"Infected.vs.Mock","png",sep=".")), plot=plot, height=8, width=6, dpi=300)
}

# Figure S4C
# UMAP plot showing the expression level of PRSS2 of SARS-CoV-2 infected human islets.
plot <- MyFeaturePlot(tobj=panc, tgenes=c('PRSS2'), tcells=cells.infected, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=1, tlegend=NULL)
ggsave(file.path(figdir, paste("UMAP","exp","PRSS2","Infected","png",sep='.')), plot=plot, width=10, height=8, dpi=300)

# Figure S6C
# UMAP and violin plots showing the expression levels of SARS-CoV-2 genes of control or 10 μM trans-ISRIB treated human islets at 24 hpi
# HPAP066_SARS-CoV-2_ISIRB.vs.HPAP066_SARS-CoV-2
cells.Healthy_8.infected <- rownames(subset(FetchData(panc, vars=c('Name')), Name == 'Healthy_8_SARS-CoV-2'))
cells.Healthy_8.treated <- rownames(subset(FetchData(panc, vars=c('Name')), Name == 'Healthy_8_SARS-CoV-2_ISIRB'))
for (gene in cov2.genes){
  print(gene)
  # UMAP - infected
  plot <- MyFeaturePlot(tobj=panc, tgenes=c(gene), tcells=cells.Healthy_8.infected, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=1, tlegend=NULL)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"Healthy_8_SARS-CoV-2","png",sep='.')), plot=plot, width=10, height=8, dpi=300)
  # UMAP - treated
  plot <- MyFeaturePlot(tobj=panc, tgenes=c(gene), tcells=cells.Healthy_8.treated, tassay="RNA", treduction.name="umap", tlowcolor=myExpLowColor, thighcolor=myExpHighColor, txlim=c(x.lower,x.upper), tylim=c(y.lower,y.upper), tncol=1, tlegend=NULL)
  ggsave(file.path(figdir, paste("UMAP","exp",gene,"Healthy_8_SARS-CoV-2_ISIRB","png",sep='.')), plot=plot, width=10, height=8, dpi=300)
  # Violin plot
  plot <- MyExpViolin(tobj=panc, tgene=gene, tgroup_by='ident', tgroup_order=0:8, tsplit_by='Name', tsplit_order=c('Healthy_8_SARS-CoV-2','Healthy_8_SARS-CoV-2_ISIRB'), tcolor_by='ident', tcolor_order=0:8, tcolor=my.cluster.color.1, tcells=c(cells.Healthy_8.infected,cells.Healthy_8.treated), tassay="RNA", tncol=1)
  ggsave(file.path(figdir, paste("Violin","exp",gene,"Healthy_8_SARS-CoV-2_ISIRB.vs.Healthy_8_SARS-CoV-2","png",sep=".")), plot=plot, height=8, width=6, dpi=300)
}

# Figure 3A
# Scoring the chemokine and cytokine expression levels in mock versus SARS-CoV-2 infected human islets at 24 hpi. 
chemokine.genes <- c('CCL3','MMP9','CCL11','IL1RN','CXCL1','IL4','CSF2','CCL28','CCL8','CXCL2','CCL13','CXCL16','CCL5','CCL4','CCL20','TNF','CXCL5','CXCL10','IL1B','IL6','CCL2','CCL7')
# collect cells from infected v.s. mock (merge two donors) only
mock.infected.cells <- rownames(subset(FetchData(panc, vars=c('Condition')), Condition %in% c('Mock','Infected')))
# plot
plot <- myScatterCdt(tobj=panc, tgenes=chemokine.genes, tgroup_by='Condition', tgroup_x='Mock', tgroup_y='Infected', tassay="RNA", tlog=F, tlabel=T, tcells=mock.infected.cells, tdowncolor="#2166ac", tupcolor="#b2182b")
ggsave(file.path(figdir, 'scatter.score.chemokine_genes.SARS-CoV-2_infected.vs.Mock.all_cells.png'), plot=plot, width=5, height=5, dpi=300)

# Figure 3C
# DE between Infected v.s. Mock on beta cells (combine both donors)
tinfo <- FetchData(panc, vars=c('ident','Condition'))
tinfo$DE.cate <- paste(paste0('C',tinfo$ident), tinfo$Condition, sep='-')
# add new labels to seurat object
DE.cate <- tinfo$DE.cate
names(DE.cate) <- rownames(tinfo)
Idents(panc) <- DE.cate
# perform DE analysis
tmarkers <- FindMarkers(panc, ident.1='C2-Infected', ident.2='C2-Mock', logfc.threshold=0.1, min.pct=0.1, test.use='wilcox')
write.table(as.data.frame(tmarkers), file=file.path(infodir, "DE.Infected.vs.Mock.all.beta_cells.wilcox.minpct_0.1.logfc.threshold_0.1.txt"), quote=FALSE, na="NA", sep="\t", col.names=NA)
# Volcano plot
plot <- myVolcanoPlot(file.path(infodir, "DE.Infected.vs.Mock.all.beta_cells.wilcox.minpct_0.1.logfc.threshold_0.1.txt"), tx='p_val', ty='avg_logFC', tcutFC=0.15, tcutP=15, tlabel.genes=NULL, txlabel='LogFoldChange', tylabel='-LogP-value', tupper=150, talpha=0.8, tcolor.up='firebrick3', tcolor.down='steelblue3', tcolor.other='gray60')
ggsave(file.path(figdir, "volcano.DE.Infected.vs.Mock.all.beta_cells.wilcox.minpct_0.1.logfc.threshold_0.1.png"), plot=plot, width=7, height=6, dpi=300)
# assign original clusters back
Idents(panc) <- panc$final.clust

# Figure S4D
# Bar plot illustrating mean expression difference of PPSS2 between mock and SARS-CoV-2 infected cells in different clusters of human islets
plots <- myMeanExpBar(tobj=panc, tgene='PRSS2', ecut=1.7, tfunc=mean, tclust_by='ident', tclust_order=0:8, tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=mock.infected.cells, tlowcolor="royalblue3", thighcolor="orangered3")
ggsave(file.path(figdir, paste('Bar.exp.diff.percent.mean_1.7','PRSS2','SARS-CoV-2_infected.vs.Mock.across_clusters.png',sep='.')), plot=plots[[2]], width=6, height=5, dpi=300)

# collect beta cells from infected v.s. mock (merge two donors) only
mock.infected.beta.cells <- rownames(subset(FetchData(panc, vars=c('ident', 'Condition')), ident %in% c(2) & Condition %in% c('Mock','Infected')))

# Figure 3E
# Dot plot illustrating gene expression levels involved in interferon signaling pathway in mock versus SARS-CoV-2 infected human islets at 24 hpi. 
isg.markers <- c('ADAR','BST2','CGAS','CD74','DDIT4','DDX58','DDX60','EIF2AK2','GBP2','IFI6','IFIH1','IFIT1','IFIT2','IFIT3','IFIT5',
                 'IFITM3','IRF7','ISG15','ISG20','MAP3K14','MOV10','MX1','MX2','NAMPT','NT5C3A','OAS3','JADE2','PML','RSAD2',
                 'RTP4','TRIM5','SUN2','ZC3HAV1')
plot <- myDotPlot.2(tobj=panc, tgenes=isg.markers, tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=mock.infected.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.isg_markers.SARS-CoV-2_infected.vs.Mock.beta_cells.png'), plot=plot, width=16, height=4, dpi=300)

# Figure 4A
# Dot plot illustrating expression level of INS in mock versus SARS-CoV-2 infected human islets at 24 hpi (MOI=1).
plot <- myDotPlot.2(tobj=panc, tgenes=c('INS'), tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=mock.infected.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.beta_markers.SARS-CoV-2_infected.vs.Mock.beta_cells.png'), plot=plot, width=4, height=4, dpi=300)

# Figure 4B
# Dot plot illustrating expression level of alpha cell markers, including GCG, KLHL41, RFX6, SMARCA1, TM4SF4, and RGS4 in mock versus SARS-CoV-2 infected human islets at 24 hpi (MOI=1).
plot <- myDotPlot.2(tobj=panc, tgenes=c('GCG','KLHL41','RFX6','SMARCA1','TM4SF4','RGS4'), tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=mock.infected.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.alpha_markers.SARS-CoV-2_infected.vs.Mock.beta_cells.png'), plot=plot, width=6, height=4, dpi=300)

# Figure 4C
# Dot plot illustrating expression level of acinar cell markers, including PRSS1, PRSS2, CPA1, CPA2, CPB1, SPINK1, and OLFM4 in mock versus SARS-CoV-2 infected human islets at 24 hpi (MOI=1).
plot <- myDotPlot.2(tobj=panc, tgenes=c('PRSS1','PRSS2','CPA1','CPA2','CPB1','SPINK1','OLFM4'), tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=mock.infected.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.acinar_markers.SARS-CoV-2_infected.vs.Mock.beta_cells.png'), plot=plot, width=6.5, height=4, dpi=300)

# Figure 5I
# Dot plot illustrating expression of cell stress associated genes in mock versus SARS-CoV-2 infected human islets at 24 hpi (MOI=1).
# read in genes
stress.genes <- as.character(read.table(file.path(sourcedir, "genesets", "stress_genes.txt"), header=F, check.names=F, stringsAsFactors=F)$V1)
plot <- myDotPlot.2(tobj=panc, tgenes=stress.genes, tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=mock.infected.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.stress_genes.SARS-CoV-2_infected.vs.Mock.beta_cells.png'), plot=plot, width=28, height=4, dpi=300)

# Figure S4E
# Dot plot illustrating expression level of acinar cell markers, including PRSS1, PRSS2, CELA3A, CELA3B, CELA2A in beta cell cluster of mock versus SARS-CoV-2 infected human islets at 24 hpi (MOI=1).
plot <- myDotPlot.2(tobj=panc, tgenes=c('PRSS1','PRSS2','CELA3A','CELA3B','CELA2A'), tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=mock.infected.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.acinar_markers.2.SARS-CoV-2_infected.vs.Mock.beta_cells.png'), plot=plot, width=6, height=4, dpi=300)

# Figure S7H
# Dot plot illustrating expression level of ALDH1A3 in mock versus SARS-CoV-2 infected human islets at 24 hpi (MOI=1)
plot <- myDotPlot.2(tobj=panc, tgenes='ALDH1A3', tgroup_by='Condition', tgroup_order=c('Mock','Infected'), tassay="RNA", tcells=mock.infected.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.ALDH1A3.SARS-CoV-2_infected.vs.Mock.beta_cells.png'), plot=plot, width=3, height=4, dpi=300)

# collect beta cells from DMSO v.s. treated (Healthy_8 donor) only
infected.treated.beta.cells <- rownames(subset(FetchData(panc, vars=c('ident', 'Name')), ident %in% c(2) & Name %in% c('Healthy_8_SARS-CoV-2','Healthy_8_SARS-CoV-2_ISIRB')))

# Figure 7A
# Dot plot illustrating expression level of INS in beta cells of control or 10 μM trans-ISRIB treated human islets at 24 hpi (MOI=1). 
plot <- myDotPlot.2(tobj=panc, tgenes=c('INS'), tgroup_by='Name', tgroup_order=c('Healthy_8_SARS-CoV-2','Healthy_8_SARS-CoV-2_ISIRB'), tassay="RNA", tcells=infected.treated.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.beta_markers.Healthy_8_SARS-CoV-2_ISIRB.vs.Healthy_8_SARS-CoV-2.beta_cells.png'), plot=plot, width=4, height=4, dpi=300)

# Figure 7B
# Dot plot illustrating expression level of alpha cell markers, including GCG, KLHL41, RFX6, SMARCA1, TM4SF4, and RGS4, in beta cells of control or 10 μM trans-ISRIB treated human islets at 24 hpi (MOI=1). 
plot <- myDotPlot.2(tobj=panc, tgenes=c('GCG','KLHL41','RFX6','SMARCA1','TM4SF4','RGS4'), tgroup_by='Name', tgroup_order=c('Healthy_8_SARS-CoV-2','Healthy_8_SARS-CoV-2_ISIRB'), tassay="RNA", tcells=infected.treated.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.alpha_markers.Healthy_8_SARS-CoV-2_ISIRB.vs.Healthy_8_SARS-CoV-2.beta_cells.png'), plot=plot, width=6, height=4, dpi=300)

# Figure 7C
# Dot plot illustrating expression level of acinar cell markers, including PRSS1, PRSS2, CPA1, CPA2, CPB1, SPINK1, and OLFM4, in beta cells of control or 10 μM trans-ISRIB treated human islets at 24 hpi (MOI=1). 
plot <- myDotPlot.2(tobj=panc, tgenes=c('PRSS1','PRSS2','CPA1','CPA2','CPB1','SPINK1','OLFM4'), tgroup_by='Name', tgroup_order=c('Healthy_8_SARS-CoV-2','Healthy_8_SARS-CoV-2_ISIRB'), tassay="RNA", tcells=infected.treated.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.acinar_markers.Healthy_8_SARS-CoV-2_ISIRB.vs.Healthy_8_SARS-CoV-2.beta_cells.png'), plot=plot, width=6.5, height=4, dpi=300)

# Figure 7O
# Dot plot illustrating expression levels of cell stress associated genes of control or 10 μM trans-ISRIB treated human islets at 24 hpi (MOI=1).
plot <- myDotPlot.2(tobj=panc, tgenes=stress.genes, tgroup_by='Name', tgroup_order=c('Healthy_8_SARS-CoV-2','Healthy_8_SARS-CoV-2_ISIRB'), tassay="RNA", tcells=infected.treated.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.stress_genes.Healthy_8_SARS-CoV-2_ISIRB.vs.Healthy_8_SARS-CoV-2.beta_cells.png'), plot=plot, width=28, height=4, dpi=300)

# Figure S7I
# Dot plot illustrating expression level of ALDH1A3 in control or 10 uM trans-ISRIB treated human islets at 24 hpi (MOI=1)
plot <- myDotPlot.2(tobj=panc, tgenes='ALDH1A3', tgroup_by='Name', tgroup_order=c('Healthy_8_SARS-CoV-2','Healthy_8_SARS-CoV-2_ISIRB'), tassay="RNA", tcells=infected.treated.beta.cells, tlowcolor="#2166ac", tmidcolor="#f7f7f7", thighcolor="#b2182b")
ggsave(file.path(figdir, 'Dot.exp.ALDH1A3.Healthy_8_SARS-CoV-2_ISIRB.vs.Healthy_8_SARS-CoV-2.beta_cells.png'), plot=plot, width=4.5, height=4, dpi=300)
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

sessionInfo()
