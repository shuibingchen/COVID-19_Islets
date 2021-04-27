# my_functions.R
# customized functions for processing data and plotting with Seurat3
# Author: Tuo Zhang
# Date: 10/09/2019
# Version: 1.0
# 

library(Seurat, quietly=T, warn.conflicts=F)
library(dplyr, quietly=T, warn.conflicts=F)
library(ggplot2, quietly=T, warn.conflicts=F)
library(Matrix, quietly=T, warn.conflicts=F)
library(scater, quietly=T, warn.conflicts=F)
library(reshape2, quietly=T, warn.conflicts=F)
library(pheatmap, quietly=T, warn.conflicts=F)
library(magrittr, quietly=T, warn.conflicts=F)
library(cowplot, quietly=T, warn.conflicts=F)


# read in raw counts data from a given sample
my.Read10X <- function(tcountdir, tprefix=NULL){
  # directory exists?
  if (! dir.exists(tcountdir)){
    print(paste("input raw counts folder does NOT exist:", tcountdir, sep=" "))
    return(NULL)
  }
  # file exists?
  tmat.file <- paste(tcountdir, "matrix.mtx.gz", sep="/")
  tfnames.file <- paste(tcountdir, "features.tsv.gz", sep="/")
  tbnames.file <- paste(tcountdir, "barcodes.tsv.gz", sep="/")
  for (tf in c(tmat.file, tfnames.file, tbnames.file)){
    if (! file.exists(tf)){
      print(paste("input file does NOT exist:", tf))
      return(NULL)
    }
  }
  # extract counts matrix
  cat(paste("Loading UMI counts table from", tcountdir, "..."))
  tmat <- readMM(paste(tcountdir, "matrix.mtx.gz", sep="/"))
  tfnames <- read.delim(paste(tcountdir, "features.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  tbnames <- read.delim(paste(tcountdir, "barcodes.tsv.gz", sep="/"), header=FALSE, stringsAsFactors=FALSE)
  # update column names (cell ids)
  if (is.null(tprefix)){
    colnames(tmat) = tbnames$V1
  }
  else{
    colnames(tmat) = paste(tprefix, tbnames$V1, sep="_")
  }
  # replace rowname (Ensembl id) by gene symbol
  # in case gene symbol is not unique, append the _EnsemblID after it
  # missing gene symbol will be replaced by EnsemblID
  rownames(tmat) <- uniquifyFeatureNames(ID=tfnames$V1, names=tfnames$V2)
  cat(" done.","\n")
  return(tmat)
}

# merge raw read counts table collected from multiple samples
my.MergeMatrix <- function(tmats){
  cat("Merge raw UMI counts ")
  tfunc <- function(x,y){
    if (class(x) != "data.frame")
      x <- as.data.frame(as.matrix(x))
    if (class(y) != "data.frame")
      y <- as.data.frame(as.matrix(y))
    tres <- merge(x, y, by=0, all=T)
    rownames(tres) <- tres$Row.names
    tres <- tres[,-1]
    cat(".")
    return(tres)
  }
  tmerged <- Reduce(f=tfunc, x=tmats)
  # fill na with 0
  tmerged[is.na(tmerged)] <- 0
  # convert to a compact matrix
  tmerged <- as(as.matrix(tmerged), "dgCMatrix")
  cat(" done.")
  return(tmerged)
}

# running original DE test, without considering conservations in each sample
myDETest <- function(tobj, tk, tmethod, tplot, tinfodir, tfigdir, tassay=NULL, tsuf=NULL){
  # DE test
  tmarkers <- FindMarkers(object=tobj, ident.1=tk, min.pct=0.25, test.use=tmethod, assay=tassay)
  # file name
  tdesp <- ".origin.markers"
  if (! is.null(tsuf)){
    tdesp <- paste("",tsuf,"origin","markers",sep=".")
  }
  # write marker genes to file
  write.table(as.data.frame(tmarkers), file=paste(tinfodir, paste("C", tk, ".", tmethod, tdesp, ".txt", sep=""), sep="/"), quote=FALSE, na="", sep="\t", col.names=NA)
  # select top 4 positive markers
  genes.viz <- head(rownames(subset(tmarkers, avg_logFC>0)),4)
  print(genes.viz)
  # visualize markers with a violin plot
  if(tplot){
    print("violin plot 1")
    png(file=paste(tfigdir, paste("VlnPlot.C", tk, ".", tmethod, tdesp, ".top4.png", sep=""), sep="/"), width=6400, height=6400, units="px", pointsize=6, res=600)
    print(VlnPlot(object=tobj, features=genes.viz, ncol=2))
    dev.off()
    # view which cells express a given gene (red is high expression) on a UMAP plot
    print("UMAP plot")
    png(file=paste(tfigdir, paste("gene.exp.UMAP.C", tk, ".", tmethod, tdesp, ".top4.png", sep=""),sep="/"), width=4000, height=3200, units="px", pointsize=6, res=600)
    print(FeaturePlot(object=tobj, features=genes.viz, reduction="umap", pt.size=1, cols=c("grey","red")))
    dev.off()
  }
}

# plot cells colored by an annotation (DimPlot)
myDimPlot <- function(tobj, treduct, tcate, torder=NULL, tsuffix, tcells=NULL, tcolor=NULL, tlabel=FALSE, tsplit=FALSE, txlim=NULL, tylim=NULL,
                      tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20){
  tdataToPlot <- data.frame()
  tlabel.pos <- data.frame()
  tg <- ggplot()
  if (treduct == "tsne"){
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("tSNE_1","tSNE_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("tSNE_1","tSNE_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      tdataToPlot <- rbind(tdataToPlot, data.frame(tSNE_1=NA, tSNE_2=NA, Category=torder))
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=tSNE_1, y=tSNE_2, color=Category))
    tg <- tg + ggtitle(paste("tSNE","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(tSNE_1, tSNE_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  } else if (treduct == "umap") {
    if (is.null(tcells)){
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate))
    } else {
      tdataToPlot <- FetchData(tobj, vars=c("UMAP_1","UMAP_2",tcate), cells=tcells)
    }
    colnames(tdataToPlot) <- c("UMAP_1","UMAP_2","Category")
    if (!is.null(torder)){
      # add fake rows to make sure each category is considered
      tdataToPlot <- rbind(tdataToPlot, data.frame(UMAP_1=NA, UMAP_2=NA, Category=torder))
      # reorder categories
      tdataToPlot$Category <- factor(tdataToPlot$Category, levels=torder)
    }
    tg <- ggplot(tdataToPlot, aes(x=UMAP_1, y=UMAP_2, color=Category))
    tg <- tg + ggtitle(paste("UMAP","plots","by",tsuffix, sep=" "))
    tlabel.pos <- aggregate(cbind(UMAP_1, UMAP_2) ~ Category, data=tdataToPlot, FUN=median)
    colnames(tlabel.pos) <- c("Category","X","Y")
  }
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (! is.null(tcolor)){
    tg <- tg + scale_color_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (tsplit == TRUE){
    tg <- tg + facet_wrap(~Category)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Category), color="black")
  }
  return(tg)
}

# make heatmap plot
myHeatmap <- function(tobj, tgenes, tcells=NULL, toutfile, tassay="RNA", tgroup_by='ident', tgroup_order=NULL, tmax=3, tmin=-3, tcolor=colorRampPalette(c("royalblue3", "white", "orangered3"))(100), tanncolor=NULL, twidth=12800, theight=16000, tunits="px", tfontsize_row=8, tres=600){
  # valid genes?
  tgenes <- intersect(tgenes, rownames(tobj))
  if (length(tgenes) > 0){
    print(paste(length(tgenes),'genes to plot.',sep=' '))
  } else {
    print('No valid genes found.')
    return(NULL)
  }
  # prepare cell annotations
  anncols <- FetchData(tobj, cells=tcells, vars=tgroup_by)
  if (! is.null(tgroup_order)){
    anncols[,tgroup_by] <- factor(anncols[,tgroup_by], levels=tgroup_order)
  }
  colnames(anncols) <- c('group')
  anncols <- anncols[with(anncols, order(group)),,drop=F]
  # split columns by group
  tgaps_col <- cumsum(as.vector(table(anncols$group)))
  # prepare matrix for pheatmap
  texp <- tobj[[tassay]]@scale.data[tgenes, rownames(anncols)]
  # set minimum/maximum display value (clip all values that are out of range)
  texp <- apply(texp, 2, function(x) ifelse(x > tmax, tmax, x))
  texp <- apply(texp, 2, function(x) ifelse(x < tmin, tmin, x))
  # plot
  png(file=toutfile, width=twidth, height=theight, units=tunits, pointsize=8, res=tres)
  if(is.null(tanncolor)){
    print(pheatmap(texp, scale="none", color=tcolor, cluster_rows=F, cluster_cols=F, annotation_col=anncols, show_colnames=F, show_rownames=T, gaps_col=tgaps_col, fontsize_row=tfontsize_row))
  } else {
    annColors = list(group=tanncolor)
    print(pheatmap(texp, annotation_color=annColors, scale="none", color=tcolor, cluster_rows=F, cluster_cols=F, annotation_col=anncols, show_colnames=F, show_rownames=T, gaps_col=tgaps_col, fontsize_row=tfontsize_row))
  }
  dev.off()
}

# highlight expression of a set of genes (FeaturePlot)
MyFeaturePlot <- function(tobj, tgenes, tcells=NULL, tassay="RNA", treduction.name="umap", txlim=NULL, tylim=NULL, tbreaks=NULL, tlimits=NULL, tlowcolor='gray80', thighcolor='red2', tncol=2, tlegend=NULL, tptsize=2, talpha=0.7){
  # genes valid?
  tgenes.valid <- intersect(tgenes, rownames(tobj))
  if (is.null(tgenes.valid)){
    cat("No valid genes found, do nothing!")
    return(NULL)
  }
  # assay valid?
  if (! tassay %in% names(tobj)){
    cat(paste("Not a valid assay:",tassay,sep=" "))
    return(NULL)
  }
  # extract gene expression
  texp <- as.matrix(tobj[[tassay]]@data[tgenes.valid, ,drop=F])
  # get coordinates
  tvars <- c("UMAP_1","UMAP_2")
  if (treduction.name == "tsne"){
    tvars <- c("tSNE_1","tSNE_2")
  }
  tdata <- FetchData(object=tobj, vars=tvars)
  colnames(tdata) <- c("X","Y")
  # plot
  tplots <- list()
  tk <- 1
  for (tgene in tgenes){
    # merge data for plotting
    tdata.merged <- merge(tdata, t(texp[tgene,,drop=F]), by=0, all=T)
    rownames(tdata.merged) <- tdata.merged$Row.names
    tdata.merged <- tdata.merged[,-1]
    colnames(tdata.merged) <- c("X","Y","Expression")
    # subset cells?
    if (! is.null(tcells)){
      tdata.merged <- tdata.merged[tcells,,drop=F]
    }
    # reorder cells by expression of the given gene
    tdata.merged <- tdata.merged[with(tdata.merged, order(Expression)),]
    # plot
    if (max(tdata.merged$Expression) > 0){ # expressed in at least one cell
      # plot (rename x and y axis)
      tg <- ggplot(tdata.merged, aes(x=X, y=Y, color=Expression))
      tg <- tg + geom_point(shape=19, size=tptsize, alpha=talpha)
      if (! is.null(tbreaks)){
        tg <- tg + scale_color_gradient(low=tlowcolor, high=thighcolor, breaks=tbreaks, limits=tlimits)
      } else {
        tg <- tg + scale_color_gradient(low=tlowcolor, high=thighcolor)
      }
      if(! is.null(txlim)){
        tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
      }
      tg <- tg + ggtitle(tgene)
      if (treduction.name == "tsne"){
        tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
      } else {
        tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
      }
      tg <- tg + theme_bw()
      tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
      tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
      # add to list
      tplots[[tk]] <- tg
      tk <- tk + 1
    } else { # no expressions at all
      tg <- ggplot(tdata.merged, aes(x=X, y=Y))
      tg <- tg + geom_point(color="gray80", shape=19, size=tptsize, alpha=talpha)
      if(! is.null(txlim)){
        tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
      }
      tg <- tg + ggtitle(tgene)
      if (treduction.name == "tsne"){
        tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
      } else {
        tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
      }
      tg <- tg + theme_bw()
      tg <- tg + theme(plot.title=element_text(hjust=0.5, size=18, face="bold"))
      tg <- tg + theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_text(size=16,face="bold"))
      # add to list
      tplots[[tk]] <- tg
      tk <- tk + 1
    }
  }
  # combine plots with Seurat::CombinePlots
  tcombined <- CombinePlots(tplots, ncol=tncol, legend=tlegend)
  return(tcombined)
}

# draw violin plot highlighting expression of a gene in selected cells in each selected group
MyExpViolin <- function(tobj, tgene, tgroup_by, tgroup_order=NULL, tsplit_by=NULL, tsplit_order=NULL, tcolor_by=NULL, tcolor_order=NULL, tcolor=NULL, tcells=NULL, tassay="RNA", tncol=2){
  # gene valid?
  if (! tgene %in% rownames(tobj[[tassay]]@data)){
    cat(paste("Gene",tgene,"unavailable.\n",sep=" "))
    return(NULL)
  }
  # assign cells
  if (is.null(tcells)){
    tcells <- colnames(tobj[[tassay]]@data)
  }
  # extract expression of selected genes
  tdata <- data.frame(gene = as.vector(as.matrix(tobj[[tassay]]@data)[tgene, tcells]))
  rownames(tdata) <- tcells
  # add color info
  if (! is.null(tcolor_by)){
    tinfo <- FetchData(tobj, vars=c(tcolor_by), cells=tcells)
    tdata$color <- as.vector(tinfo[tcells,tcolor_by])
  }
  if (! is.null(tcolor_order)){
    tdata$color <- factor(tdata$color, levels=tcolor_order)
  } else {
    tcolor_order <- unique(tdata$color)
  }
  # add split info
  if (! is.null(tsplit_by)){
    tinfo <- FetchData(tobj, vars=c(tsplit_by), cells=tcells)
    tdata$split <- as.vector(tinfo[tcells,tsplit_by])
    if (! is.null(tsplit_order)){
      tdata$split <- factor(tdata$split, levels=tsplit_order)
    }
  }
  # add group info
  tinfo <- FetchData(tobj, vars=c(tgroup_by), cells=tcells)
  tdata$group <- as.vector(tinfo[tcells,tgroup_by])
  # missing groups in any sample?
  if (! is.null(tsplit_by)){
    all.groups <- unique(as.vector(tdata$group))
    for (tsmp in unique(tdata$split)){
      tmiss.groups <- setdiff(all.groups, unique(as.vector(subset(tdata, split %in% c(tsmp))[,'group'])))
      if (length(tmiss.groups) > 0){
        # add fake rows to make sure all groups are considered
        tdummy <- data.frame(gene=NA, color=NA, split=tsmp, group=tmiss.groups)
        tdata <- rbind(tdata, tdummy[, colnames(tdata)])
      }
    }
  }
  # reorder groups
  if (! is.null(tgroup_order)){
    # reorder groups
    tdata$group <- factor(tdata$group, levels=tgroup_order)
  }
  # plot
  g <- ggplot(tdata, aes(x=group, y=gene))
  if (is.null(tcolor_by)){
    g <- g + geom_violin(scale="width", alpha=0.7)
    g <- g + geom_jitter(shape=16, position=position_jitter(width=0.2, height=0), alpha=0.75)
  } else {
    g <- g + geom_violin(aes(fill=color), scale="width", alpha=0.7)
    g <- g + geom_jitter(aes(color=color), shape=16, position=position_jitter(width=0.2, height=0), alpha=0.75)
  }
  if (! is.null(tsplit_by)){
    g <- g + facet_wrap(~split, ncol=tncol)
  }
  if (! is.null(tcolor)){
    g <- g + scale_color_manual(values=tcolor)
    g <- g + scale_fill_manual(values=tcolor)
  }
  g <- g + stat_summary(fun.data=median_hilow, geom="pointrange", color="gray75")
  g <- g + ggtitle(tgene)
  g <- g + theme_bw()
  g <- g + theme(legend.position="none")
  g <- g + theme(axis.title=element_blank(), axis.text.x=element_text(size=18, angle=45, hjust=1, face="bold"))
  g <- g + theme(axis.text.y=element_text(size=18, face="bold"), plot.title=element_text(size=20, hjust=0.5, face="bold"))
  return(g)
}

# scatter plot showing a set of genes in two different conditions
myScatterCdt <- function(tobj, tgenes, tgroup_by, tgroup_x, tgroup_y, tassay="RNA", tlog=FALSE, tmin, tmax, tlabel=FALSE, tcells=NULL, tdowncolor="royalblue3", tupcolor="orangered3"){
  # check genes
  tgenes <- intersect(tgenes, rownames(tobj))
  if (length(tgenes) == 0){
    print('No valid genes found.')
    return(NULL)
  }
  # get cluster info
  tdata <- FetchData(tobj, cells=tcells, vars=c(tgroup_by))
  colnames(tdata) <- c('group')
  # add expression
  if (! is.null(tcells)){
    texp <- as.data.frame(t(as.matrix(tobj[[tassay]]@data[tgenes, tcells, drop=F])))
    tdata <- merge(tdata, texp, by=0, all=T)
  } else {
    texp <- as.data.frame(t(as.matrix(tobj[[tassay]]@data[tgenes, ,drop=F])))
    tdata <- merge(tdata, texp, by=0, all=T)
  }
  rownames(tdata) <- tdata$Row.names
  tdata <- tdata[,-1]
  # define functions: Percent * MeanExp
  calScore <- function(x){
    tpercent <- sum(x>0)/length(x)*100
    tmeanexp <- 0
    if (sum(x>0) > 0){
      tmeanexp <- mean(x[x>0])
    }
    return(tpercent*tmeanexp)
  }
  # summary data for plotting
  tscore <- as.data.frame(tdata %>% group_by(group) %>% summarise_all(calScore))
  rownames(tscore) <- tscore$group
  tscore <- tscore[, -which(colnames(tscore)=='group')]
  tscore.forPlot <- as.data.frame(t(tscore))
  
  # add color
  tscore.forPlot$color <- ifelse(tscore.forPlot[,tgroup_y] > tscore.forPlot[,tgroup_x], 'up', 'down')
  # add label
  tscore.forPlot$label <- rownames(tscore.forPlot)
  # max value
  tscore.max <- max(tscore)
  tup <- ceiling(tscore.max)
  if (ceiling(tscore.max) - 0.5 > tscore.max){
    tup <- ceiling(tscore.max) - 0.5
  }
  # make scatter plot
  g <- ggplot(tscore.forPlot, aes_string(x=tgroup_x, y=tgroup_y, color='color', label='label'))
  g <- g + geom_abline(intercept=0, slope=1, color='gray75', linetype='dotted')
  g <- g + geom_point(shape=19)
  if (tlabel){
    g <- g + geom_text_repel()
  }
  g <- g + scale_color_manual(values=c('down'=tdowncolor, 'up'=tupcolor))
  if (tlog){
    g <- g + coord_fixed(xlim=c(tmin, tmax), ylim=c(tmin, tmax))
    g <- g + scale_x_log10(labels = function(x) format(x, scientific = FALSE))
    g <- g + scale_y_log10(labels = function(x) format(x, scientific = FALSE))
  } else {
    g <- g + coord_fixed(xlim=c(0, tup), ylim=c(0, tup))
  }
  g <- g + theme_classic()
  g <- g + theme(axis.line=element_blank())
  g <- g + theme(panel.border=element_rect(color="black", fill=NA, size=1))
  g <- g + theme(axis.text=element_text(color='black', size=14), axis.title=element_text(color='black', size=16))
  g <- g + theme(legend.position='none')
  return(g)
}

# make volcano plot on a set of DE genes
myVolcanoPlot <- function(tdefile, tx='p_val', ty='avg_logFC', tcutFC=0.25, tcutP=20, tlabel.genes=NULL, txlabel='LogFoldChange', tylabel='-LogP-value', tupper=300, talpha=0.8, tcolor.up='firebrick3', tcolor.down='steelblue3', tcolor.other='gray60'){
  # read in DE results
  tde <- read.table(tdefile, sep='\t', header=T, check.names=F, stringsAsFactors=F, row.names=1)
  # arrange data
  tdataToPlot <- tde[,c(tx,ty)]
  tdataToPlot$gene <- rownames(tde)
  colnames(tdataToPlot) <- c('pval','logFC','gene')
  # log transform p-value
  tdataToPlot$logPval <- -log10(tdataToPlot$pval)
  # color
  tdataToPlot$color <- with(tdataToPlot, ifelse(logPval > tcutP & logFC > tcutFC, 'up', ifelse(logPval > tcutP & logFC < -tcutFC, 'down', 'other')))
  # label
  if (is.null(tlabel.genes)){
    tdataToPlot$label <- with(tdataToPlot, ifelse(color %in% c('up','down'), gene, ''))
  } else{
    tdataToPlot$label <- with(tdataToPlot, ifelse(gene %in% tlabel.genes, gene, ''))
  }
  # any 0 p-values thus Inf logPval? modify to the upperlimit
  tdataToPlot$logPval[tdataToPlot$logPval > tupper] <- tupper
  # plot
  tg <- ggplot(tdataToPlot, aes(x=logFC, y=logPval, color=color, label=label))
  tg <- tg + geom_point(shape=19, size=2, alpha=talpha)
  tg <- tg + scale_color_manual(values=c('up'=tcolor.up,'down'=tcolor.down, 'other'=tcolor.other))
  tg <- tg + geom_text_repel()
  tg <- tg + xlab(txlabel) + ylab(tylabel)
  tg <- tg + theme_classic()
  tg <- tg + theme(legend.position='none')
  tg <- tg + theme(axis.title=element_text(size=18, color='black'), axis.text=element_text(size=16, color='black'))
  return(tg)
}

# normalize by mean
normalizeByMean <- function(tvals){
  # normalize a vector of two elements a,b: 2a/(a+b)-1, 2b/(a+b)-1 (i.e. a-b/a+b, b-a/a+b)
  tmean = mean(tvals)
  if (tmean != 0){
    return(tvals/tmean-1)
  } else {
    return(tvals)
  }
}

# normalize by z-score
normalizeByZScore <- function(tvals){
  if (sd(tvals) == 0){
    return(rep(0, length(tvals)))
  } else {
    return(scale(tvals)[,1])
  }
}

# dot plot (recolor dot after normalize expression across samples)
myDotPlot.2 <- function(tobj, tgenes, tgroup_by, tgroup_order=NULL, tassay="RNA", tcells=NULL, tlowcolor="royalblue3", tmidcolor="white", thighcolor="orangered3"){
  # check genes
  tgenes <- intersect(tgenes, rownames(tobj))
  if (length(tgenes) == 0){
    print('No valid genes found.')
    return(NULL)
  }
  # get cluster info
  tdata <- FetchData(tobj, cells=tcells, vars=c(tgroup_by))
  colnames(tdata) <- c('group')
  # add expression
  if (! is.null(tcells)){
    texp <- as.data.frame(t(as.matrix(tobj[[tassay]]@data[tgenes, tcells, drop=F])))
    tdata <- merge(tdata, texp, by=0, all=T)
  } else {
    texp <- as.data.frame(t(as.matrix(tobj[[tassay]]@data[tgenes, ,drop=F])))
    tdata <- merge(tdata, texp, by=0, all=T)
  }
  rownames(tdata) <- tdata$Row.names
  tdata <- tdata[,-1]
  # define functions
  calPercent <- function(x){
    return(sum(x>0)/length(x)*100)
  }
  calMeanExp <- function(x){
    tmeanexp <- 0
    if (sum(x>0) > 0){
      tmeanexp <- mean(x[x>0])
      ##tmeanexp <- log1p(mean(expm1(x[x>0])))
    }
    return(tmeanexp)
  }
  # summary data for plotting
  # percent of cells expressing a gene
  tdata.percent <- melt(as.data.frame(tdata %>% group_by(group) %>% summarise_all(calPercent)))
  colnames(tdata.percent) <- c('Group','Gene','Percent')
  # mean expression on expressed cells
  tdata.meanexp <- melt(as.data.frame(tdata %>% group_by(group) %>% summarise_all(calMeanExp)))
  colnames(tdata.meanexp) <- c('Group','Gene','MeanExp')
  # merge
  tdata.forDotPlot <- merge(tdata.percent, tdata.meanexp, by=c('Group','Gene'), all=T)

  # normalize mean expression across samples
  tdata.forDotPlot.v2 <- NULL
  if (length(unique(tdata.forDotPlot$Group)) > 2){
    tdata.forDotPlot.v2 <- as.data.frame(tdata.forDotPlot %>% group_by(Gene) %>% mutate(NormMeanExp=normalizeByZScore(MeanExp)))
  } else {# only two groups, use mean normalization
    tdata.forDotPlot.v2 <- as.data.frame(tdata.forDotPlot %>% group_by(Gene) %>% mutate(NormMeanExp=normalizeByMean(MeanExp)))
  }
  
  # set group order
  if (! is.null(tgroup_order)){
    tdata.forDotPlot.v2$Group <- factor(tdata.forDotPlot.v2$Group, levels=tgroup_order)
  }

  # make expression dot plot 
  g <- ggplot(tdata.forDotPlot.v2, aes(x=Gene, y=Group, color=NormMeanExp, size=Percent))
  g <- g + geom_point(shape=19)
  g <- g + scale_color_gradient2(low=tlowcolor,mid=tmidcolor, high=thighcolor)
  g <- g + scale_x_discrete(position="top")
  g <- g + scale_radius(breaks=c(0,25,50,75,100), labels=c('0','25','50','75','100'), limits=c(0,100), range=c(0,10))
  g <- g + theme_classic()
  g <- g + theme(axis.line=element_blank())
  g <- g + theme(panel.border=element_rect(color="black", fill=NA, size=1))
  g <- g + theme(axis.text.x=element_text(angle=45, hjust=0, face="bold"), axis.title.x=element_blank())
  g <- g + theme(axis.text.y=element_text(face="bold"), axis.title.y=element_blank())
  return(g)
}

# bar plot highlighting changes in mean expression per cluster among conditions
# borrow codes for function myDotPlot.3
myMeanExpBar <- function(tobj, tgene, ecut=0, tfunc=mean, tclust_by='ident', tclust_order=NULL, tgroup_by, tgroup_order, tassay="RNA", tcells=NULL, tlowcolor="royalblue3", thighcolor="orangered3"){
  # check genes
  if (! tgene %in% rownames(panc)){
    print('No a valid gene.')
    return(NULL)
  }
  # double check on group order, has to be two groups only
  if (length(tgroup_order) != 2){
    print('Should be two groups only!')
    print(paste(tgroup_order, collapse=', '))
    return(NULL)
  }
  # get cluster info
  tdata <- FetchData(tobj, cells=tcells, vars=c(tclust_by, tgroup_by))
  colnames(tdata) <- c('cluster','group')
  # add expression
  if (! is.null(tcells)){
    texp <- as.data.frame(t(as.matrix(tobj[[tassay]]@data[tgene, tcells, drop=F])))
    tdata <- merge(tdata, texp, by=0, all=T)
  } else {
    texp <- as.data.frame(t(as.matrix(tobj[[tassay]]@data[tgene, ,drop=F])))
    tdata <- merge(tdata, texp, by=0, all=T)
  }
  rownames(tdata) <- tdata$Row.names
  tdata <- tdata[,-1]
  # define functions
  # mean
  calMeanExp <- function(x){
    tmeanexp <- 0
    if (sum(x >= ecut) > 0){
      tmeanexp <- tfunc(x[x >= ecut])
    }
    return(tmeanexp)
  }
  # summary data for plotting
  # mean expression on expressed cells
  tdata <- as.data.frame(tdata %>% group_by(group, cluster) %>% summarise_at(tgene, list(MeanExp=calMeanExp)))
  colnames(tdata) <- c('Group','Cluster','MeanExp')
  
  # calculate difference
  # get two groups seperately and manually merge them
  grp.1 <- tgroup_order[1]
  grp.2 <- tgroup_order[2]
  tdata.1 <- tdata %>% filter(Group == grp.1) %>% select(-Group) %>% rename(c("MeanExp"=grp.1))
  tdata.2 <- tdata %>% filter(Group == grp.2) %>% select(-Group) %>% rename(c("MeanExp"=grp.2))
  tdata.forPlot <- merge(tdata.1, tdata.2, by='Cluster')
  # mean
  tdata.forPlot$AbsDiff <- tdata.forPlot[,grp.2] - tdata.forPlot[,grp.1]
  tdata.forPlot$PercentDiff <- round((tdata.forPlot[,grp.2] - tdata.forPlot[,grp.1]) / tdata.forPlot[,grp.1] * 100, 2)
  tdata.forPlot$Sign <- ifelse(tdata.forPlot$AbsDiff > 0, "increase", "decrease")
  
  # set cluster order
  if (! is.null(tclust_order)){
    tdata.forPlot$Cluster <- factor(tdata.forPlot$Cluster, levels=tclust_order)
  }
  
  # make absolute expression changes plot 
  plot1 <- ggplot(tdata.forPlot, aes(x=Cluster, y=AbsDiff, fill=Sign))
  plot1 <- plot1 + geom_bar(stat="identity")
  plot1 <- plot1 + geom_hline(yintercept=0)
  plot1 <- plot1 + scale_fill_manual(values=c('decrease'=tlowcolor, 'increase'=thighcolor))
  plot1 <- plot1 + ylab("Mean expression difference (absolute)")
  plot1 <- plot1 + theme_classic()
  plot1 <- plot1 + theme(axis.line=element_blank())
  plot1 <- plot1 + theme(panel.border=element_rect(color="black", fill=NA, size=1))
  plot1 <- plot1 + theme(axis.text.x=element_text(angle=45, hjust=0, face="bold"), axis.title.x=element_blank())
  plot1 <- plot1 + theme(axis.text.y=element_text(face="bold"), legend.position="none")
  # make absolute expression changes plot 
  plot2 <- ggplot(tdata.forPlot, aes(x=Cluster, y=PercentDiff, fill=Sign))
  plot2 <- plot2 + geom_bar(stat="identity")
  plot2 <- plot2 + geom_hline(yintercept=0)
  plot2 <- plot2 + scale_fill_manual(values=c('decrease'=tlowcolor, 'increase'=thighcolor))
  plot2 <- plot2 + ylab("Mean expression difference (%)")
  plot2 <- plot2 + theme_classic()
  plot2 <- plot2 + theme(axis.line=element_blank())
  plot2 <- plot2 + theme(panel.border=element_rect(color="black", fill=NA, size=1))
  ##plot2 <- plot2 + theme(axis.text.x=element_text(angle=45, hjust=0, face="bold"), axis.title.x=element_blank())
  plot2 <- plot2 + theme(axis.text.x=element_text(size=16), axis.title.x=element_blank())
  ##plot2 <- plot2 + theme(axis.text.y=element_text(face="bold"), legend.position="none")
  plot2 <- plot2 + theme(axis.text.y=element_text(size=16), legend.position="none")
  return(list(plot1, plot2))
}

# plot cells colored by an annotation (DimPlot)
# similar to myDimPlot, except that split by another label (e.g. color by ident but split by sample)
# this version is a more generalized function 'myDimPlot'
myDimPlot3 <- function(tobj, treduct, tgroup_by, tgroup_order=NULL, thighlight=NULL, tsuffix, tcells=NULL, tcolor=NULL, tlabel=FALSE, tsplit_by=NULL, tsplit_order=NULL, txlim=NULL, tylim=NULL,
                       tncol=1, tptshape=19, tptsize=2, talpha=0.6, tltsize=18, tatlsize=20){
  # set coordinates variable name
  vars.reduct <- c("UMAP_1","UMAP_2")
  if (treduct == "tsne"){
    vars.reduct <- c("tSNE_1","tSNE_2")
  }
  # extract coordinates + group
  tdataToPlot <- FetchData(tobj, cells=tcells, vars=c(vars.reduct, tgroup_by))
  colnames(tdataToPlot) <- c("Dim_1","Dim_2","Group")
  # update group order if available
  if (!is.null(tgroup_order)){
    # add fake rows to make sure each group is considered
    tmiss.cates <- setdiff(tgroup_order, unique(tdataToPlot$Group))
    if (length(tmiss.cates) > 0){
      tdataToPlot <- rbind(tdataToPlot, data.frame(Dim_1=NA, Dim_2=NA, Group=tmiss.cates))
    }
    # reorder categories
    tdataToPlot$Group <- factor(tdataToPlot$Group, levels=tgroup_order)
  }
  # extract split
  if (! is.null(tsplit_by)){
    tsp <- FetchData(tobj, cells=tcells, vars=c(tsplit_by))
    tdataToPlot$Split <- tsp[rownames(tdataToPlot), c(tsplit_by)]
    # update split order if available
    if (! is.null(tsplit_order)){
      tdataToPlot$Split <- factor(tdataToPlot$Split, levels=tsplit_order)
    }
  }
  # reorder cells that needs to highlight (draw those cell points later)
  if (!is.null(thighlight)){
    tdataToPlot <- rbind(subset(tdataToPlot, ! Group %in% thighlight), subset(tdataToPlot, Group %in% thighlight))
  }
  # prepare group labeling
  tlabel.pos <- aggregate(cbind(Dim_1, Dim_2) ~ Group, data=tdataToPlot, FUN=median)
  colnames(tlabel.pos) <- c("Group","X","Y")
  # plot
  tg <- ggplot(tdataToPlot, aes(x=Dim_1, y=Dim_2, color=Group))
  tg <- tg + ggtitle(paste(toupper(treduct),"plots","by",tsuffix, sep=" "))
  # set range on coordinates
  if (! is.null(txlim)){
    tg <- tg + coord_cartesian(xlim=txlim, ylim=tylim)
  }
  tg <- tg + theme_bw() + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), legend.justification=c(0,0), legend.title=element_blank())
  tg <- tg + theme(legend.key=element_blank()) + theme(legend.text=element_text(size=tltsize))
  tg <- tg + theme(axis.text=element_blank(), axis.title=element_text(size=tatlsize,face="bold"))
  tg <- tg + theme(axis.ticks=element_blank())
  tg <- tg + theme(plot.title=element_text(hjust=0.5))
  tg <- tg + geom_point(shape=tptshape, size=tptsize, alpha=talpha)
  if (! is.null(tcolor)){
    tg <- tg + scale_color_manual(values=tcolor)
  }
  if (treduct == "tsne"){
    tg <- tg + xlab("tSNE 1") + ylab("tSNE 2")
  } else {
    tg <- tg + xlab("UMAP 1") + ylab("UMAP 2")
  }
  if (! is.null(tsplit_by)){
    tg <- tg + facet_wrap(~Split, ncol=tncol)
  }
  if (tlabel == TRUE){
    tg <- tg + geom_text(data=tlabel.pos,aes(x=X, y=Y, label=Group), color="black")
  }
  return(tg)
}
