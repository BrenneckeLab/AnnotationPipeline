args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

# available variables
print(ls())

###################################################################################################

library(gridExtra)
library(tidyverse)
library(reshape2)
library(RColorBrewer)
library(cowplot)
library(plotly)
library(Cairo)
library("tximport")
library("readr")
library("DESeq2")
library(gridExtra)
library(edgeR)
library(pheatmap)
library(Glimma)

theme_set(theme_cowplot())
###################################################################################################
#define Glimma Volcano function with modified 0 pvalue behavior
extractGroups <- function(cdata)
{
  if ("group" %in% colnames(cdata))
  {
    return(cdata[, "group"])
  } else
  {
    return(1)
  }
}


buildXYData <- function(
  table,
  status,
  main,
  display.columns,
  anno,
  counts,
  xlab,
  ylab,
  status.cols,
  sample.cols,
  groups,
  transform.counts)
{
  if (is.null(counts)) {
    counts <- -1
    level <- NULL
  } else {
    # df format for serialisation
    if (transform.counts != "none") {
      if (!all.equal(counts, round(counts))) {
        warning("count transform requested but not all count values are integers.")
      }

      if (transform.counts == "logcpm") {
        counts <- edgeR::cpm(counts, log=TRUE)
      } else if (transform.counts == "cpm") {
        counts <- edgeR::cpm(counts, log=FALSE)
      } else if (transform.counts == "rpkm") {
        if (is.null(anno$length)) {
          stop("no 'length' column in gene annotation, rpkm cannot be computed")
        }

        if (!is.numeric(anno$length)) {
          stop("'length' column of gene annotation must be numeric values")
        }
        counts <- edgeR::rpkm(counts, gene.length = anno$length)
      }
    }

    counts <- data.frame(counts)
    #if (is.null(groups)) stop("If counts arg is supplied, groups arg must be non-null.")
    if (is.null(groups)) {
      groups <- factor("group")
    } else {
      if (ncol(counts) != length(groups)) stop("Length of groups must be equal to the number of columns in counts.\n")
    }

    level <- levels(groups)
    groups <- data.frame(group=groups)
    groups <- cbind(groups, sample=colnames(counts))
  }

  if (length(status)!=nrow(table)) stop("Status vector
     must have the same number of genes as the main arguments.")

  table <- cbind(table, status=as.vector(status))
  if (!is.null(anno))
  {
    colnames(anno) <- gsub("symbol", "symbol", colnames(anno), ignore.case=TRUE)
    table <- cbind(table, anno)
  }

  if (is.null(display.columns)) {
    display.columns <- colnames(table)
  } else {
    # if it's specified, make sure at least x, y, gene are displayed in the table and tooltips
    if (!(xlab %in% display.columns)) display.columns <- c(display.columns, xlab)
    if (!(ylab %in% display.columns)) display.columns <- c(display.columns, ylab)
    if (!("gene" %in% display.columns)) display.columns <- c("gene", display.columns)
  }

  table <- data.frame(index=0:(nrow(table)-1), table)

  if (length(status.cols) != 3) stop("status.cols
          arg must have exactly 3 elements for [downreg, notDE, upreg]")

  xData <- list(data=list(x=xlab,
                          y=ylab,
                          table=table,
                          cols=display.columns,
                          counts=counts,
                          groups=groups,
                          levels=level,
                          expCols=colnames(groups),
                          statusColours=status.cols,
                          sampleColours= if (is.null(sample.cols)) {-1} else {sample.cols},
                          samples=colnames(counts),
                          title=main))
  return(xData)
}

glimmaXYWidget <- function(xData, width, height, html)
{
  widget <- htmlwidgets::createWidget(
    name = 'glimmaXY',
    xData,
    width = width,
    height = height,
    package = 'Glimma',
    elementId = NULL,
    sizingPolicy = htmlwidgets::sizingPolicy(defaultWidth=width, defaultHeight=height, browser.fill=TRUE, viewer.suppress=TRUE)
  )
  if (is.null(html))
  {
    return(widget)
  }
  else
  {
    message("Saving widget...")
    htmlwidgets::saveWidget(widget, file=html)
    message(html, " generated.")
  }
}

glimmaVolcano.DESeqDataSet.pValFixed  <- function(
  x,
  counts=DESeq2::counts(x),
  groups=extractGroups(colData(x)),
  status=NULL,
  anno=NULL,
  display.columns = NULL,
  status.cols=c("dodgerblue", "silver", "firebrick"),
  sample.cols=NULL,
  transform.counts = c("logcpm", "cpm", "rpkm", "none"),
  main="Volcano Plot",
  xlab="logFC",
  ylab="negLog10PValue",
  html=NULL,
  width = 920,
  height = 920,
  ...)
{
  transform.counts <- match.arg(transform.counts)
  res.df <- as.data.frame(DESeq2::results(x))

  # filter out genes that have missing data
  complete_genes <- complete.cases(res.df)
  res.df <- res.df[complete_genes, ]
  x <- x[complete_genes, ]

  # extract status if it is not given
  if (is.null(status))
  {
    status <- ifelse(
      res.df$padj < 0.05,
      ifelse(res.df$log2FoldChange < 0, -1, 1),
      0
    )
  }
  else
  {
    if (length(status)!=length(complete_genes)) stop("Status vector
      must have the same number of genes as the main arguments.")
    status <- status[complete_genes]
  }

  # create initial table with logFC and -log10(pvalue) features
  res.df$pvalue[res.df$pvalue==0]=1e-280
  table <- data.frame(signif(res.df$log2FoldChange, digits=4),
                      signif(-log10(res.df$pvalue), digits=4))
  colnames(table) <- c(xlab, ylab)
  table <- cbind(table, logCPM=signif(log(res.df$baseMean + 0.5), digits=4),
                        AdjPValue=signif(res.df$padj, digits=4))
  table <- cbind(gene=rownames(x), table)
  xData <- buildXYData(table, status, main, display.columns, anno, counts, xlab, ylab, status.cols, sample.cols, groups, transform.counts)
  return(glimmaXYWidget(xData, width, height, html))
}

###################################################################################################
setwd( TMP)

rawFOLDER=FOLDER
FOLDER=paste(FOLDER,"/gene-expression/",sep="")

 if (DGE == "Y"){
   INFO=paste(rawFOLDER, "DGE/DGE_info.txt", sep="")
   TX2GENE <- read.table(file.path(TMP, "transcript_to_gene.txt"), stringsAsFactors = FALSE, header= FALSE)
 }else{
   if (Nexec == 1){
    INFO=paste(TMP, "DGE_info_all.txt", sep="")
   }else{
     INFO=paste(TMP, "DGE_info_all.CDS.txt", sep="")
    FOLDER=paste(FOLDER,"/gene-expression-CDS/",sep="")
   }
   TX2GENE <- read.table(file.path(TMP, "transcript_to_gene.txt"), stringsAsFactors = FALSE, header= FALSE)
 }

###################################################################################################
#read in info
SAMPLES <- read.table(INFO, header=TRUE, as.is=TRUE)
rownames(SAMPLES) <- SAMPLES$sample

FILES <- SAMPLES$path
names(FILES) <- SAMPLES$sample
colnames(TX2GENE)=c("TXNAME","GENEID")



#---------------------------------------------------------------------------------------------------------
#prepare in data 
TXI <- tximport(FILES, type="salmon", tx2gene=TX2GENE)

#prepare data for edgeR
CTS <- TXI$counts
normMat <- TXI$length

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- CTS/normMat

#create edgeR object
y <- DGEList(normCts)

#extract unnormalized log transformed RPK numbers
lcpm <- cpm(y, log=TRUE)

#reate RPK plot
p=ggplot(melt(lcpm),aes(x=Var2, y=value) )+
  labs(title="length-normalized RPK",y="RPK [reads per kb]")+
  geom_boxplot()+
  scale_y_log10(limits=c(0.5,50))+
  annotation_logticks(sides = "l")  +
  theme(
    axis.text.x = element_text(angle = 90),
  )

#calculate TMM offsets
y <- calcNormFactors(y)
y$samples$norm.factors
## [1] 0.0577 6.0829 1.2202 1.1648 1.1966 1.0466 1.1505 1.2543 1.1090

#extract GeTMM normalized counts
lcpm <- cpm(y, log=TRUE)

#create GeTMM plot
pnorm=ggplot(melt(lcpm),aes(x=Var2, y=value) )+
  labs(title="GeTMM-normalized",y="GeTMM [Gene length \ncorrected trimmed mean of M-values]")+
  geom_boxplot()+
  scale_y_log10(limits=c(1,30))+
  annotation_logticks(sides = "l")  +
  theme(
    axis.text.x = element_text(angle = 90),
  )

#arrange and output plot
g=grid.arrange(p,pnorm, nrow=1, ncol = 2)
ggsave(file=paste(FOLDER,  "/normaliztion.png", sep=""), g,height=7, width=10)


#calculate sample distance of normalized GeTMM counts for clustering heatmap
sampleDists = dist(t(lcpm))
sampleDistMatrix = as.matrix(sampleDists)

colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#plot heatmap using pheatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         filename=paste(FOLDER,  "/sample-clustering.png", sep=""))

#export count table 
cpm <- cpm(y, log=FALSE)

writeTABLE = as.data.frame(cpm) %>% rownames_to_column('gene_id')
write.table(writeTABLE, paste(FOLDER,  "/GeTMM_normalized.txt", sep=""), sep = "\t",  col.names=TRUE,  row.names=FALSE, quote=FALSE)

#stop if no DGE analysis is performed
if(DGE == "N"){
  q()
}

#---------------------------------------------------------------------------------------------------------

GENOlist=unlist(strsplit(GENO, "~"))
FOLDER=paste(rawFOLDER,"DGE/", sep="")


#perform DEseq2 analysis for all samples
#prepare in data
TXI <- tximport(FILES, type="salmon", tx2gene=TX2GENE)


DDSTXI <- DESeqDataSetFromTximport(TXI,
                                   colData = SAMPLES,
                                   design = ~ condition)


keep <- rowSums(counts(DDSTXI)) >= 10
DDSTXI <- DDSTXI[keep,]


DDSTXI$condition <- relevel(DDSTXI$condition, ref = refGENO)

DDSTXI <- DESeq(DDSTXI)

htmlwidgets::saveWidget(glimmaMDS(DDSTXI), paste(FOLDER,  "/sample_clustering.html", sep=""))

nROWvst = nrow(DDSTXI)
if(nROWvst>1000){
  vsd <- vst(DDSTXI, blind=FALSE)
}else{
  vsd <- varianceStabilizingTransformation(DDSTXI, blind=FALSE)
}

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
  pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         filename=paste(FOLDER,  "/sample-clustering.png", sep=""))



#do analysis for each genotype
for (refGENO in GENOlist) {
  GENOlist=GENOlist[GENOlist != refGENO]
  
  for (currGENO in GENOlist) {
    X=c(refGENO,currGENO)
    print(X)

    #subset input to current samples
    subSAMPLES = subset(SAMPLES, condition %in% c(refGENO, currGENO))
    subFILES = subset(FILES, names(FILES) %in% subSAMPLES$sample )

    #read in data
    TXI <- tximport(subFILES, type="salmon", tx2gene=TX2GENE)


    #construct DEseq data
    DDSTXI <- DESeqDataSetFromTximport(TXI,
                                   colData = subSAMPLES,
                                   design = ~ condition)


    #pre-filter input
    keep <- rowSums(counts(DDSTXI)) >= 10
    DDSTXI <- DDSTXI[keep,]

    DDSTXI$group <- factor(paste0(DDSTXI$condition))

    #define correct genotype as reference
    DDSTXI$condition <- relevel(DDSTXI$group, ref = refGENO)

    #perform analysis
    DDSTXI <- DESeq(DDSTXI)


    #create results table 
    RESname=resultsNames(DDSTXI)[2]
    RES <- results(DDSTXI, name=RESname,
      independentFiltering=TRUE, alpha=0.05, pAdjustMethod="BH", parallel=TRUE)

    #log fold change shrinkage
    RESlfc <- lfcShrink(DDSTXI, coef=2, type="apeglm")
    
    #sort for highest pval
    RESordered <- RES[order(RESlfc$padj),]
    
    #export datatable
    writeTABLE = as.data.frame(RESlfc) %>% rownames_to_column('gene_id')
    write.table(writeTABLE, 
          file=paste(FOLDER, refGENO, "_vs_", currGENO, "/dataTable.txt", sep="") , 
          sep = "\t",  col.names=TRUE,  row.names=FALSE,
          quote=FALSE)

    png(paste(FOLDER, refGENO, "_vs_", currGENO, "/data-plot.png", sep=""))
      DESeq2::plotMA(RESlfc, ylim=c(-8,8))
    dev.off()

    htmlwidgets::saveWidget(glimmaMA(DDSTXI), paste(FOLDER, refGENO, "_vs_", currGENO, "/volcano_counts_vs_foldChange.html", sep=""))  
    htmlwidgets::saveWidget(glimmaVolcano.DESeqDataSet.pValFixed(DDSTXI), paste(FOLDER, refGENO, "_vs_", currGENO, "/volcano_foldChange_vs_pValAdj.html", sep=""))  
    htmlwidgets::saveWidget(glimmaMDS(DDSTXI), paste(FOLDER, refGENO, "_vs_", currGENO, "/sample_clustering.html", sep=""))  

  }
}



