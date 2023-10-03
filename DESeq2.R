setwd("C:\\Users\\fanta\\Desktop\\DGE\\DESeq")

library(DESeq2)
library(pheatmap)
library(ggplot2)

#####

# Define functions

createDESeq <- function(countsFile, colsFile) {
  counts <- read.delim(countsFile, row.name="gene_id")
  cols <- read.delim(colsFile, row.names=1)
  
  ddsObj <- DESeqDataSetFromMatrix(countData=counts, colData=cols, design=~condition)
  ddsObj$condition <- relevel(ddsObj$condition, ref="gametophyte")
  ddsObj <- DESeq(ddsObj)
  
  return(ddsObj)
}

resLFC <- function(ddsObj) {
  res_GA <- results(ddsObj, lfcThreshold=2, altHypothesis="greaterAbs", alpha = 0.05)
  res_LFC <- lfcShrink(ddsObj, coef="condition_sporophyte_vs_gametophyte", type="apeglm", res=res_GA)
  summary(res_LFC)
  
  return(res_LFC)
}

sporoUpreg <- function(ddsRes) {
  vec <- c()
  for (i in  1:length(ddsRes$log2FoldChange)) {
    if (ddsRes$log2FoldChange[i] > 0 & is.na(ddsRes$log2FoldChange[i])==FALSE & is.na(ddsRes$padj[i])==FALSE & ddsRes$padj[i] < 0.05) {
      vec <- append(vec, names(ddsRes$log2FoldChange[i]))
    }
  }
  return(vec)
}

gametoUpreg <- function(ddsRes) {
  vec <- c()
  for (i in  1:length(ddsRes$log2FoldChange)) {
    if (ddsRes$log2FoldChange[i] < 0 & is.na(ddsRes$log2FoldChange[i])==FALSE & is.na(ddsRes$padj[i])==FALSE & ddsRes$padj[i] < 0.05) {
      vec <- append(vec, names(ddsRes$log2FoldChange[i]))
    }
  }
  return(vec)
}

makeHeatmap <- function(ddsObj) {
  ntd <- normTransform(ddsObj)
  select <- order(rowMeans(counts(ddsObj)), decreasing = TRUE) [1:50]
  df <- as.data.frame(colData(ddsObj) [,c("species", "condition")])
  phm <- assay(ntd)
  pheatmap(phm[select,], cutree_cols = 2, annotation_col=df, cluster_rows = TRUE)
}

makeHeatmap2 <- function(ddsObj) {
  vsd<-varianceStabilizingTransformation(ddsObj)
  sampleDists <- dist(t(assay(vsd)))
  sampleDistMatrix <- as.matrix(sampleDists)
  rownames(sampleDistMatrix) <- paste(vsd$treatment)
  colnames(sampleDistMatrix) <- paste(vsd$condition)
  ph <- pheatmap(sampleDistMatrix,
           clustering_distance_rows=sampleDists, 
           clustering_distance_cols=sampleDists)
  #ph <-  ph + scale_fill_continuous(limits = c(0,700), type="viridis")
  ph
}

#####

# Create DESeq objects

physco_dds <- createDESeq("physco_counts.txt", "physco_sample_info.txt")
antho_dds <- createDESeq("antho_counts.txt", "antho_sample_info.txt")
arab_dds <- createDESeq("arab_counts.txt", "arab_sample_info.txt")
polypod_dds <- createDESeq("polypod_counts.txt", "polypod_sample_info.txt")
pinus_dds <- createDESeq("pinus_counts.txt", "pinus_sample_info.txt")

#------------LFC Shrinkage------------#

#LFC shrinkage; reduces variability of lowly expressed genes
#LFC > 0 (Up): num genes upregulated in sporo compared to gameto
#LF < 0 (Down): num genes downregulated in sporo compared to gameto

physco_LFC <- resLFC(physco_dds)

# out of 25719 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1790, 7%
# LFC < 0 (down)     : 1395, 5.4%
# outliers [1]       : 49, 0.19%
# low counts [2]     : 2952, 11%

physco_sporo_up <- sporoUpreg(physco_LFC)
physco_gameto_up <- gametoUpreg(physco_LFC)
# got list of genes; not the same format as fasta headers. same as gene name in .gff3 files

antho_LFC <- resLFC(antho_dds)

# USING GENES
# out of 20986 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 742, 3.5%
# LFC < 0 (down)     : 292, 1.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 2009, 9.6%

#USING CDS
# out of 23377 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 515, 2.2%
# LFC < 0 (down)     : 212, 0.91%
# outliers [1]       : 0, 0%
# low counts [2]     : 6544, 28%

antho_sporo_up <- sporoUpreg(antho_LFC)
antho_gameto_up <- gametoUpreg(antho_LFC)

arab_LFC <- resLFC(arab_dds)

# out of 23868 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4309, 18%
# LFC < 0 (down)     : 1743, 7.3%
# outliers [1]       : 3, 0.013%
# low counts [2]     : 339, 1.4%

arab_sporo_up <- sporoUpreg(arab_LFC)
arab_gameto_up <- gametoUpreg(arab_LFC)

polypod_LFC <- resLFC(polypod_dds)

# of 75009 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 452, 0.6%
# LFC < 0 (down)     : 857, 1.1%
# outliers [1]       : 1622, 2.2%
# low counts [2]     : 40430, 54%

polypod_sporo_up <- sporoUpreg(polypod_LFC)
polypod_gameto_up <- gametoUpreg(polypod_LFC)

pinus_LFC <- resLFC(pinus_dds)
# out of 59524 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 2321, 3.9%
# LFC < 0 (down)     : 498, 0.84%
# outliers [1]       : 2172, 3.6%
# low counts [2]     : 33729, 57%

pinus_sporo_up <- sporoUpreg(pinus_LFC)
pinus_gameto_up <- gametoUpreg(pinus_LFC)

#------------Create graphs------------#

drawLines <- function() abline(h=c(-2,2),col="dodgerblue",lwd=2)

#shows LFC > 2
#par(mfrow=c(2,2))
plotMA(physco_LFC, ylim=c(-15,15), main="L2FC vs base mean: Physcomitrium patens (moss)"); drawLines()
plotMA(antho_LFC, ylim=c(-15,15), main="L2FC vs base mean: Anthoceros agrestis (hornwort)"); drawLines()
plotMA(arab_LFC, ylim=c(-15,15), main="L2FC vs base mean: Arabidopsis thaliana (angiosperm)"); drawLines()
plotMA(polypod_LFC, ylim=c(-15,15), main="L2FC vs base mean: Polypodium amorphum (fern)"); drawLines()
plotMA(pinus_LFC, ylim=c(-15,15), main="L2FC vs base mean: Pinus_tabuliformis (gymnosperm)"); drawLines()

#####

# Gene heatmap

# Shows heatmap of different genes
makeHeatmap(antho_dds)
makeHeatmap(physco_dds)
makeHeatmap(arab_dds)
makeHeatmap(polypod_dds)
makeHeatmap(pinus_dds)

# Compares overall gene profiles
makeHeatmap2(antho_dds)
makeHeatmap2(physco_dds)
makeHeatmap2(arab_dds)
makeHeatmap2(polypod_dds)
png(file="C:\\Users\\fanta\\Desktop\\DGE\\pinus_heatmap2.png",
    width=800, height=800)
makeHeatmap2(pinus_dds)
dev.off()

# Principal Component Analysis
# plotPCA(vsd, intgroup=c("condition"))

#####