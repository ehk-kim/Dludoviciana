#setwd("C:\\Users\\fanta\\Desktop\\DGE\\DESeq")

library(DESeq2)
library(pheatmap)
library(ggplot2)
library(dplyr)

#####

# Define functions

createDESeq <- function(countsFile, colsFile) {
  counts <- read.delim(countsFile, row.name="gene_id")
  cols <- read.delim(colsFile, row.names=1)
  
  ddsObj <- DESeqDataSetFromMatrix(countData=counts, colData=cols, design=~condition)
  # ddsObj$condition <- relevel(ddsObj$condition, ref="gametophyte")
  ddsObj <- DESeq(ddsObj)
  
  res <- results(ddsObj, alpha = 0.05, lfcThreshold = 2, contrast = c("condition", "gametophyte", "sporophyte"))
  
  # smallestGroupSize <- 2
  # keep <- rowSums(counts(ddsObj) >= 10) >= smallestGroupSize
  # ddsObj <- ddsObj[keep,]
  res
  
  return(ddsObj)
}

resLFC <- function(ddsObj) {
  #ddsObj$condition <- factor(ddsObj$condition, levels = c("gametophyte","sporophyte"))
  #ddsObj$condition <- relevel(ddsObj$condition, ref = "sporophyte")
  res_GA <- results(ddsObj, lfcThreshold=2, altHypothesis="greaterAbs", alpha = 0.05, contrast = c("condition", "gametophyte", "sporophyte"))
  #res_GA <- results(ddsObj, lfcThreshold=2, altHypothesis="greaterAbs", alpha = 0.05)
  res_LFC <- lfcShrink(ddsObj, type = "ashr", res=res_GA, contrast = c("condition", "gametophyte", "sporophyte"))
  #res_LFC <- lfcShrink(ddsObj, type="apeglm", res=res_GA)
  summary(res_LFC)
  
  return(res_LFC)
}

gametoUpreg <- function(ddsRes) {
  ddsRes <-subset(ppat_LFC, padj < 0.05)
  ddsRes1 <- as.data.frame(ddsRes)
  vec <- ddsRes1 %>% filter(log2FoldChange >= 2)
  vec1 <- rownames(vec)
  
  # vec <- c()
  # for (i in  1:length(ddsRes$log2FoldChange)) {
  #   if (ddsRes$log2FoldChange[i] >= 2 & is.na(ddsRes$log2FoldChange[i])==FALSE & is.na(ddsRes$padj[i])==FALSE & ddsRes$padj[i] < 0.05) {
  #     vec <- append(vec, names(ddsRes$log2FoldChange[i]))
  #   }
  # }
  return(vec1)
}

gametoDownreg <- function(ddsRes) {
  ddsRes <-subset(ppat_LFC, padj < 0.05)
  ddsRes1 <- as.data.frame(ddsRes)
  vec <- ddsRes1 %>% filter(log2FoldChange <= -2)
  vec1 <- rownames(vec)
  # vec <- c()
  # for (i in  1:length(ddsRes$log2FoldChange)) {
  #   if (ddsRes$log2FoldChange[i] <= -2 & is.na(ddsRes$log2FoldChange[i])==FALSE & is.na(ddsRes$padj[i])==FALSE & ddsRes$padj[i] < 0.05) {
  #     vec <- append(vec, names(ddsRes$log2FoldChange[i]))
  #   }
  # }
  return(vec1)
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

pamo_counts <- read.delim("pamo_counts.txt")
pamo_test <- pamo_counts %>%
                group_by(gene_id) %>%
                summarize(gameto_2 = sum(gameto_2),
                          gameto_4 = sum(gameto_4),
                          gameto_6 = sum(gameto_6),
                          sporo_1 = sum(sporo_1),
                          sporo_3 = sum(sporo_3),
                          sporo_5 = sum(sporo_5))
pamo_test <- as.data.frame(pamo_test)
write.table(pamo_test, file = "pamo_counts_genelvl.txt", quote = F, row.names = F, sep = "\t")
test <- read.delim("pamo_counts_genelvl.txt")

dlud_counts <- read.delim("dlud_counts.txt")
dlud_genelvl <- dlud_counts %>%
  group_by(gene_id) %>%
  summarize(gameto_206 = sum(gameto_206),
            gameto_199 = sum(gameto_199),
            gameto_204 = sum(gameto_204),
            sporo_200 = sum(sporo_200),
            sporo_204 = sum(sporo_204),
            sporo_206 = sum(sporo_206))
dlud_genelvl <- as.data.frame(dlud_genelvl)
write.table(dlud_genelvl, file = "dlud_counts_genelvl.txt", quote = F, row.names = F, sep = "\t")

ppat_dds <- createDESeq("ppat_counts.txt", "ppat_sample_info.txt")
aagr_dds <- createDESeq("aagr_counts.txt", "aagr_sample_info.txt")
atha_dds <- createDESeq("atha_counts.txt", "atha_sample_info.txt")
pamo_dds <- createDESeq("pamo_counts_genelvl.txt", "pamo_sample_info.txt")
ptab_dds <- createDESeq("ptab_counts.txt", "ptab_sample_info.txt")
dlud_dds <- createDESeq("dlud_counts_genelvl.txt", "dlud_sample_info.txt")

#------------LFC Shrinkage------------#

#LFC shrinkage; reduces variability of lowly expressed genes
#LFC > 0 (Up): num genes upregulated in sporo compared to gameto
#LF < 0 (Down): num genes downregulated in sporo compared to gameto
write(aagr_sporo_up, file = "aagr_sporo_upreg.txt", sep = "\n")

ppat_LFC <- resLFC(ppat_dds)

# out of 28457 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 1932, 6.8%
# LFC < 0 (down)     : 1567, 5.5%
# outliers [1]       : 60, 0.21%
# low counts [2]     : 3283, 12%

ppat_gameto_up <- gametoUpreg(ppat_LFC)
ppat_gameto_down <- gametoDownreg(ppat_LFC)
# got list of genes; not the same format as fasta headers. same as in .gff3 files

aagr_LFC <- resLFC(aagr_dds)

# out of 23379 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 767, 3.3%
# LFC < 0 (down)     : 311, 1.3%
# outliers [1]       : 0, 0%
# low counts [2]     : 2253, 9.6%

aagr_sporo_up <- sporoUpreg(aagr_LFC)
aagr_gameto_up <- gametoUpreg(aagr_LFC)

atha_LFC <- resLFC(atha_dds)

# out of 25982 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 4281, 16%
# LFC < 0 (down)     : 1971, 7.6%
# outliers [1]       : 3, 0.012%
# low counts [2]     : 2510, 9.7%

atha_sporo_up <- sporoUpreg(atha_LFC)
atha_gameto_up <- gametoUpreg(atha_LFC)

pamo_LFC <- resLFC(pamo_dds)

# out of 43913 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 536, 1.2%
# LFC < 0 (down)     : 1107, 2.5%
# outliers [1]       : 600, 1.4%
# low counts [2]     : 11068, 25%

pamo_sporo_up <- sporoUpreg(pamo_LFC)
pamo_gameto_up <- gametoUpreg(pamo_LFC)

ptab_LFC <- resLFC(ptab_dds)

# out of 88670 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 7022, 7.9%
# LFC < 0 (down)     : 1873, 2.1%
# outliers [1]       : 3645, 4.1%
# low counts [2]     : 6877, 7.8%

ptab_sporo_up <- sporoUpreg(ptab_LFC)
ptab_gameto_up <- gametoUpreg(ptab_LFC)

dlud_LFC <- resLFC(dlud_dds)

# out of 186128 with nonzero total read count
# adjusted p-value < 0.05
# LFC > 0 (up)       : 571, 0.31%
# LFC < 0 (down)     : 2035, 1.1%
# outliers [1]       : 26826, 14%
# low counts [2]     : 71531, 38%

dlud_sporo_up <- sporoUpreg(dlud_LFC)
dlud_gameto_up <- gametoUpreg(dlud_LFC)

# NOTES---------
dlud_vsd <- vst(dlud_dds)
sigs <- subset(dlud_LFC, padj < 0.05)
sigs2 <- rownames(sigs)
vsd_subset <- dlud_vsd[rownames(dlud_vsd) %in% sigs2, ]
summary(vsd_subset)
obj <- plotPCA(vsd_subset, intgroup=c("species", "condition"))

q <- ggplot_build(obj)
df <- as.data.frame(q$plot$data)
df2 <- df %>%
  dplyr::select(PC1, PC2, species, condition, species)

ggplot(data = df2, mapping = aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3, alpha = 0.75) + theme_classic()

#------------Create graphs------------#

drawLines <- function() abline(h=c(-2,2),col="dodgerblue",lwd=2)

#shows LFC > 2
#par(mfrow=c(2,2))
plotMA(ppat_LFC, ylim=c(-15,15), main="L2FC vs base mean: Physcomitrium patens (moss)"); drawLines()
plotMA(aagr_LFC, ylim=c(-15,15), main="L2FC vs base mean: Anthoceros agrestis (hornwort)"); drawLines()
plotMA(atha_LFC, ylim=c(-15,15), main="L2FC vs base mean: Arabidopsis thaliana (angiosperm)"); drawLines()
plotMA(pamo_LFC, ylim=c(-15,15), main="L2FC vs base mean: Polypodium amorphum (fern)"); drawLines()
plotMA(ptab_LFC, ylim=c(-15,15), main="L2FC vs base mean: Pinus_tabuliformis (gymnosperm)"); drawLines()
plotMA(dlud_LFC, ylim=c(-15,15), main="L2FC vs base mean: Dryopteris_ludoviciana (fern)"); drawLines()

#####

# Gene heatmap

# Shows heatmap of different genes
png(file="C:\\Users\\fanta\\Desktop\\DGE\\Ppat_heatmap1.png",
    width=800, height=800)
makeHeatmap(antho_dds)
dev.off()
makeHeatmap(physco_dds)
makeHeatmap(arab_dds)
makeHeatmap(polypod_dds)
makeHeatmap(pinus_dds)
makeHeatmap(dlud_dds)

# Compares overall gene profiles
makeHeatmap2(antho_dds)
makeHeatmap2(physco_dds)
makeHeatmap2(arab_dds)
makeHeatmap2(polypod_dds)
makeHeatmap2(pinus_dds)
makeHeatmap2(dlud_dds)

#####

#ppat_vsd <- vst(ppat_dds)
#ppat_df <- read.csv()

aagr_vsd <- vst(aagr_dds)
aagr_df <- read.csv("aagr_homeodomain.csv", header = FALSE)
aagr_df <- as.data.frame(aagr_df$V2, row.names = aagr_df$V1)

aagr_homeoboxes <- aagr_vsd[rownames(aagr_vsd) %in% rownames(aagr_df), ]
aagr_homeoboxes_ordered <- aagr_homeoboxes[rownames(aagr_df), ]

pheatmap(assay(aagr_homeoboxes_ordered),
         scale="row", annotation_row = aagr_df, 
         show_rownames = TRUE,cluster_rows = FALSE, treeheight_col = 0,
         color = colorRampPalette(c("gray18", "dodgerblue2", "gold"))(100))

atha_vsd <- vst(atha_dds)
atha_df <- read.csv("atha_homeodomain.csv", header = FALSE)
atha_df <- as.data.frame(atha_df$V2, row.names = atha_df$V1)

atha_homeoboxes <- atha_vsd[rownames(atha_vsd) %in% rownames(atha_df), ]
atha_homeoboxes_ordered <- atha_homeoboxes[rownames(atha_df), ]

pheatmap(assay(atha_homeoboxes_ordered),
         scale="row", annotation_row = atha_df, 
         show_rownames = TRUE,cluster_rows = FALSE, treeheight_col = 0,
         color = colorRampPalette(c("gray18", "dodgerblue2", "gold"))(100),
         annotation_colors = ann_colors)

pamo_vsd <- vst(ptab_dds)
pamo_df <- read.csv("pamo_homeodomain.csv", header = FALSE)
pamo_df <- as.data.frame(pamo_df$V2, row.names = pamo_df$V1)
colnames(pamo_df) = c('Homeodomains')

pamo_homeoboxes <- pamo_vsd[rownames(pamo_vsd) %in% rownames(pamo_df), ]
pamo_homeoboxes_ordered <- pamo_homeoboxes[rownames(pamo_df), ]

pheatmap(assay(pamo_homeoboxes_ordered),
         scale="row", annotation_row = pamo_df, 
         show_rownames = TRUE,cluster_rows = FALSE, treeheight_col = 0,
         color = colorRampPalette(c("gray18", "dodgerblue2", "gold"))(100),
         annotation_colors = ann_colors)

ptab_vsd <- vst(ptab_dds)
ptab_df <- read.csv("ptab_homeodomain.csv", header = FALSE)
ptab_df <- as.data.frame(ptab_df$V2, row.names = ptab_df$V1)

ptab_homeoboxes <- ptab_vsd[rownames(ptab_vsd) %in% rownames(ptab_df), ]
ptab_homeoboxes_ordered <- ptab_homeoboxes[rownames(ptab_df), ]

pheatmap(assay(ptab_homeoboxes_ordered),
         scale="row", annotation_row = ptab_df, 
         show_rownames = TRUE,cluster_rows = FALSE, treeheight_col = 0,
         color = colorRampPalette(c("gray18", "dodgerblue2", "gold"))(100))

dlud_vsd <- vst(dlud_dds)
dlud_df <- read.csv("Dlud_homeodomain_1.csv", header = FALSE)
dlud_df <- as.data.frame(dlud_df$V2, row.names = dlud_df$V1)

dlud_homeoboxes <- dlud_vsd[rownames(dlud_vsd) %in% rownames(dlud_df), ]
dlud_homeoboxes_ordered <- dlud_homeoboxes[rownames(dlud_df), ]

pheatmap(assay(dlud_homeoboxes_ordered),
         scale="row", annotation_row = dlud_df, 
         show_rownames = TRUE,cluster_rows = FALSE, treeheight_col = 0,
         color = colorRampPalette(c("gray18", "dodgerblue2", "gold"))(100))