setwd("C:\\Users\\fanta\\Desktop\\DGE\\TopGO")

BiocManager::install("topGO")
library(topGO)

#####

#-----Function definitions-----#

topGO_object <- function(filePath, geneUniverse, geneID2GO, GOdescription) {
  sp <- read.delim(filePath, header=FALSE, sep='\t')
  sp$V1 <- gsub(">","",as.character(sp$V1))
  sp <- as.character(sp$V1)
  
  # Set gene universe
  geneList <- factor(as.integer(geneUniverse %in% sp))
  names(geneList) <- geneUniverse
  
  # Create topGO object
  sp_GO <- new("topGOdata", description=GOdescription, ontology="BP", allGenes=geneList, annot=annFUN.gene2GO, gene2GO=geneID2GO)
  
  return(sp_GO)
}

GO_results <- function(sp_GO, outFile) {
  resFisher <- runTest(sp_GO, algorithm="classic", statistic="fisher")
  allGO <- usedGO(sp_GO)
  all_res <- GenTable(sp_GO, classicFisher=resFisher, orderBy="resultFisher", ranksOf="classicFisher", topNodes=length(allGO))
  
  p.adj <- round(p.adjust(all_res$classicFisher,method="BH"),digits = 4)
  res <- cbind(all_res,p.adj)
  res <- res[order(res$p.adj),]
  
  write.table(res[1:50,], outFile, sep=",", quote=FALSE, row.names=FALSE)
}

#####

#-----Analysis-----#

# Physcomitrella patens

# Creating annotation file for P. patens
ppat_annot_raw <-  read.delim("Ppat_V3.3_annotation.txt", header=TRUE, sep='\t')
ppat_annot <- subset(ppat_annot_raw, select=c(locusName, GO))
ppat_annot <- ppat_annot[!(ppat_annot$GO==""),]
write.table(ppat_annot, file="ppat_annot.txt", row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE)

# Importing files for topGO
geneID2GO_ppat <- readMappings(file="ppat_annot.txt")
geneUniverse_ppat <- names(geneID2GO_ppat)

# Create TopGO object
ppat_spor_GO <- topGO_object(filePath = "C:\\Users\\fanta\\Desktop\\DGE\\gogetter\\Ppat_sporo_upreg.txt",
                             geneUniverse = geneUniverse_ppat,
                             geneID2GO = geneID2GO_ppat,
                             GOdescription = "Physcomitrium patens upregulated genes in sporophyte")

ppat_gam_GO <- topGO_object(filePath = "C:\\Users\\fanta\\Desktop\\DGE\\gogetter\\Ppat_gameto_upreg.txt",
                             geneUniverse = geneUniverse_ppat,
                             geneID2GO = geneID2GO_ppat,
                             GOdescription = "Physcomitrium patens upregulated genes in gametophyte")

# Complete statistical analysis
GO_results(sp_GO = ppat_spor_GO, outFile = "Ppat_sporo_topGO.csv")
GO_results(sp_GO = ppat_gam_GO, outFile = "Ppat_gam_topGO.csv")

#-------------------------------------------------------------------------------

# Anthoceros agrestis

# Importing files for topGO
geneID2GO_aagr <- readMappings(file="aagr_annot.txt")
geneUniverse_aagr <- names(geneID2GO_aagr)

# Create TopGO object
aagr_spor_GO <- topGO_object(filePath = "C:\\Users\\fanta\\Desktop\\DGE\\gogetter\\Aagr_sporo_upreg.txt",
                             geneUniverse = geneUniverse_aagr,
                             geneID2GO = geneID2GO_aagr,
                             GOdescription = "Anthoceros agrestis upregulated genes in sporophyte")

aagr_gam_GO <- topGO_object(filePath = "C:\\Users\\fanta\\Desktop\\DGE\\gogetter\\Aagr_gameto_upreg.txt",
                            geneUniverse = geneUniverse_aagr,
                            geneID2GO = geneID2GO_aagr,
                            GOdescription = "Anthoceros agrestis upregulated genes in gametophyte")

# Complete statistical analysis
GO_results(sp_GO = aagr_spor_GO, outFile = "Aagr_sporo_topGO.csv")
GO_results(sp_GO = aagr_gam_GO, outFile = "Aagr_gam_topGO.csv")

#-------------------------------------------------------------------------------

# Polypodium amorphum

# Importing files for topGO
geneID2GO_pamo <- readMappings(file="pamo_annot.txt")
geneUniverse_pamo <- names(geneID2GO_pamo)

# Create TopGO object
pamo_spor_GO <- topGO_object(filePath = "Pamo_sporo_upreg.txt",
                             geneUniverse = geneUniverse_pamo,
                             geneID2GO = geneID2GO_pamo,
                             GOdescription = "Polypodium amorphum upregulated genes in sporophyte")

pamo_gam_GO <- topGO_object(filePath = "Pamo_gameto_upreg.txt",
                            geneUniverse = geneUniverse_pamo,
                            geneID2GO = geneID2GO_pamo,
                            GOdescription = "Polypodium amorphum upregulated genes in gametophyte")

# Complete statistical analysis
GO_results(sp_GO = pamo_spor_GO, outFile = "Pamo_sporo_topGO.csv")
GO_results(sp_GO = pamo_gam_GO, outFile = "Pamo_gam_topGO.csv")
#-------------------------------------------------------------------------------

# Pinus tabuliformis

# Importing files for topGO
geneID2GO_ptab <- readMappings(file="Ptab_annotation.txt")
geneUniverse_ptab <- names(geneID2GO_ptab)

# Create TopGO object
ptab_spor_GO <- topGO_object(filePath = "C:\\Users\\fanta\\Desktop\\DGE\\gogetter\\ptab_sporo_upreg.txt",
                             geneUniverse = geneUniverse_ptab,
                             geneID2GO = geneID2GO_ptab,
                             GOdescription = "Anthoceros agrestis upregulated genes in sporophyte")

ptab_gam_GO <- topGO_object(filePath = "C:\\Users\\fanta\\Desktop\\DGE\\gogetter\\ptab_gameto_upreg.txt",
                            geneUniverse = geneUniverse_ptab,
                            geneID2GO = geneID2GO_ptab,
                            GOdescription = "Anthoceros agrestis upregulated genes in gametophyte")

# Complete statistical analysis
GO_results(sp_GO = ptab_spor_GO, outFile = "ptab_sporo_topGO.csv")
GO_results(sp_GO = ptab_gam_GO, outFile = "ptab_gam_topGO.csv")

#-------------------------------------------------------------------------------

# Arabidopsis thaliana

# Creating annotation file for A. thaliana
atha_annot_raw <-  read.delim("Atha_TAIR10_annotation.txt", header=FALSE, sep='\t')
atha_annot <- subset(atha_annot_raw, select=c(V2, V10))
atha_annot <- atha_annot[!(atha_annot$V10==""),]
write.table(atha_annot, file="atha_annot.txt", row.names=FALSE, quote=FALSE, sep='\t', col.names=FALSE)

# Importing files for topGO
geneID2GO_atha <- readMappings(file="atha_annot.txt")
geneUniverse_atha <- names(geneID2GO_atha)

# Create TopGO object
atha_spor_GO <- topGO_object(filePath = "C:\\Users\\fanta\\Desktop\\DGE\\gogetter\\Atha_sporo_upreg.txt",
                             geneUniverse = geneUniverse_atha,
                             geneID2GO = geneID2GO_atha,
                             GOdescription = "Anthoceros agrestis upregulated genes in sporophyte")

atha_gam_GO <- topGO_object(filePath = "C:\\Users\\fanta\\Desktop\\DGE\\gogetter\\Atha_gameto_upreg.txt",
                            geneUniverse = geneUniverse_atha,
                            geneID2GO = geneID2GO_atha,
                            GOdescription = "Anthoceros agrestis upregulated genes in gametophyte")

# Complete statistical analysis
GO_results(sp_GO = atha_spor_GO, outFile = "atha_sporo_topGO.csv")
GO_results(sp_GO = atha_gam_GO, outFile = "atha_gam_topGO.csv")
