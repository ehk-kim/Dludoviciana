setwd("C:\\Users\\fanta\\Desktop\\DGE\\DESeq")

library(dplyr)
library(ggplot2)

# Functions
# NOTE: need to change the file name depending on the species

getUniques <- function(countsFile) {
  sp_spor <- countsFile %>%
    filter(gam == 0 & spor != 0)
  
  sp_gam <- countsFile %>%
    filter(gam != 0, spor == 0)
  
  write(row.names(sp_spor), file = "ptab_sporo_unique.txt", sep = "\n")
  write(row.names(sp_gam), file = "ptab_gameto_unique.txt", sep = "\n")
}

#####

# Arabidopsis thaliana
# Sporophyte only: 4782 (18.41%)
# Gametophyte only: 1025 (3.95%)

athaCounts <- read.delim("atha_counts.txt", row.name="gene_id")
athaCounts$gam <- rowSums(athaCounts[,1:2])
athaCounts$spor <- rowSums(athaCounts[,3:5])

getUniques(athaCounts)

# Anthoceros agrestis
# Sporophyte only: 2923 (12.50%)
# Gametophyte only: 326 (1.39%)

aagrCounts <- read.delim("aagr_counts.txt", row.name="gene_id")
aagrCounts$gam <- rowSums(aagrCounts[,1:2])
aagrCounts$spor <- rowSums(aagrCounts[,3:4])

getUniques(aagrCounts)

# Dryopteris ludoviciana
# Sporophyte only: 22953 (12.33%)
# Gametophyte only: 119819 (64.37%)

dludCounts <- read.delim("dlud_counts_genelvl.txt", row.name="gene_id")
dludCounts$gam <- rowSums(dludCounts[,1:3])
dludCounts$spor <- rowSums(dludCounts[,4:6])

getUniques(dludCounts)

# Polypodium amorphum
# Sporophyte only: 1604 (3.65%)
# Gametophyte only: 1255 (2.86%)

pamoCounts <- read.delim("pamo_counts_genelvl.txt", row.name="gene_id")
pamoCounts$gam <- rowSums(pamoCounts[,1:3])
pamoCounts$spor <- rowSums(pamoCounts[,4:6])

getUniques(pamoCounts)

# Physcomitrium patens
# Sporophyte only: 1273 (4.47%)
# Gametophyte only: 1330 (4.67%)

ppatCounts <- read.delim("ppat_counts.txt", row.name="gene_id")
ppatCounts$gam <- rowSums(ppatCounts[,1:3])
ppatCounts$spor <- rowSums(ppatCounts[,4:6])

getUniques(ppatCounts)

# Pinus tabuliformis
# Sporophyte only: 2991 (3.37%)
# Gametophyte only: 289 (0.33%)

ptabCounts <- read.delim("ptab_counts.txt", row.name="gene_id")
ptabCounts$gam <- rowSums(ptabCounts[,1:4])
ptabCounts$spor <- rowSums(ptabCounts[,5:7])

getUniques(ptabCounts)
