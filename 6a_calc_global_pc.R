# Global Methylation PCA

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

## Remove gene-associated CpGs:
ANNOT <- load_450k_annot()
ANNOT <- subset(ANNOT, ! is.na(UCSC_RefGene_Name))

## Load cohort datasets excluding CpGs analyzed:
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
tcga450k <- gbm450[rownames(gbm450) %in% ANNOT$Name, ]
rm(gbm450, gbmGENE, patients)

load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
dkfz450k <- gbm450[rownames(gbm450) %in% ANNOT$Name, ]
rm(gbm450, gbmGENE, patients)

# --------------------- Select Most Variable Coding & Intersect ---------------------
PROP <- 0.2

tcga450k <- select_most_var(tcga450k, PROP)
dkfz450k <- select_most_var(dkfz450k, PROP)

sele_sites <- intersect(rownames(tcga450k), rownames(dkfz450k))
length(sele_sites) #47225

## Export for LIMMA & input for KEGG (separate analyses):
# write.table(sele_sites, file=paste0(OUT_DIR, "global_dnam_CpGs.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

# --------------------- PCA --------------------- 
calc_pca <- function(dnam) {
  #'@describeIn Chen 2024 IJCDSE
  pca <- prcomp(t(dnam), center=TRUE, scale.=TRUE)
  print( summary(pca)$importance[c(2,3), 1:3] )
  
  resPca <- pca$x[ , c(1,2)]
  resPca <- apply(resPca, 2, StandardScale)
  return(as.data.frame(resPca))
}

## Subset for PCA here:
tcga450k <- tcga450k[sele_sites, ]
dkfz450k <- dkfz450k[sele_sites, ]

## Execute wrapper:
pca_tcga <- calc_pca(tcga450k)
pca_tcga$Cohort <- "TCGA"

pca_dkfz <- calc_pca(dkfz450k)
pca_dkfz$Cohort <- "DKFZ"

# write.csv(rbind(pca_tcga, pca_dkfz), paste0(OUT_DIR, "global_pca_by_cohort.csv"), row.names=TRUE, quote=FALSE)
