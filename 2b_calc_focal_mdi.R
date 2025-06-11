# Calculate Focal Methylation Dysregulation

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

calc_mdi <- function(betas, ctrlBetas) {
  #'@description Wrapper to calculate Epigenetic Dysregulation (Methylation Dysregulation Index) proposed by Salas et al. Epigenet.
  #'@describedIn Chen et al. Int J Canc
  
  stopifnot( identical(rownames(betas), rownames(ctrlBetas)) ) #checkpoint
  
  absDeviats <- abs(betas - matrixStats::rowMedians(ctrlBetas))
  
  mdi <- data.frame(
    row.names = NULL,
    Accession = colnames(betas),
    MDI = matrixStats::colMeans2(absDeviats, na.rm=TRUE)
  )
  return(mdi)
}

## Normal Control data:
load(paste0(OUT_DIR, "healthy_brain_dnam.RData"))
ctrl450k <- pfcGENE
rm(pfc450, pfcGENE, patients)

## Load DNAm matrices for GBM:
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
tcga450k <- gbmGENE
rm(gbm450, gbmGENE, patients)

load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
dkfz450k <- gbmGENE
rm(gbm450, gbmGENE, patients)

## Execute wrapper:
mdi_tcga <- calc_mdi(tcga450k, ctrl450k)
mdi_tcga$Cohort <- "TCGA"

mdi_dkfz <- calc_mdi(dkfz450k, ctrl450k)
mdi_dkfz$Cohort <- "DKFZ"

write.csv(rbind(mdi_tcga, mdi_dkfz), paste0(OUT_DIR,"mdi_by_cohort.csv"), row.names=FALSE, quote=FALSE)
