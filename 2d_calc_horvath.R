# Calculate Horvath Methylation Age

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

library(wateRmelon)

calc_horvath <- function(betas) {
  #'@description Wrapper for wateRmelon implementation of Horvath Age
  res <- wateRmelon::agep(betas)
  res$Accession <- rownames(res)
  rownames(res) <- NULL
  return(res[ , c(3,1)])
}

## Load DNAm matrices for GBM:
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
tcga450k <- gbm450
rm(gbm450, gbmGENE, patients)

load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
dkfz450k <- gbm450
rm(gbm450, gbmGENE, patients)

## Execute wrapper (ref. Chen et al. 2019)
horvath_tcga <- calc_horvath(tcga450k)
horvath_dkfz <- calc_horvath(dkfz450k)

write.csv(rbind(horvath_tcga, horvath_dkfz), file=paste0(OUT_DIR,"horvath_by_cohort.csv"), row.names=FALSE, quote=FALSE)

## Verify NO overlap btwn Signature & target loci:
