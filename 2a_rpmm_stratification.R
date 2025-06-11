# Recursively Partitioned Mixture Model (RPMM) Clustering of Locus DNAm Data

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(RPMM)
library(seriation)

get_rpmm_samp_order <- function(rpmmClusters, Y_inv) {
  #'@description Computes sample orders for heatmap visualization
  #'@describeIn  Chen et al. Int J Canc 2024
  #'@param rpmmClusters data.frame w/ row.names=SampleID & 1 column named "RPMM"
  #'@param Y_inv Input matrix for RPMM computation
  
  # Extract sample order 1 RPMM cluster at a time:
  ordered_samps <- c()
  for(cl in names(table(rpmmClusters$RPMM)) ) {
    samps <- rownames(rpmmClusters)[rpmmClusters$RPMM == cl]
    msub <- Y_inv[samps, ]
    s_i <- seriate(msub, margin=2)
    so_i <- get_order(s_i)
    ordered_samps <- c(ordered_samps, names(so_i))
  }
  
  # Organize results:
  res <- data.frame(
    Accession = ordered_samps,
    RPMMSampleOrder = 1:length(ordered_samps)
  )
  return(res)
}

run_custom_rpmm <- function(dnam, max_level=2) {
  #'@description Wrapper to run RPMM & extract cluster info
  #'@describeIn Chen et al. Int J Canc 202
  
  # Execute standard RPMM:
  Y_inv <- t(dnam)
  res <- blcTree(Y_inv, verbose=1, maxlevel=max_level)
  print(res)
  plot(res)
  
  cls <- blcTreeLeafClasses(res)
  res_rpmm <- data.frame(row.names=rownames(Y_inv), RPMM=cls)
  
  # Extract sample order:
  samp_order <- get_rpmm_samp_order(res_rpmm, Y_inv)
  
  res_rpmm <- merge(samp_order, res_rpmm, by.x="Accession", by.y="row.names")
  
  # Collapse to maximum cluster leve of choice:
  res_rpmm$Cluster <- substr(res_rpmm$RPMM, 2, max_level)
  return(res_rpmm)
}


## Load cohort datasets:
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
tcgaGENE <- gbmGENE
rm(gbm450, gbmGENE, patients)

load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
dkfzGENE <- gbmGENE
rm(gbm450, gbmGENE, patients)

## Execute wrapper:
tcga_rpmm <- run_custom_rpmm(tcgaGENE, 2)
tcga_rpmm$Cohort <- "TCGA"

dkfz_rpmm <- run_custom_rpmm(dkfzGENE, 2)
dkfz_rpmm$Cohort <- "DKFZ"

# write.csv(rbind(tcga_rpmm, dkfz_rpmm), paste0(OUT_DIR, "rpmm_by_cohort.csv"), row.names=FALSE, quote=FALSE)
