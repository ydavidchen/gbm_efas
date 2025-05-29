# Recursively Partitioned Mixture Model (RPMM) Clustering of Locus DNAm Data

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(RPMM)

get_rpmm_samp_order <- function(rpmmClusters, Y_inv, sampleIdName) {
  #'@description Retrieves sample orders fo heat map visualization
  #'@describeIn  Chen et al. Int J Canc 2024
  #'@param rpmmClusters data.frame & 1 column named "RPMM"
  #'@param Y_inv Input matrix for RPMM computation
  
  sampOrder <- c()
  
  for(cl in names(table(rpmmClusters$RPMM)) ) {
    samps <- rownames(rpmmClusters)[rpmmClusters$RPMM == cl]
    clu <- Y_inv[rownames(Y_inv) %in% samps, ]
    s_i <- seriation::seriate(clu, margin=2)
    so_i <- seriation::get_order(s_i)
    sampOrder <- c(sampOrder, samps[so_i])
  }
  
  res <- data.frame(
    sampleId = sampOrder,
    RPMMSampleOrder = 1:length(sampOrder)
  )
  colnames(res)[1] <- sampleIdName
  return(res)
}

run_custom_rpmm <- function(dnam, maxlevel=2, sampleIdName="Accession") {
  #'@description Wrapper to run RPMM & extract cluster info
  #'@describeIn Chen et al. Int J Canc 202
    #'
  Y_inv <- t(dnam)
  res <- blcTree(Y_inv, verbose=1, maxlevel=maxlevel)
  print(res)
  plot(res)
  
  cls <- blcTreeLeafClasses(res)
  res_rpmm <- data.frame(table(RPMM=cls, Sample_ID=rownames(Y_inv)))
  res_rpmm <- subset(res_rpmm, Freq != 0)
  res_rpmm$Freq <- NULL
  res_rpmm <- res_rpmm[ , c(2,1)]
  colnames(res_rpmm)[1] <- sampleIdName
  
  ## Extract sample order:
  rownames(res_rpmm) <- as.character(res_rpmm[ , sampleIdName]) #for sample order helper
  samp_order <- get_rpmm_samp_order(res_rpmm, Y_inv, sampleIdName)
  res_rpmm <- merge(res_rpmm, samp_order, by=sampleIdName)
  res_rpmm$Cluster <- substr(res_rpmm$RPMM, 2, 2) #just L or R
  return(res_rpmm)
}

## Load cohort datasets:
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
tcgaGENE <- gbmGENE

load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
dkfzGENE <- gbmGENE

## Execute wrapper:
tcga_rpmm <- run_custom_rpmm(tcgaGENE, 2)
tcga_rpmm$Cohort <- "TCGA"

dkfz_rpmm <- run_custom_rpmm(dkfzGENE, 2)
dkfz_rpmm$Cohort <- "DKFZ"

# write.csv(rbind(tcga_rpmm, dkfz_rpmm), paste0(OUT_DIR, "rpmm_by_cohort.csv"), row.names=FALSE, quote=FALSE)
