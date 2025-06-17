# Recursively Partitioned Mixture Model (RPMM) Clustering of Locus DNAm Data

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(RPMM)
library(seriation)

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
