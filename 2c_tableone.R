# Table One by Cluster

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

library(tableone)

COVARS <-  c("age", "sexF", "mMGMT")

cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))

## TCGA: 
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
pdat_tcga <- patients
rm(gbm450, gbmGENE, patients)

pdat_tcga <- merge(pdat_tcga, cls, by="Accession")
CreateTableOne(COVARS, "Cluster", pdat_tcga)

## DKFZ:
load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
pdat_dkfz <- patients
rm(gbm450, gbmGENE, patients)

pdat_dkfz <- merge(pdat_dkfz, cls, by="Accession")
CreateTableOne(COVARS, "Cluster", pdat_dkfz)
