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
t1_tcga <- CreateTableOne(COVARS, "Cluster", pdat_tcga, testNormal=t.test, testExact=fisher.test)
t1_tcga
# write.csv(print(t1_tcga, showAllLevels=TRUE), paste0(OUT_DIR,"tableone_tcga.csv"))

## DKFZ:
load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
pdat_dkfz <- patients
rm(gbm450, gbmGENE, patients)

pdat_dkfz <- merge(pdat_dkfz, cls, by="Accession")
t1_dkfz <- CreateTableOne(COVARS, "Cluster", pdat_dkfz, testNormal=t.test, testExact=fisher.test)
t1_dkfz
# write.csv(print(t1_dkfz, showAllLevels=TRUE), paste0(OUT_DIR, "tableone_dkfz.csv"))
