# MGMT Promoter Methylation Status
# Notes:
# 1. Obtain the necessary input CpGs from scratch, i.e. non-filtered datasets
# 2. Both CpGs needed for prediction are available.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(mgmtstp27)

PROBES <- c("cg12981137","cg12434587")

wrapper_mgmt <- function(dnam, probes=PROBES) {
  dnam <- dnam[probes, ]
  dnam <- lumi::beta2m(dnam)
  dnam <- as.data.frame(t(dnam))
  res <- MGMTpredict(dnam)
  res$mMGMT <- res$state == "M"
  rownames(res) <- NULL
  colnames(res)[1] <- "Accession"
  return(res)
}

tcga450k <- load_450k_from_csv(paste0(TCGA_DIR,"GBM/jhu-usc.edu_GBM_HumanMethylation450.betaValue.tsv"))
dkfz450k <- load_450k_from_csv(paste0(DKFZ_DIR, "GSE103659.txt"))

## Execute wrapper:
mgmt_tcga <- wrapper_mgmt(tcga450k)
mgmt_tcga$Accession <- substr(mgmt_tcga$Accession, 1, 16)
mgmt_tcga$Cohort <- "TCGA"

mgmt_dkfz <- wrapper_mgmt(dkfz450k)
mgmt_dkfz$Cohort <- "DKFZ"

# write.csv(rbind(mgmt_tcga, mgmt_dkfz), file=paste0(OUT_DIR,"MGMT_by_cohort.csv"), row.names=FALSE, quote=FALSE)
