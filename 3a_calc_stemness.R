# Calculate Methylation-based Stemness Index (mDNAsi)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

calc_stemness <- function(mat, wts) {
  #'@description 
    #' Wrapper to implement f(x) = w^T x as a Matrix Operation plus add'l processing
    #' Missng value in `mat` is IGNORED by zero'ing out the signature's contribution
  #'@param mat data matrix where rows=CpGs, columns=sampless
  #'@param wts Weight matrix with up to 219 rows of CpGs & 1 column
  #'@describeIn Malta et al. 2018 Cell STAR Methods
  
  probes_shared <- intersect(rownames(mat), rownames(wts))
  cat("Numbers of probes available:", length(probes_shared), "of", nrow(wts))
  
  mat <- mat[probes_shared, ]
  if(anyNA(mat)) mat[is.na(mat)] <- 0 #zero out contribution of any missing probes
  
  wts <- wts[probes_shared, , drop=FALSE]
  
  stopifnot(identical(rownames(mat), rownames(wts)))
  
  res <- t(wts) %*% mat
  res <- as.data.frame(t(res))
  
  res$mDNAsi <- MinmaxScale(res$Weight)
  return(res)
}

## Published mDNAsi weights:
WEIGHTS <- readxl::read_excel(paste0(DIR, "TCGA_MISC/Cell2018_pancancer_stemness/DNAmethylation_and_RNAexpression_Stemness_Signatures.xlsx"), sheet=1, skip=1)
WEIGHTS <- as.data.frame(WEIGHTS)
WEIGHTS$Weight <- as.numeric(WEIGHTS$Weight)
rownames(WEIGHTS) <- WEIGHTS$`Probe ID`
WEIGHTS$`Probe ID` <- NULL
WEIGHTS <- data.matrix(WEIGHTS)

## Shared loci:
cpg_shared <- read.table(paste0(OUT_DIR,"CpGs_all_shared.csv"))$V1
sum(cpg_shared %in% rownames(WEIGHTS)) #212

## Load DNAm matrices for GBM:
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
tcga450k <- gbm450
rm(gbm450, gbmGENE, patients)

load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
dkfz450k <- gbm450
rm(gbm450, gbmGENE, patients)

## Execute wrapper:
csc_tcga <- calc_stemness(tcga450k, WEIGHTS)
csc_tcga$Cohort <- "TCGA"

csc_dkfz <- calc_stemness(dkfz450k, WEIGHTS)
csc_dkfz$Cohort <- "DKFZ"

# write.csv(rbind(csc_tcga, csc_dkfz), paste0(OUT_DIR, "stemness_by_cohort.csv"), row.names=TRUE, quote=FALSE)

## Verify NO OVERLAP btwn CpGs analyzed vs. signature:
mGENE <- read.csv(paste0(OUT_DIR,"CpGs_GENE.csv")) #sorted by coord
sum(mGENE$Name %in% rownames(WEIGHTS))
