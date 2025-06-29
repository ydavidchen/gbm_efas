# Differential Methylation

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(limma)

SELE_VARS <- c("Accession", "Cluster", "age", "sexF", "mMGMT")

cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))

sele_sites <- read.csv(paste0(OUT_DIR, "global_dnam_CpGs.txt"), col.names=0)$X0

ANNOT <- load_450k_annot()
ANNOT <- subset(ANNOT, Name %in% sele_sites)
rownames(ANNOT) <- ANNOT$Name #for limma output
ANNOT$Name <- NULL

run_limma <- function(betavals, design_mat, varname="ClusterR", padj=0.05, lfc=0) {
  #'@description Wrapper to run LIMMA using M-values
  #'@param betavals Input matrix of beta-values, which will be transformed
  stopifnot( identical(colnames(betavals), rownames(design_mat)) )
  
  print("Transforming to M-values")
  mvals <- lumi::beta2m(betavals)
  
  fit <- lmFit(mvals, design=design_mat)
  fit <- eBayes(fit)
  
  diff_sites <- topTable(fit, number=Inf, coef=varname, genelist=ANNOT, adjust.method="bonferroni", sort.by="p")
  diff_sites$isSig <- diff_sites$adj.P.Val < padj
  
  diff_sites$direction <- ""
  diff_sites$direction[diff_sites$isSig & diff_sites$logFC > lfc] <- "UP"
  diff_sites$direction[diff_sites$isSig & diff_sites$logFC < -lfc] <- "DOWN"
  
  return(diff_sites)
}

# ---------------------------- TCGA ----------------------------
load(paste0(OUT_DIR, "tcgagbm_dnam.RData")) #covariates not complete
stopifnot(identical(patients$Accession, colnames(gbm450))) #important checkpoint
pdata_tcga <- merge(patients, cls, by="Accession")
pdata_tcga <- pdata_tcga[ , SELE_VARS]
pdata_tcga <- pdata_tcga[complete.cases(pdata_tcga), ]
tcga450k <- gbm450[ , pdata_tcga$Accession]
stopifnot(identical(pdata_tcga$Accession, colnames(tcga450k)))
rm(gbm450, gbmGENE, patients)

tcga450k <- tcga450k[sele_sites, ]

dm_tcga <- model.matrix(~ Cluster + age + sexF + mMGMT, data=pdata_tcga)
rownames(dm_tcga) <- pdata_tcga$Accession

ewas_tcga <- run_limma(tcga450k, dm_tcga)
table(ewas_tcga$isSig)
table(ewas_tcga$direction)
100*mean(ewas_tcga$isSig)

# write.csv(ewas_tcga, file=paste0(OUT_DIR,"ewas_tcga.csv"), row.names=TRUE, quote=FALSE)

# ---------------------------- DKFZ ----------------------------
load(paste0(OUT_DIR, "dkfzgbm_dnam.RData")) #data complete
stopifnot(identical(patients$Accession, colnames(gbm450))) #important checkpoint
pdata_dkfz <- merge(patients, cls, by="Accession")
pdata_dkfz <- pdata_dkfz[ , SELE_VARS]
dkfz450k <- gbm450
stopifnot(identical(pdata_dkfz$Accession, colnames(dkfz450k)))
rm(gbm450, gbmGENE, patients)

dkfz450k <- dkfz450k[sele_sites, ]

dm_dkfz <- model.matrix(~ Cluster + age + sexF + mMGMT, data=pdata_dkfz)
rownames(dm_dkfz) <- pdata_dkfz$Accession

ewas_dkfz <- run_limma(dkfz450k, dm_dkfz)
table(ewas_dkfz$isSig)
table(ewas_dkfz$direction)
100*mean(ewas_dkfz$isSig)

# write.csv(ewas_dkfz, file=paste0(OUT_DIR,"ewas_dkfz.csv"), row.names=TRUE, quote=FALSE)

# ---------------------------- Common Output Prep ----------------------------
cpgs_down <- intersect(
  rownames(ewas_tcga)[ewas_tcga$direction=="DOWN"], 
  rownames(ewas_dkfz)[ewas_dkfz$direction=="DOWN"]
)
length(cpgs_down) #1182
# write.table(cpgs_down, paste0(OUT_DIR, "ewas_combined_down_cpgs.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

## Input files for WebGestalt KEGG:
genes_down <- parse_distinct_refgenes(ANNOT[cpgs_down, "UCSC_RefGene_Name"])
length(genes_down) #904
# write.table(genes_down, paste0(OUT_DIR, "ewas_combined_down_genes.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)

genes_bg <- parse_distinct_refgenes(ANNOT$UCSC_RefGene_Name)
length(genes_bg)
# write.table(genes_bg, paste0(OUT_DIR, "ewas_background_genes.txt"), row.names=FALSE, col.names=FALSE, quote=FALSE)
