# Heatmap Visualization

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)

ANNOT_COLORS <- list(
  Cluster = c(L="lightgray", R="black"),
  RPMM = c(LL="darkgreen", LR="lightgreen", RL="magenta", RR="purple"),
  Txn = c(TSS="black", UTR="darkgray", Body="lightgray"),
  CGI = c(Island="black", Shore="darkgray", Shelf="lightgray")
)

wrapper_hm <- function(dnam, gaps_row=NULL, annot_row=NULL) {
  #'@param dnam rows=CpGs ordered by position, cols=samples ordered by RPMM
  pheatmap(
    t(dnam),
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    gaps_row = gaps_row,
    annotation_row = annot_row,
    annotation_colors = ANNOT_COLORS,
    color = HEAT_COLS,
    fontsize = 10,
    border_color = NA,
    show_rownames = FALSE,
    show_colnames = TRUE
  )
}

## Load cohort datasets:
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
tcgaGENE <- gbmGENE
rm(gbm450, gbmGENE, patients)

load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
dkfzGENE <- gbmGENE
rm(gbm450, gbmGENE, patients)

## RPMM results w/ sample order:
cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))
cls$RPMM <- gsub("^r", "", cls$RPMM, ignore.case=FALSE)
rownames(cls) <- cls$Accession #for heatmap

cls_tcga <- subset(cls, Cohort == "TCGA")
cls_tcga <- cls_tcga[order(cls_tcga$RPMMSampleOrder), ]
tcgaGENE <- tcgaGENE[ , cls_tcga$Accession]

cls_dkfz <- subset(cls, Cohort == "DKFZ")
cls_dkfz <- cls_dkfz[order(cls_dkfz$RPMMSampleOrder), ]
dkfzGENE <- dkfzGENE[ , cls_dkfz$Accession]

## Main visualization:
wrapper_hm(tcgaGENE, 73, cls_tcga[,"Cluster",drop=FALSE]) #c(11, 74, 124)-1
wrapper_hm(dkfzGENE, 115, cls_dkfz[,"Cluster",drop=FALSE]) #c(63, 116, 149)-1

## Supplemental visualization:
wrapper_hm(tcgaGENE, c(11, 74, 124)-1, cls_tcga[,c("RPMM","Cluster")])
wrapper_hm(dkfzGENE, c(63, 116, 149)-1, cls_dkfz[,c("RPMM","Cluster")])

