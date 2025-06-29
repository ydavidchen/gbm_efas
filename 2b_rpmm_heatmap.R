# Heatmap Visualization

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)

mGENE <- read.csv(paste0(OUT_DIR,"CpGs_GENE.csv"), row.names=1)
mGENE$Context <- gsub("[N|S]_", "", mGENE$Relation_to_UCSC_CpG_Island)

ANNOT_COLORS <- list(
  Cluster = c(L="lightgray", R="black"),
  RPMM = c(LL="darkgreen", LR="lightgreen", RL="magenta", RR="purple"),
  Context = c(Island="black", Shore="darkgray", Shelf="lightgray", OpenSea="lightblue")
)

wrapper_hm <- function(dnam, gaps_row=NULL, annot_row=NULL) {
  #'@param dnam rows=CpGs ordered by position, cols=samples ordered by RPMM
  pheatmap(
    t(dnam),
    cluster_rows = FALSE, 
    cluster_cols = FALSE,
    gaps_row = gaps_row,
    annotation_row = annot_row,
    annotation_col = mGENE[ , "Context", drop=FALSE],
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

## Main figure:
table(cls_tcga$Cluster)[1] #61
table(cls_dkfz$Cluster)[1] #109

wrapper_hm(tcgaGENE, 61, cls_tcga[,"Cluster",drop=FALSE]) 
wrapper_hm(dkfzGENE, 109, cls_dkfz[,"Cluster",drop=FALSE]) 

## Supplemental figure with terminal RPMM solutions:
cumsum(table(cls_tcga$RPMM)[1:3]) #31  61 128 
cumsum(table(cls_dkfz$RPMM)[1:3]) #60 109 147

wrapper_hm(tcgaGENE, c(31, 61, 128), cls_tcga[,c("RPMM","Cluster")])
wrapper_hm(dkfzGENE, c(60, 109, 147), cls_dkfz[,c("RPMM","Cluster")])
