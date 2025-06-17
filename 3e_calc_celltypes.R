#  Brain-specific TME & TIME Deconvolution

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(MDBrainT)
library(nnls) #dependency
library(pheatmap)

## Load annotation & reference data:
data(brain_sig_matrix)
data(celltype_anno)
cpg_shared <- read.table(paste0(OUT_DIR,"CpGs_all_shared.csv"))$V1 #excludes gene of interest & MGMT

## Get universe & subset reference data:
sum(celltype_anno$id %in% cpg_shared)
100 * mean(celltype_anno$id %in% cpg_shared)

univ <- intersect(celltype_anno$id, cpg_shared)

brain_sig_matrix <- brain_sig_matrix[univ, ]
celltype_anno <- celltype_anno[univ, ]

## Load & subset datasets:
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
tcga450k <- gbm450[univ, ]
rm(gbm450, gbmGENE, patients)

load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
dkfz450k <- gbm450[univ, ]
rm(gbm450, gbmGENE, patients)

## Execute wrapper:
wrapper <- function(dnam) {
  res <- TMEdeconvolute(
    sig_matrix = brain_sig_matrix,
    cell_annotation = celltype_anno,
    mixture_file = as.data.frame(dnam)
  )
  
  ph <- pheatmap(
    res[-nrow(res), ],
    color = HEAT_COLS,
    fontsize = 10,
    border_color = NA,
    show_rownames = TRUE,
    show_colnames = FALSE,
    cluster_rows = FALSE
  )
  print(ph)
  
  return(res)
}

celltypes_tcga <- wrapper(tcga450k)
celltypes_dkfz <- wrapper(dkfz450k)

# write.csv(cbind(celltypes_tcga, celltypes_dkfz), paste0(OUT_DIR,"celltypes_by_cohort.csv"), row.names=TRUE, quote=FALSE)
