# DNAm EDA

# TODO-Clinical:
# 1. (both) MGMT27stp
# 2. sex prediction before chrXY dropping: TCGA (eval), DKFZ (apply)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)
library(survival)
library(survminer)

# ---------- USER CHOICE----------
COHORT <- "TCGA"
# COHORT <- "DKFZ"
# --------------------------------

# RPMM results:
cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))

if(COHORT == "TCGA") {
  load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
  cls <- merge(patients, cls, by.x="sample", by.y="Accession") #TCGA
} else if (COHORT == "DKFZ") {
  load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
  cls <- merge(patients, cls, by="Accession") #DKFZ
}
rm(gbm450)

table(cls$Cluster)

# ------------------------- Part I. Cluster Structure ------------------------- 
## Heatmap annotation:
if(COHORT == "TCGA") {
  hm_annot_samps <- data.frame(
    row.names = patients$sample,
    IDH = ifelse(patients$IDH, "Yes", "No"),
    SexF = ifelse(patients$gender=="FEMALE", "Yes", "No")
  )
} else if (COHORT == "DKFZ") {
  hm_annot_samps <- data.frame(
    row.names = patients$Accession,
    mMGMT = ifelse(patients$mgmt_status == "methylated", "Yes", "No")
  )
}
hm_annot_samps$AgeHi <- ifelse(patients$age > median(patients$age, na.rm=TRUE), "Yes", "No")

hm_cols <- list(
  Txn = c(TSS="black", UTR="darkgray", Body="lightgray"),
  CGI = c(Island="black", Shore="darkgray", Shelf="lightgray")
  
)
hm_cols[["AgeHi"]] <- hm_cols[["SexF"]] <- hm_cols[["IDH"]] <- hm_cols[["mMGMT"]] <- hm_cols[["cerebellum"]] <- BINARY_COLORS

mat <- t(gbmFASN)
mat <- mat[cls$RPMMSampleOrder, ]

## UHC:
GAP_COL <- list(TCGA=c(78), DKFZ=c(116))
pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cutree_rows = 3,
  clustering_method = CL_PARAMS[1],
  clustering_distance_rows = CL_PARAMS[2],
  gaps_row = GAP_COL[[COHORT]],
  # annotation_row = hm_annot_samps,
  # annotation_col = hm_loci_annot,
  annotation_colors = hm_cols,
  color = HEAT_COLS,
  fontsize = 10,
  border_color = NA,
  show_rownames = FALSE,
  show_colnames = TRUE
)

# ------------------------- Part II. Survival with Admin Censoring ------------------------- 
MO_CENSOR <- 60

## Overall Survival:
cls$os_months[cls$os_months > MO_CENSOR] <- MO_CENSOR
cls$os_status[cls$os_months > MO_CENSOR] <- 0

binaryos <- survfit(Surv(os_months, os_status)~Cluster, data=cls)
ggsurvplot(binaryos, cls, palette=SURV_COLS, ggtheme=THEME_SURV, pval=FALSE, risk.table=FALSE)

# TODO: Cox Model

## Progression Free Survival (EDA / Supplemental):
cls$pfs_months[cls$pfs_months > MO_CENSOR] <- MO_CENSOR
cls$pfs_status[cls$pfs_status > MO_CENSOR] <- 0

binaryps <- survfit(Surv(pfs_months, pfs_status)~Cluster, data=cls)
ggsurvplot(binaryps, cls, palette=SURV_COLS, ggtheme=THEME_SURV, pval=FALSE, risk.table=FALSE)
