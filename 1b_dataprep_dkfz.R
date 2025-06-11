# DKFZ-GBM Data Prep

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

cpg_shared <- read.table(paste0(OUT_DIR,"CpGs_all_shared.csv"))$V1
mGENE <- read.csv(paste0(OUT_DIR,"CpGs_GENE.csv")) #already sorted by coord

## Clinical covariates:
patients <- read.csv(paste0(DKFZ_DIR, "GSE103659_clinical.csv"), na.strings=STRNAS, stringsAsFactors=FALSE, strip.white=TRUE)
patients <- subset(patients, histology == "Glioblastoma")
patients$sample_title <- as.integer(gsub("SAMPLE ", "", patients$sample_title))
rownames(patients) <- patients$histology <- patients$gradeWHO <- NULL

## Add in predicted sexes from IDAT files:
biol_sexes <- read.csv(paste0(OUT_DIR, "predicted_sex_GSE103659.csv"))
patients <- merge(patients, biol_sexes[,c(1,4)], by="Accession")
patients$sexF <- patients$sex == "F"

## Add in newly predicted MGMT status:
mgmt <- read.csv(paste0(OUT_DIR,"MGMT_by_cohort.csv"))[ , c("Accession","mMGMT")]
patients <- merge(patients, mgmt, by="Accession")
table(patients$mgmt_status=="methylated", patients$mMGMT) #compare w/ Check with previously published
patients$mgmt_status <- NULL

## Load all CpGs:
gbm450 <- load_450k_from_csv(paste0(DKFZ_DIR, "GSE103659.txt"))
gbm450 <- gbm450[cpg_shared, ]
dim(gbm450)

## Standardize-Match samples:
patients <- subset(patients, Accession %in% colnames(gbm450))
gbm450 <- gbm450[ , colnames(gbm450) %in% patients$Accession]
identical(colnames(gbm450), patients$Accession) #checkpoint

## Impute missing, if needed:
100 * mean(is.na(gbm450))
gbm450 <- impute::impute.knn(gbm450)$data

## Subset DNAm to Gene of Interest:
gbmGENE <- subset(gbm450, rownames(gbm450) %in% mGENE$Name)
mGENE <- subset(mGENE, Name %in% rownames(gbmGENE))
gbmGENE <- gbmGENE[match(mGENE$Name, rownames(gbmGENE)), ]
identical(mGENE$Name, rownames(gbmGENE)) #checkpoint

save(
  list = c("patients", "gbm450", "gbmGENE"),
  file = paste0(OUT_DIR, "dkfzgbm_dnam.RData"),
  compress = TRUE
)
