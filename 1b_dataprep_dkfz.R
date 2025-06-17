# DKFZ-GBM Data Prep

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

cpg_shared <- read.table(paste0(OUT_DIR,"CpGs_all_shared.csv"))$V1
mGENE <- read.csv(paste0(OUT_DIR,"CpGs_GENE.csv")) #sorted by coord & shared across cohorts

## Clinical covariates:
patients <- read.csv(paste0(DKFZ_DIR, "GSE103659_clinical.csv"), na.strings=STRNAS, stringsAsFactors=FALSE, strip.white=TRUE)
patients <- subset(patients, histology == "Glioblastoma")
patients$sample_title <- as.integer(gsub("SAMPLE ", "", patients$sample_title))
rownames(patients) <- patients$histology <- patients$gradeWHO <- NULL

## Add in inferred sexes from IDAT files:
biol_sexes <- read.csv(paste0(OUT_DIR, "predicted_sex_GSE103659.csv"))
patients <- merge(patients, biol_sexes[,c(1,4)], by="Accession")
patients$sexF <- patients$sex == "F"

## Add in predicted MGMT status:
mgmt <- read.csv(paste0(OUT_DIR,"MGMT_by_cohort.csv"))[ , c("Accession","mMGMT")]
patients <- merge(patients, mgmt, by="Accession")
table(patients$mgmt_status=="methylated", patients$mMGMT) #compare w/ Check with previously published
patients$mgmt_status <- NULL

## Load all CpGs & standardize-Match samples:
gbm450 <- load_450k_from_csv(paste0(DKFZ_DIR, "GSE103659.txt"))
gbm450 <- gbm450[ , colnames(gbm450) %in% patients$Accession]
patients <- subset(patients, Accession %in% colnames(gbm450))
stopifnot( identical(colnames(gbm450), patients$Accession) )

## Impute missing:
100 * mean(is.na(gbm450))
gbm450 <- impute::impute.knn(gbm450)$data

## Subset to focal set:
gbmGENE <- subset(gbm450, rownames(gbm450) %in% mGENE$Name)
gbmGENE <- gbmGENE[match(mGENE$Name, rownames(gbmGENE)), ]
stopifnot(identical(mGENE$Name, rownames(gbmGENE))) #checkpoint

## Subset to shared global set WITHOUT gene of interest & MGMT:
gbm450 <- gbm450[cpg_shared, ]
stopifnot( sum(rownames(gbm450) %in% rownames(gbmGENE)) == 0 ) #anti-dataleakage checkpoint

save(
  list = c("patients", "gbm450", "gbmGENE"),
  file = paste0(OUT_DIR, "dkfzgbm_dnam.RData"),
  compress = TRUE
)
