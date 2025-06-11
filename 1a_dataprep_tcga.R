# TCGA-GBM Data Prep

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

cpg_shared <- read.table(paste0(OUT_DIR,"CpGs_all_shared.csv"))$V1
mGENE <- read.csv(paste0(OUT_DIR,"CpGs_GENE.csv")) #already sorted by coord

# --------------- Part I: Patient- & Sample-level Clinical Metadata ---------------
## Patient-level annotations:
patients <- read.csv(paste0(TCGA_DIR,"GBM/nationwidechildrens.org_GBM_bio.patient.tsv"), sep="\t", na.strings=c("NA",""))
colnames(patients)[1] <- "patient"
patients$age <- as.integer(patients$age_at_initial_pathologic_diagnosis)
patients <- patients[ , c("patient","age","gender")]
patients$sexF <- patients$gender == "FEMALE"

## Sample-level annotations:
samples <- read.csv(paste0(TCGA_DIR,"GBM/nationwidechildrens.org_GBM_bio.sample.tsv"), sep="\t")
colnames(samples)[1] <- "Accession"
samples <- subset(samples, sample_type %in% c("Primary Tumor","Recurrent Tumor"))
samples <- samples[ , c("Accession","sample_type")]
samples$patient <- substr(samples$Accession, 1, 12)
patients <- merge(samples, patients, by="patient")

## Survival data BACKUP: DO NOT DELETE
# patients$SAMPLE_ID <- substr(patients$Accession, 1, 15) #for cbio merge
# cbio <- read.csv(paste0(CBIO_DIR,"gbm_tcga_pan_can_atlas_2018_clinical_data.tsv"), sep="\t", check.names=FALSE)
# cbio$os_months <- as.numeric(cbio$`Overall Survival (Months)`)
# cbio$os_status <- as.integer(gsub(":.*", "", cbio$`Overall Survival Status`))
# cbio <- cbio[ , c("Sample ID","os_months","os_status")]
# patients <- merge(patients, cbio, by.x="SAMPLE_ID", by.y="Sample ID")
# patients$SAMPLE_ID <- NULL

## Survival data DO NOT REMOVE
tcgacdr <- readxl::read_excel(paste0(DIR,"TCGA_MISC/Liu2018CellTCGACDR.xlsx"), sheet=1)
tcgacdr <- as.data.frame(tcgacdr)
tcgacdr <- subset(tcgacdr, type=="GBM")
tcgacdr$os_months <- tcgacdr$OS.time / 30.436 #convert from days to months
tcgacdr$os_status <- as.integer(tcgacdr$OS)
tcgacdr <- tcgacdr[ , c("bcr_patient_barcode","os_status","os_months")]
patients <- merge(patients, tcgacdr, by.x="patient", by.y="bcr_patient_barcode", all.x=TRUE)

## MGMT status predictions
mgmt <- read.csv(paste0(OUT_DIR,"MGMT_by_cohort.csv"))[ , c("Accession","mMGMT")]
patients <- merge(patients, mgmt, by="Accession", all.x=TRUE)

anyDuplicated(patients$Accession)

# --------------- Part II. CpG Matrices ---------------
gbm450 <- load_450k_from_csv(paste0(TCGA_DIR,"GBM/jhu-usc.edu_GBM_HumanMethylation450.betaValue.tsv"))
gbm450 <- gbm450[cpg_shared, ]

## Standardize sample names:
anyDuplicated(substr(colnames(gbm450), 1, 16))
colnames(gbm450) <- substr(colnames(gbm450), 1, 16)
gbm450 <- gbm450[ , colnames(gbm450) %in% patients$Accession]
dim(gbm450)

patients <- subset(patients, Accession %in% colnames(gbm450))
identical(patients$Accession, colnames(gbm450)) #checkpoint: if not, match

## Check missing values in covariates
# tmp <- patients[ , c("age","sexF","os_status")]
# tmp <- tmp[complete.cases(tmp), ]

## Impute missing:
mean(is.na(gbm450)) * 100
gbm450 <- impute::impute.knn(gbm450)$data

## Subset DNAm for Gene of Interest:
gbmGENE <- subset(gbm450, rownames(gbm450) %in% mGENE$Name)
mGENE <- subset(mGENE, Name %in% rownames(gbm450))
identical(mGENE$Name, rownames(gbmGENE)) #checkpoint: if not, match
gbmGENE <- gbmGENE[match(mGENE$Name, rownames(gbmGENE)), ]

save(
  list = c("patients", "gbm450", "gbmGENE"),
  file = paste0(OUT_DIR, "tcgagbm_dnam.RData"),
  compress = TRUE
)
