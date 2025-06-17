# TCGA-GBM Data Prep

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

cpg_shared <- read.table(paste0(OUT_DIR,"CpGs_all_shared.csv"))$V1
mGENE <- read.csv(paste0(OUT_DIR,"CpGs_GENE.csv")) #sorted by coord & shared across cohorts

# --------------- Part I: Patient- & Sample-level Clinical Metadata ---------------
## Patient-level annotations:
patients <- read.csv(paste0(TCGA_DIR,"GBM/nationwidechildrens.org_GBM_bio.patient.tsv"), sep="\t", na.strings=STRNAS)
colnames(patients)[1] <- "patient"
patients$age <- as.integer(patients$age_at_initial_pathologic_diagnosis)
patients <- patients[ , c("patient","age","gender")]
patients$sexF <- patients$gender == "FEMALE"

## Sample-level annotations:
samples <- read.csv(paste0(TCGA_DIR,"GBM/nationwidechildrens.org_GBM_bio.sample.tsv"), sep="\t", na.strings=STRNAS)
colnames(samples)[1] <- "Accession"
samples <- subset(samples, sample_type %in% c("Primary Tumor","Recurrent Tumor"))
samples <- samples[ , c("Accession","sample_type")]
samples$patient <- substr(samples$Accession, 1, 12)
patients <- merge(samples, patients, by="patient")

## Survival data DO NOT REMOVE
tcgacdr <- readxl::read_excel(paste0(DIR,"TCGA_MISC/Liu2018CellTCGACDR.xlsx"), sheet=1)
tcgacdr <- as.data.frame(tcgacdr)
tcgacdr <- subset(tcgacdr, type=="GBM")
tcgacdr$os_months <- tcgacdr$OS.time / 30.436875
tcgacdr$os_status <- as.integer(tcgacdr$OS)
tcgacdr <- tcgacdr[ , c("bcr_patient_barcode","os_status","os_months")]
patients <- merge(patients, tcgacdr, by.x="patient", by.y="bcr_patient_barcode", all.x=TRUE)

## MGMT status predicted:
mgmt <- read.csv(paste0(OUT_DIR,"MGMT_by_cohort.csv"))[ , c("Accession","mMGMT")]
patients <- merge(patients, mgmt, by="Accession")

anyDuplicated(patients$Accession)

# --------------- Part II. CpG Matrices ---------------
gbm450 <- load_450k_from_csv(paste0(TCGA_DIR,"GBM/jhu-usc.edu_GBM_HumanMethylation450.betaValue.tsv"))

## Standardize sample names:
anyDuplicated(substr(colnames(gbm450), 1, 16)) #required
colnames(gbm450) <- substr(colnames(gbm450), 1, 16)
gbm450 <- gbm450[ , colnames(gbm450) %in% patients$Accession]
dim(gbm450)

patients <- subset(patients, Accession %in% colnames(gbm450))
stopifnot( identical(patients$Accession, colnames(gbm450)) )

## Impute missing:
mean(is.na(gbm450)) * 100
gbm450 <- impute::impute.knn(gbm450)$data

## Subset to focal set:
gbmGENE <- subset(gbm450, rownames(gbm450) %in% mGENE$Name)
gbmGENE <- gbmGENE[match(mGENE$Name, rownames(gbmGENE)), ]
stopifnot(identical(mGENE$Name, rownames(gbmGENE))) #checkpoint

## Subset to shared global set WITHOUT gene of interest & MGMT:
gbm450 <- gbm450[cpg_shared, ]
stopifnot( sum(rownames(gbm450) %in% rownames(gbmGENE)) == 0 ) #anti-dataleakage checkpoint

## Subset to Gene of Interest:
save(
  list = c("patients", "gbm450", "gbmGENE"),
  file = paste0(OUT_DIR, "tcgagbm_dnam.RData"),
  compress = TRUE
)
