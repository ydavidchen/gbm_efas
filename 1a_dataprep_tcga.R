# TCGA-GBM Data Prep

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

cpg_shared <- read.table(paste0(OUT_DIR,"CpGs_all_shared.csv"), header=FALSE)$V2
mGENE <- read.csv(paste0(OUT_DIR,"CpGs_GENE.csv")) #already sorted by coord

# --------------- Part I: Patient- & Sample-level Clinical Metadata ---------------
patients <- read.csv(paste0(TCGA_DIR,"GBM/nationwidechildrens.org_GBM_bio.patient.tsv"), sep="\t", na.strings=c("NA",""))
colnames(patients)[1] <- "patient"
patients$age <- as.integer(patients$age_at_initial_pathologic_diagnosis)
patients <- patients[ , c("patient","age","gender")]
patients$sexF <- patients$gender == "FEMALE"

samples <- read.csv(paste0(TCGA_DIR,"GBM/nationwidechildrens.org_GBM_bio.sample.tsv"), sep="\t")
colnames(samples)[1] <- "Accession"
samples <- subset(samples, sample_type %in% c("Primary Tumor","Recurrent Tumor"))
samples <- samples[ , c("Accession","sample_type")]
samples$patient <- substr(samples$Accession, 1, 12)

patients <- merge(samples, patients, by="patient")
patients$SAMPLE_ID <- substr(patients$Accession, 1, 15) #for cbio

cbio <- read.csv(paste0(CBIO_DIR,"gbm_tcga_pan_can_atlas_2018_clinical_data.tsv"), sep="\t", check.names=FALSE)
cbio$os_months <- as.numeric(cbio$`Overall Survival (Months)`)
cbio$os_status <- as.integer(gsub(":.*", "", cbio$`Overall Survival Status`))
cbio <- cbio[ , c("Sample ID","Subtype","Fraction Genome Altered","os_months","os_status")]

patients <- merge(patients, cbio, by.x="SAMPLE_ID", by.y="Sample ID")
patients$SAMPLE_ID <- NULL

mgmt <- read.csv(paste0(OUT_DIR,"MGMT_by_cohort.csv"))[ , c("Accession","mMGMT")]
patients <- merge(patients, mgmt, by="Accession")

anyDuplicated(patients$Accession)

# --------------- Part II. CpG Matrices ---------------
gbm450 <- load_450k_from_csv(paste0(TCGA_DIR,"GBM/jhu-usc.edu_GBM_HumanMethylation450.betaValue.tsv"))
gbm450 <- gbm450[cpg_shared, ]

## Impute missing:
mean(is.na(gbm450)) * 100
gbm450 <- impute::impute.knn(gbm450)$data

## Standardize sample names:
anyDuplicated(substr(colnames(gbm450), 1, 16))
colnames(gbm450) <- substr(colnames(gbm450), 1, 16)
gbm450 <- gbm450[ , colnames(gbm450) %in% patients$Accession]
dim(gbm450)

patients <- subset(patients, Accession %in% colnames(gbm450))
identical(patients$Accession, colnames(gbm450)) #checkpoint: if not, match

## Subset DNAm for Gene of Interest:
gbmGENE <- subset(gbm450, rownames(gbm450) %in% mGENE$Name)
mGENE <- subset(mGENE, Name %in% rownames(gbm450))
identical(mGENE$Name, rownames(gbmGENE)) #checkpoint: if not, match
gbmGENE <- gbmGENE[match(mGENE$Name, rownames(gbmGENE)), ]

# save(
#   list = c("patients", "gbm450", "gbmGENE"),
#   file = paste0(OUT_DIR, "tcgagbm_dnam.RData"),
#   compress = TRUE
# )
