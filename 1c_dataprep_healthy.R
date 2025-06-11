# PFCs from Non-diseased Controls

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

cpg_shared <- read.table(paste0(OUT_DIR,"CpGs_all_shared.csv"))$V1
mGENE <- read.csv(paste0(OUT_DIR,"CpGs_GENE.csv")) #sorted by coord

## Load & filter clinical data:
patients <- read.csv(paste0(HEALTHY_DIR, "GSE61431_clinical.csv"), na.strings=STRNAS, stringsAsFactors=FALSE, strip.white=TRUE)
patients <- subset(patients, diagnosis=="normal" & tissue=="frontal cortex")
patients$subject_id <- paste0("P", patients$subject_id)
patients <- patients[order(patients$subject_id), ]
rownames(patients) <- patients$sample_title <- NULL

## Load overall DNAm:
pfc450 <- load_450k_from_csv(paste0(HEALTHY_DIR,"GSE61431_dasen.csv.gz"))
pfc450 <- pfc450[ , colnames(pfc450) %in% patients$barcode]
pfc450 <- pfc450[cpg_shared, ]
dim(pfc450)

patients <- subset(patients, barcode %in% colnames(pfc450))
patients <- patients[match(colnames(pfc450), patients$barcode), ]

identical(colnames(pfc450), patients$barcode)
colnames(pfc450) <- patients$subject_id #only possible if 1 disease & 1 tissue type selected!

## Impute missing, if needed:
100 * mean(is.na(pfc450))
pfc450 <- impute::impute.knn(pfc450)$data

## DNAm for Gene of Interest:
pfcGENE <- subset(pfc450, rownames(pfc450) %in% mGENE$Name)
pfcGENE <- pfcGENE[match(mGENE$Name, rownames(pfcGENE)), ]

save(
  list = c("patients", "pfc450", "pfcGENE"),
  file = paste0(OUT_DIR, "healthy_brain_dnam.RData"),
  compress = TRUE
)
