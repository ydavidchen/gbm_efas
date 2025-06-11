# Biological Sex Prediction from IDAT Files
# Notes: Only DKFZ tumors needed

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(minfi)

BASENAME <- paste0(DKFZ_DIR,"GSE103659_RAW/")

# ------------ Prepare Sample Sheet for Selected Samples ------------
idats <- list.files(paste0(DKFZ_DIR,"GSE103659_RAW/"))
idats <- gsub("_Grn.idat", "", idats, fixed=TRUE)
idats <- gsub("_Red.idat", "", idats, fixed=TRUE)
idats <- unique(idats)
idats <- as.data.frame(idats)
idats <- splitstackshape::cSplit(idats, "idats", "_", drop=TRUE)
colnames(idats) <- c("Accession", "Slide", "Array")

patients <- read.csv(paste0(DKFZ_DIR, "GSE103659_clinical.csv"), na.strings=STRNAS, stringsAsFactors=FALSE, strip.white=TRUE)
patients <- subset(patients, histology == "Glioblastoma")

idats <- subset(idats, Accession %in% patients$Accession)
# write.csv(idats, paste0(BASENAME, "gbm_samples.csv"), row.names=FALSE, quote=FALSE)

# ------------ Load DNAm Data Object ------------
TARGETS <- read.metharray.sheet(BASENAME)

RGSet <- read.metharray.exp(base=BASENAME, targets=TARGETS)

MSet <- preprocessRaw(RGSet)

RSet <- ratioConvert(MSet, what="both", keepCN=TRUE)

GRset <- mapToGenome(RSet)

predictedSex <- getSex(GRset, cutoff=-2)#default params

res <- data.frame(
  fullsampleid = predictedSex@rownames,
  sex = predictedSex$predictedSex 
)

res <- splitstackshape::cSplit(res, "fullsampleid", "_")
res <- res[ , c(2,3,4,1)]
colnames(res)[1:3] <- c("Accession", "Slide", "Array")

# write.csv(res, paste0(OUT_DIR, "predicted_sex_GSE103659.csv"), row.names=FALSE, quote=FALSE)
