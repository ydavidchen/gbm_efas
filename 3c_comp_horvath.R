# Compare Horvath Epigenetic Clock & Age by Cluster

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(reshape2)

cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))

horvath <- read.csv(paste0(OUT_DIR,"horvath_by_cohort.csv"))

## TCGA: 
load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
pdat_tcga <- patients
rm(gbm450, gbmGENE, patients)
pdat_tcga <- merge(pdat_tcga, cls, by="Accession")
pdat_tcga <- merge(pdat_tcga, horvath, by="Accession")

## DKFZ:
load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
pdat_dkfz <- patients
rm(gbm450, gbmGENE, patients)
pdat_dkfz <- merge(pdat_dkfz, cls, by="Accession")
pdat_dkfz <- merge(pdat_dkfz, horvath, by="Accession")

## Chronological Age vs. DNAmA
COL_SELE <- c("Accession","Cohort", "Cluster", "age", "horvath.age")
df <- rbind(pdat_tcga[,COL_SELE], pdat_dkfz[,COL_SELE])
mDf <- melt(df)
mDf$Cohort <- forcats::fct_rev(mDf$Cohort)
mDf$variable <- gsub(".age", "", mDf$variable, fixed=TRUE)
mDf$variable <- stringr::str_to_title(mDf$variable)
mDf$variable <- gsub("Age", "Chronological", mDf$variable)

ggplot(mDf, aes(variable, value, color=Cluster)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.25)) +
  facet_wrap(~ Cohort) + 
  ylab("Age Value") +
  THEME_BOX

## Statistical test:
t.test(age ~ Cluster, data=subset(df, Cohort=="TCGA")) #for tableone
t.test(horvath.age ~ Cluster, data=subset(df, Cohort=="TCGA"))

## Statistical test:
t.test(age ~ Cluster, data=subset(df, Cohort=="DKFZ")) #for tableone
t.test(horvath.age ~ Cluster, data=subset(df, Cohort=="DKFZ"))

## TODO: Age Acceleration

