# Compare Horvath Epigenetic Clock & Age by Cluster

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(reshape2)

# ----------------------------- Preparation -----------------------------
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

## Shared object for tests:
COL_SELE <- c("Accession","Cohort", "sexF", "Cluster", "age", "horvath.age","mMGMT")
df <- rbind(pdat_tcga[,COL_SELE], pdat_dkfz[,COL_SELE])

wrapper <- function(df) {
  mod <- glm(horvath.age ~ (Cluster=="R")+sexF+mMGMT, data=df, family=gaussian)
  print(summary(mod))
  print( rbind(cbind(adjDiff=coef(mod), confint(mod))) )
}

# ----------------------------- Uni & Multivariate Tests -----------------------------
## TCGA:
t.test(horvath.age ~ Cluster, data=subset(df, Cohort=="TCGA"))
wrapper(subset(df, Cohort=="TCGA"))

## DKFZ
t.test(horvath.age ~ Cluster, data=subset(df, Cohort=="DKFZ"))
wrapper(subset(df, Cohort=="DKFZ"))

# ----------------------------- Data Visualization -----------------------------
## Combined Boxplots:
mDf <- melt(df)
mDf$Cohort <- forcats::fct_rev(mDf$Cohort)
mDf$variable <- gsub("horvath.age", "Epigenetic", mDf$variable, fixed=TRUE)
mDf$variable <- gsub("age", "Chronological", mDf$variable)

ggplot(mDf, aes(variable, value, color=Cluster)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.25)) +
  facet_wrap(~ Cohort) + 
  ylab("Age") +
  THEME_BOX
