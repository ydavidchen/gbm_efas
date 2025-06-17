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
COL_SELE <- c("Accession","Cohort", "sexF", "Cluster", "age", "horvath.age")
df <- rbind(pdat_tcga[,COL_SELE], pdat_dkfz[,COL_SELE])

# ----------------------------- Uni & Multivariate Tests -----------------------------
## TCGA:
t.test(age ~ Cluster, data=subset(df, Cohort=="TCGA")) #show in Table 1

t.test(horvath.age ~ Cluster, data=subset(df, Cohort=="TCGA"))
summary(glm(horvath.age ~ Cluster+sexF+mMGMT, data=subset(df, Cohort=="TCGA")))

## DKFZ
t.test(age ~ Cluster, data=subset(df, Cohort=="DKFZ")) #show in Table 1

t.test(horvath.age ~ Cluster, data=subset(df, Cohort=="DKFZ"))
summary(glm(horvath.age ~ Cluster+sexF+mMGMT, data=subset(df, Cohort=="DKFZ")))

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

## Age Acceleration:
ggplot(df, aes(age, horvath.age, color=Cluster)) +
  geom_point() +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~ Cohort) + 
  THEME_SCATTER
