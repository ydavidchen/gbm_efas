# Compare Focal Epigenetic Dysregulation

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

## Load computed profiles:
cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))
mdi <- read.csv(paste0(OUT_DIR, "mdi_by_cohort.csv"))
mdi$Cohort <- NULL

## Visualization & testing:
df <- merge(cls, mdi, by="Accession")

ggplot(df, aes(Cluster, 100*MDI)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.25) +
  facet_wrap(~ Cohort) +
  ylab(paste("% Methyl. Dysregulation,", GENE)) +
  THEME_BOX

t.test(MDI ~ Cluster, data=subset(df, Cohort=="TCGA"))
t.test(MDI ~ Cluster, data=subset(df, Cohort=="DKFZ"))
