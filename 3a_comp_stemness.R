# Compare Stemness Index between Clusters

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

## Load computed profiles:
cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))
stem <- read.csv(paste0(OUT_DIR, "stemness_by_cohort.csv"))
stem$Cohort <- NULL

df <- merge(cls, stem, by.x="Accession", by.y="X")

ggplot(df, aes(Cluster, mDNAsi)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.25) +
  facet_wrap(~ Cohort) + 
  THEME_BOX

t.test(mDNAsi ~ Cluster, data=subset(df, Cohort=="TCGA"))
t.test(mDNAsi ~ Cluster, data=subset(df, Cohort=="DKFZ"))

