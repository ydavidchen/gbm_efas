# Compare Stemness between Clusters

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

# ---------- USER CHOICE: TOGGLE SWITCH ----------
COHORT <- "TCGA"
# COHORT <- "DKFZ"
# ------------------------------------------------

csc <- read.csv(paste0(OUT_DIR, "stemness_by_cohort.csv"))
cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))

csc <- subset(csc, Cohort == COHORT)
cls <- subset(cls, Cohort == COHORT)
csc$Cohort <- cls$Cohort <- NULL

if(COHORT=="TCGA") csc$X <- substr(csc$X, 1, 16)

df <- merge(cls, csc, by.x="Accession", by.y="X")

ggplot(df, aes(Cluster, mDNAsi)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.25) +
  THEME_BOX

t.test(mDNAsi ~ Cluster, data=df)
wilcox.test(mDNAsi ~ Cluster, data=df)
