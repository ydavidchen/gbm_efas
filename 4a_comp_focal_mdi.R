# Compare Focal Epigenetic Alteration

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

## Visualization:
cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))
mdi <- read.csv(paste0(OUT_DIR, "mdi_by_cohort.csv"))
mdi$Cohort <- NULL

df <- merge(cls, mdi, by="Accession")
df$Cohort <- forcats::fct_rev(df$Cohort)

ggplot(df, aes(Cluster, 100*MDI)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.25) +
  facet_wrap(~ Cohort) +
  ylab(paste("% Epi Alteration at", GENE)) +
  THEME_BOX

## Univariate tests:
t.test(MDI ~ Cluster, data=subset(df, Cohort=="TCGA"))
t.test(MDI ~ Cluster, data=subset(df, Cohort=="DKFZ"))

## Multivariate tests:
wrapper_mult <- function(cohort) {
  if(cohort == "TCGA") {
    load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
  } else if (cohort == "DKFZ") {
    load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
  }
  
  mdf <- merge(patients, df, by="Accession")
  n <- nrow(mdf)
  mdf <- mdf[ , c("MDI","Cluster","age","sexF","mMGMT")]
  mdf <- mdf[complete.cases(mdf), ]
  mod <- glm(MDI ~ (Cluster=="R")+age+sexF+mMGMT, data=mdf, family=gaussian)
  n_sub <- nrow(mdf)
  
  cat("-------------- GLM on ", cohort, "n =", n_sub, "of", n, "--------------")
  print(summary(mod))
  print( rbind(cbind(adjDiff=coef(mod), confint(mod))) )
  
}

wrapper_mult("TCGA")
wrapper_mult("DKFZ")
