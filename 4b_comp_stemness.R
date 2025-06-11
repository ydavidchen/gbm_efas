# Compare Stemness Index between Clusters
# Note sample size reduction in TCGA due to missing covariates

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

## Visualization:
cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))

stem <- read.csv(paste0(OUT_DIR, "stemness_by_cohort.csv"))
stem$Cohort <- NULL

df <- merge(cls, stem, by.x="Accession", by.y="X")
df$Cohort <- forcats::fct_rev(df$Cohort)

ggplot(df, aes(Cluster, mDNAsi)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.25) +
  facet_wrap(~ Cohort) + 
  THEME_BOX

## Univariate:
t.test(mDNAsi ~ Cluster, data=subset(df, Cohort=="TCGA"))
t.test(mDNAsi ~ Cluster, data=subset(df, Cohort=="DKFZ"))

## Multivariate:
wrapper_mult <- function(cohort) {
  if(cohort == "TCGA") {
    load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
  } else if (cohort == "DKFZ") {
    load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
  }
  
  mdf <- merge(patients, df, by="Accession")
  n <- nrow(mdf)
  mdf <- mdf[ , c("mDNAsi","Cluster","age","sexF")]
  mdf <- mdf[complete.cases(mdf), ]
  mod <- glm(mDNAsi ~ (Cluster=="R")+age+sexF, data=mdf, family=gaussian)
  n_sub <- nrow(mdf)
  
  cat("-------------- GLM on ", cohort, "n =", n_sub, "of", n, "--------------")
  print(summary(mod))
  print( rbind(cbind(adjDiff=coef(mod), confint(mod))) )
  
}

wrapper_mult("TCGA")
wrapper_mult("DKFZ")
