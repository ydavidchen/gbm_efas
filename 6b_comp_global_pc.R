# Compare Global PC Groups

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

## Global methylation in 2D:
glb_pc <- read.csv(paste0(OUT_DIR, "global_pca_by_cohort.csv"))[ ,1:3]
glb_pc$pcGroup <- glb_pc$PC2 > glb_pc$PC1

cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))
cls <- merge(cls, glb_pc, by.x="Accession", by.y="X")
cls$Cohort <- forcats::fct_rev(cls$Cohort)

ggplot(cls, aes(PC1, PC2, color=Cluster)) +
  geom_point(size=3, alpha=0.75) +
  facet_wrap(~ Cohort) + 
  THEME_SCATTER

## Univariate tests:
ctab_tcga <- table(cls$Cluster[cls$Cohort=="TCGA"], cls$pcGroup[cls$Cohort=="TCGA"])
ctab_tcga #inspect!
fisher.test(ctab_tcga)

ctab_dkfz <- table(cls$Cluster[cls$Cohort=="DKFZ"], cls$pcGroup[cls$Cohort=="DKFZ"])
ctab_dkfz #inspect!
fisher.test(ctab_dkfz)

## Multivariate tests:
wrapper_mult <- function(cohort) {
  if(cohort == "TCGA") {
    load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
  } else if (cohort == "DKFZ") {
    load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
  }
  
  mdf <- merge(patients, cls, by="Accession")
  n <- nrow(mdf)
  mdf <- mdf[ , c("pcGroup","Cluster","age","sexF","mMGMT")]
  mdf <- mdf[complete.cases(mdf), ]
  mod <- glm((Cluster=="R")~pcGroup+age+sexF+mMGMT, data=mdf, family=binomial)
  n_sub <- nrow(mdf)
  
  cat("-------------- GLM on ", cohort, "n =", n_sub, "of", n, "--------------")
  print(summary(mod))
  print( exp( rbind(cbind(aOR=coef(mod), confint(mod))) ) )
  
}

wrapper_mult("TCGA")
wrapper_mult("DKFZ")
