# Compare Tumor Heterogeneity & Specific Immune Cell types

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)
library(matrixStats)

cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))
rownames(cls) <- cls$Accession

celltypes <- read.csv(paste0(OUT_DIR,"celltypes_by_cohort.csv"), row.names=1, check.names=FALSE)
rownames(celltypes)

eda_wrapper <- function(mat) {
  pheatmap(
    mat,
    annotation_col = cls[ , c("Cohort","Cluster")],
    color = HEAT_COLS,
    fontsize = 10,
    border_color = NA,
    show_rownames = TRUE,
    show_colnames = FALSE,
    cluster_rows = FALSE
  )
}

# ----------------------- Tumor Purity -----------------------
purity <- as.data.frame(t(celltypes["Tumor", ]))
purity <- merge(cls, purity, by.x="Accession", by.y="row.names")

ggplot(purity, aes(Cluster, Tumor)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.25) +
  facet_wrap(~ Cohort) + 
  THEME_BOX

# ----------------------- Tumor MicroEnv (TME) Fraction -----------------------
mat_tme <- celltypes[c("Astrocyte","Dura","Microglia","Neuron","Oligodendrocyte","Endothelium"), ]
eda_wrapper(mat_tme[ , 1:153])
eda_wrapper(mat_tme[ , 154:ncol(mat_tme)])

dfn <- as.data.frame(colSums2(as.matrix(mat_tme)))
colnames(dfn) <- "TME"
dfn <- merge(cls, dfn, by.x="Accession", by.y="row.names")
dfn$Cohort <- forcats::fct_rev(dfn$Cohort)

ggplot(dfn, aes(Cluster, 100*TME)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.25) +
  facet_wrap(~ Cohort) + 
  ylab("Microenvir. Cell %") +
  THEME_BOX

## Univariate tests:
t.test(TME ~ Cluster, data=subset(dfn, Cohort=="TCGA"))
t.test(TME ~ Cluster, data=subset(dfn, Cohort=="DKFZ"))

## Multivariate tests:
wrapper_mult <- function(cohort) {
  if(cohort == "TCGA") {
    load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
  } else if (cohort == "DKFZ") {
    load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
  }
  
  mdf <- merge(patients, dfn, by="Accession")
  n <- nrow(mdf)
  mdf <- mdf[ , c("TME","Cluster","age","sexF","mMGMT")]
  mdf <- mdf[complete.cases(mdf), ]
  mod <- glm(TME ~ (Cluster=="R")+age+sexF+mMGMT, data=mdf, family=gaussian)
  n_sub <- nrow(mdf)
  
  cat("-------------- GLM on ", cohort, "n =", n_sub, "of", n, "--------------")
  print(summary(mod))
  print( rbind(cbind(adjDiff=coef(mod), confint(mod))) )
  
}

wrapper_mult("TCGA")
wrapper_mult("DKFZ")

# ----------------------- Tumor Immune MicroEnv (TME) Fraction -----------------------
## Total Immune infiltration:
tils <- celltypes[c("CD4Tcells","CD8Tcells","Bcells","NKcells","Monocytes","Granulocytes","Neutrophils"), ]
eda_wrapper(tils[ , 1:153])
eda_wrapper(tils[ , 154:ncol(tils)])

dft <- as.data.frame(colSums(tils))
colnames(dft) <- "Immune"
dft <- merge(cls, dft, by.x="Accession", by.y="row.names")
dft$Cohort <- forcats::fct_rev(dft$Cohort)

ggplot(dft, aes(Cluster, Immune)) +
  geom_boxplot(outlier.shape=NA) +
  geom_jitter(width=0.25) +
  facet_wrap(~ Cohort) + 
  THEME_BOX
