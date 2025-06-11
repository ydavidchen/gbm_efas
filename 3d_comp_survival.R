# Overall Survival: Kaplan-Meier & CoxPH

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(pheatmap)
library(survival)
library(survminer)

MO_CENSOR <- 60

# ---------- USER CHOICE----------
COHORT <- "TCGA"
# COHORT <- "DKFZ"
# --------------------------------

cls <- read.csv(paste0(OUT_DIR,"rpmm_by_cohort.csv"))

if(COHORT == "TCGA") {
  load(paste0(OUT_DIR, "tcgagbm_dnam.RData"))
} else if (COHORT == "DKFZ") {
  load(paste0(OUT_DIR, "dkfzgbm_dnam.RData"))
  
}
cls <- merge(patients, cls, by="Accession")
rm(gbm450, gbmGENE, patients)

## Administrative censoring:
cls$os_months[cls$os_months > MO_CENSOR] <- MO_CENSOR
cls$os_status[cls$os_months > MO_CENSOR] <- 0 #not dead

## Kaplan-Meier: 
os_univar <- survfit(Surv(os_months, os_status)~Cluster, data=cls)
os_univar #w/ median os & 95%CI by cluster
median(os_univar) #more refined
surv_pvalue(os_univar) #log-rank test
ggsurvplot(os_univar, cls, palette=SURV_COLS, ggtheme=THEME_SURV, pval=FALSE, risk.table=FALSE)

## Cox Model:
os_multi <- coxph(Surv(os_months, os_status)~Cluster+age+sexF, data=cls)
summary(os_multi)
exp(cbind(HR=coef(os_multi), confint(os_multi))) #HR
