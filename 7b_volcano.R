# Volcano Plots for Differential Methylation

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")
library(ggrepel)

THEME_VOLCANO <- theme_classic() +
  theme(axis.text.x=element_text(color="black",size=16),axis.title.x=element_text(size=21, color="black"), 
        axis.text.y=element_text(color="black",size=16),axis.title.y=element_text(size=21, color="black"),
        legend.title=element_blank(), legend.text=element_blank(), legend.position="none")

ewas_tcga <- read.csv(paste0(OUT_DIR,"ewas_tcga.csv"), row.names=1)
ewas_dkfz <- read.csv(paste0(OUT_DIR,"ewas_dkfz.csv"), row.names=1)

ADJ_P <- 0.05 / nrow(ewas_tcga)
ADJ_P 

wrapper_volcano <- function(df_ewas, xThreshDown=-2.5, xThreshUp=1, yThresh=2.4*-log10(ADJ_P),
                            textHeight=3, pThresh=ADJ_P, esThresh=0) {
  #'@param df_ewas Requires `adj.P.Val` column to be Bonferroni corrected!
  
  ## Operate on transformed Raw P only:
  df_ewas$adj.P.Val <- NULL 
  
  ## Use Raw P value for visualization:
  df_ewas$negLog10P <- -log10(df_ewas$P.Value)
  
  ## Label for direction of change:
  df_ewas$dir[df_ewas$P.Value < pThresh & df_ewas$logFC > esThresh] <- "Positive"
  df_ewas$dir[df_ewas$P.Value < pThresh & df_ewas$logFC < -esThresh] <- "Negative"
  df_ewas$dir[is.na(df_ewas$dir)] <- ""
  
  ## Add per-group N:
  num_pos <- sum(df_ewas$dir=="Positive", na.rm=TRUE)
  num_neg <- sum(df_ewas$dir=="Negative", na.rm=TRUE)
  
  ## Label CpGs w/ largest effect sizes:
  dmp_genes <- as.character(df_ewas$UCSC_RefGene_Name)
  dmp_genes <- strsplit(dmp_genes, split=";")
  df_ewas$Gene <- rep(NA, length(dmp_genes))
  for(k in 1:length(dmp_genes)) {
    dmp_genes[[k]] <- unique(dmp_genes[[k]])
    if(length(dmp_genes[[k]]) > 1) {
      df_ewas$Gene[k] <- paste(dmp_genes[[k]], collapse=";")
    } else if(length(dmp_genes[[k]]) == 1) {
      df_ewas$Gene[k] <- dmp_genes[[k]]
    } else if(length(dmp_genes[[k]]) == 0) {
      df_ewas$Gene[k] <- NA
    }
  }
  
  df_ewas$Label <- df_ewas$Gene
  
  ## Remove some gene-name texts using higher thresholds:
  df_ewas$Label[df_ewas$P.Value > pThresh] <- NA
  df_ewas$Label[df_ewas$logFC > xThreshDown & df_ewas$logFC < xThreshUp] <- NA 
  df_ewas$Label[df_ewas$negLog10P > yThresh] <- df_ewas$Gene[df_ewas$negLog10P > yThresh] #recover most signif
  
  ## Auto-determine text position & axis range/ticks:
  hMax <- max(df_ewas$negLog10P)
  lMax <- min(df_ewas$logFC)
  rMax <- max(df_ewas$logFC)
  
  plt <- ggplot(df_ewas, aes(x=logFC, y=negLog10P, color=dir)) +
    geom_point(aes(size=dir, alpha=dir)) + #override
    scale_color_manual(values = c("gray","royalblue","red")) +
    scale_size_manual(values=c(0.75, 1, 1)) +
    scale_alpha_manual(values=c(0.75, 1, 1)) +
    geom_text_repel(aes(label=Label), color="black", size=4) +
    geom_hline(yintercept=-log10(ADJ_P), linetype="dashed") +
    labs(x="Log2 Fold Diff.", y="-log10(P-value)") +
    annotate("text", rMax/2, hMax+0.5,  label=paste(num_pos,"HYPER \n in Cluster R"), size=6, color="red") +
    annotate("text", lMax/2, hMax+0.5, label=paste(num_neg,"HYPO \n in Cluster R"), size=6, color="royalblue") +
    THEME_VOLCANO
  
  if(esThresh != 0) plt <- plt + geom_vline(xintercept=c(-esThresh, esThresh), linetype="dashed")
  return(plt)
}

wrapper_volcano(ewas_tcga, -2.50)
wrapper_volcano(ewas_dkfz, -2.15)
