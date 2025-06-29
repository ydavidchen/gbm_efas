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

wrapper_volcano <- function(df_ewas, xThreshUp=0, xThreshDown=-2.5, textHeight=3, pThresh=0.05, esThresh=0) {
  #'@param df_ewas Requires `adj.P.Val` column to be Bonferroni corrected!
  
  df_ewas$negLog10P <- -log10(df_ewas$adj.P.Val)
  
  ## Label for direction of change:
  df_ewas$dir[df_ewas$adj.P.Val < pThresh & df_ewas$logFC > esThresh] <- "Positive"
  df_ewas$dir[df_ewas$adj.P.Val < pThresh & df_ewas$logFC < -esThresh] <- "Negative"
  df_ewas$dir[is.na(df_ewas$dir)] <- ""
  
  ## Add per-group N:
  nPos <- sum(df_ewas$dir=="Positive", na.rm=TRUE)
  nNeg <- sum(df_ewas$dir=="Negative", na.rm=TRUE)
  
  ## Label CpGs w/ large effect sizes:
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
  
  ## Remove gene-name texts below threshold:
  df_ewas$Label[(df_ewas$logFC <= xThreshUp & df_ewas$logFC >= xThreshDown) | df_ewas$adj.P.Val >= pThresh] <- NA #higher threshold for gene labeling

  ## Auto-determine text position & axis range/ticks:
  yHoriz <- -log10(pThresh / nrow(df_ewas)) #Undo Bonferroni & then transform
  hMax <- max(df_ewas$negLog10P)
  lMax <- min(df_ewas$logFC)
  rMax <- max(df_ewas$logFC)
  
  plt <- ggplot(df_ewas, aes(x=logFC, y=-log10(P.Value), color=dir)) +
    geom_point(aes(size=dir, alpha=dir)) + #override
    scale_color_manual(values = c("gray","royalblue","red")) +
    scale_size_manual(values=c(1, 1.1, 1.1)) +
    scale_alpha_manual(values=c(0.9, 1, 1)) +
    geom_text_repel(aes(label=Label), color="black", size=4) +
    geom_hline(yintercept=yHoriz, linetype="dashed") +
    labs(x="Log2 Fold Diff.", y="-log10(raw P)") +
    annotate("text", rMax/2, hMax+0.5,  label=paste(nPos,"HYPER \n in Cluster R"), size=6, color="red") +
    annotate("text", lMax/2, hMax+0.5, label=paste(nNeg,"HYPO \n in Cluster R"), size=6, color="royalblue") +
    THEME_VOLCANO
  
  if(esThresh == 0) 
    plt <- plt + geom_vline(xintercept=c(esThresh), linetype="dashed")
  else 
    plt <- plt + geom_vline(xintercept=c(-esThresh, esThresh), linetype="dashed")
  return(plt)
}

wrapper_volcano(ewas_tcga)
wrapper_volcano(ewas_dkfz)
