# CpG Selection by Intersection & Filtering

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

XREACTIVE <- fread(paste0(DIR,"annotation_files/infinium_filtering/450K/48639-non-specific-probes-Illumina450k.csv"), header=TRUE, data.table=FALSE)

ANNOT <- load_450k_annot()
ANNOT <- subset(ANNOT, chr %in% paste0("chr",1:22))
ANNOT <- subset(ANNOT, ! Name %in% XREACTIVE$TargetID)

## Shared autosomal CpGs across all datasets:
tcga450 <- load_450k_from_csv(paste0(TCGA_DIR,"GBM/jhu-usc.edu_GBM_HumanMethylation450.betaValue.tsv"))
dkfz450 <- load_450k_from_csv(paste0(DKFZ_DIR,"GSE103659.txt"))

cpg_shared <- Reduce(intersect, list(rownames(tcga450), rownames(dkfz450), ANNOT$Name))
length(cpg_shared)
rm(tcga450, dkfz450)

## Set aside CpGs from gene of interest:
ANNOT <- subset(ANNOT, Name %in% cpg_shared)

cpg_list <- subset(ANNOT, grepl(GENE, UCSC_RefGene_Name))
cpg_list <- cpg_list[order(cpg_list$pos, decreasing=TRUE), ] #gene on rev strand
nrow(cpg_list)

## Exclude CpGs in gene of interest & MGMT:
mgmt_sites <- ANNOT$Name[grepl("MGMT", ANNOT$UCSC_RefGene_Name)]
cpg_shared <- setdiff(cpg_shared, mgmt_sites)
cpg_shared <- setdiff(cpg_shared, cpg_list$Name)

## Export:
# write.table(cpg_shared, paste0(OUT_DIR,"CpGs_all_shared.csv"), row.names=FALSE, col.names=FALSE, quote=FALSE)
# write.csv(cpg_list, paste0(OUT_DIR,"CpGs_GENE.csv"), row.names=FALSE, quote=FALSE)
