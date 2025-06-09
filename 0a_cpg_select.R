# CpG Selection by Filtering & Intersection

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("utils.R")

XREACTIVE <- fread(paste0(DIR,"annotation_files/infinium_filtering/450K/48639-non-specific-probes-Illumina450k.csv"), header=TRUE, data.table=FALSE)

ANNOT <- load_450k_annot()
ANNOT <- subset(ANNOT, chr %in% paste0("chr",1:22))
ANNOT <- subset(ANNOT, ! Name %in% XREACTIVE$TargetID)

PROP_KEEP <- 0.9 #variance filter

# ------------------ Part I: Shared Autosomal Loci between Datasets ------------------
tcga450 <- load_450k_from_csv(paste0(TCGA_DIR,"GBM/jhu-usc.edu_GBM_HumanMethylation450.betaValue.tsv"))
tcga450 <- tcga450[ , substr(colnames(tcga450), 14, 15) %in% TYPES_INCLUDE]
tcga450 <- select_most_var(tcga450, PROP_KEEP) #variance filter

dkfz450 <- load_450k_from_csv(paste0(DKFZ_DIR,"GSE103659.txt")) #all tumors
dkfz450 <- select_most_var(dkfz450, PROP_KEEP) #variance filter

norm450 <- load_450k_from_csv(paste0(HEALTHY_DIR, "GSE61431_dasen.csv.gz")) #all non-tumors
# No filter intended for healthy dataset

## TODO: Keep OpenSea probes only:
cpg_shared <- Reduce(intersect, list(rownames(tcga450), rownames(dkfz450), rownames(norm450), 
                                     ANNOT$Name))
length(cpg_shared)
rm(tcga450, dkfz450, norm450)

# write.table(cpg_shared, paste0(OUT_DIR,"CpGs_all_shared.csv"), row.names=TRUE, col.names=FALSE, quote=FALSE)

# ------------------ Part II: Locus of Interest from Annotation ------------------
cpg_list <- subset(ANNOT, grepl(GENE, UCSC_RefGene_Name))
cpg_list <- cpg_list[order(cpg_list$pos, decreasing=TRUE), ] #gene on rev strand
# write.csv(cpg_list, paste0(OUT_DIR,"CpGs_GENE.csv"), row.names=FALSE, quote=FALSE)
