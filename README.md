# Epigenetic State of a Fatty Acid Synthesis (FAS) Machinery in Glioblastoma Multiforme (GBM)

## Introduction

Altered lipid metabolism is an emerging hallmark of cancer. In GBM, lipid metabolism is essential to cancer cell proliferation. Using two independent population-scale GM cohorts, this study stratifies each into _predominant subgroups R and L_ by applying unsupervised machine learning to the epigenetic profile of a high-profile candidate gene in FAS. 

Subsequently, this study compared the FAS epigenetic subgroups, one cohort at a time, and...

* Shows significant differences in _cancer stemness_
* Shows significant differences in _cellular differentiation potential_ 
* Identifies distinct _neuronal cell type proportions_
* Demonstrates a significant overlap between subgrouping and global epigenetic landscape 
* Uses EWAS to reveal biologically relevant gene sets and processes 
* Reveals differences in _overall survival_ by subgroups

## Repository Structure

The scripts prefixed with numbers, e.g. `01_`, are used for the analyses presented in the final write-up. File `utils.R` is loaded as a utility module with shared functions and global variables per R session. Some variables in this script are masked.

Scripts for data cleaning, wrangling, exploratory data analyses (EDAs), analyses not presented in the final write-up, and additional data/results I/O are not included.

## Audience Engagement

This work has been submitted for publication and under review. A preprint version of this work is available at [Preprints](https://www.preprints.org/manuscript/202509.0713).

An interactive dashboard of crucial measures, analogous to "KPIs" in business settings, can be found in this [Tableau Public dashboard](https://public.tableau.com/app/profile/david.chen8785/viz/efas/Final).  

Use the data/results, findings, and interpretations of this work responsibly and at your own risk. Cite relevant sources wherevera applicable.

