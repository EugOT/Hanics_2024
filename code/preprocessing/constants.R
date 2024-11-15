#!/usr/bin/env Rscript
# constants.R

Sys.setenv(RETICULATE_PYTHON = "/opt/python/3.8.8/bin/python")
# Load tidyverse infrastructure packages
suppressPackageStartupMessages({
  library(here)
  library(knitr)
  library(RColorBrewer)
  library(viridis)
  library(tidyverse)
  library(magrittr)
  library(stringr)
  library(skimr)
  library(future)
  library(furrr)
  library(zeallot)
  library(kableExtra)
  library(reticulate)
})
reticulate::use_condaenv("/opt/python/3.8.8/bin/python")

# Set paths
src_dir <- here("code")
preprocessing_src <- here(src_dir, "preprocessing")
analysis_src <- here(src_dir, "analysis")
analysis_dir <- here("analysis")
data_dir <- here("data")
output_dir <- here("output")
plots_dir <- here(output_dir, "figures/")
tables_dir <- here(output_dir, "tables/")

# parallelisation
n_cores <- 16

# set seed
reseed <- 42
set.seed(seed = reseed)

# ggplot2 theme
theme_set(ggmin::theme_powerpoint())


cb_fpr <- 0.001
low_cutoff_gene <- 200
high_cutoff_gene <- 2500
low_cutoff_umis <- NULL
low_cutoff_umis <- -Inf
high_cutoff_umis <- 10000
high_cutoff_pc_mt <- 20
high_cutoff_pc_ribo <- 2
high_cutoff_pc_hb <- 0.1
high_cutoff_doublet_score <- 0.33
high_cutoff_complexity <- 0.90

var_regex <- "^Hla-|^Ig[hjkl]|^Rna|^mt-|^Rp[sl]|^Hb[^(p)]|^Gm"

# TODO: update this with additional categories of QC filter!
###### COLOURS: ################################################################
qc_palette <- c(
  "Doublet" = "#5050FFFF",
  "Pass" = "#ff3700f1",
  "High_MT" = "#6BD76BFF",
  "High_Hgb" = "#C75127FF",
  "Low_Complexity" = "#3B1B53FF",
  "High_Ribo" = "#E7C76FFF",
  "Low_nFeature" = "#A9A9A9FF",
  "High_UMIs" = "#1A0099FF",
  "High_Hgb,High_Ribo,Low_nFeature" = "#660099FF",
  "High_UMIs,Doublet" = "#990080FF",
  "High_Hgb,High_Ribo,High_MT,Low_nFeature,Doublet" = "#FF1463FF",
  "High_Ribo,Low_nFeature,Doublet" = "#00D68FFF",
  "High_Hgb,High_Ribo,Low_nFeature,Doublet" = "#D60047FF",
  "High_Ribo,High_MT" = "#00CC33FF",
  "High_Hgb,High_Ribo" = "#4775FFFF",
  "High_Ribo,High_MT,Low_Complexity" = "#991A00FF",
  "High_Hgb,Low_nFeature,Low_Complexity" = "#00991AFF",
  "High_Hgb,Low_nFeature" = "#CE3D32FF",
  "High_Hgb,High_MT" = "#BA6338FF",
  "High_MT,Low_nFeature" = "#D595A7FF",
  "High_Ribo,High_MT,Doublet" = "#D58F5CFF",
  "High_Hgb,High_Ribo,High_MT,Low_nFeature" = "#CDDEB7FF",
  "High_MT,Low_Complexity" = "#5A655EFF",
  "High_Ribo,High_MT,Low_nFeature" = "#CC9900FF",
  "High_MT,Low_nFeature,Low_Complexity" = "#00CC99FF",
  "High_MT,Doublet" = "#FFC20AFF",
  "High_Hgb,High_Ribo,High_MT" = "#996600FF",
  "Low_nFeature,Low_Complexity" = "#009966FF",
  "High_Ribo,Doublet" = "#749B58FF",
  "High_Hgb,Doublet" = "#5DB1DDFF",
  "Low_nFeature,Doublet" = "#924822FF",
  "High_Hgb,High_MT,Low_nFeature" = "#7A65A5FF",
  "High_Hgb,High_Ribo,Low_Complexity" = "#612A79FF",
  "High_MT,Low_nFeature,Doublet" = "#CC9900FF",
  "High_Ribo,High_MT,Low_nFeature,Low_Complexity" = "#99CC00FF",
  "High_UMIs,Low_Complexity" = "#0099CCFF",
  "High_Ribo,High_UMIs,Low_Complexity" = "#FFD147FF",
  "High_UMIs,Low_Complexity,Doublet" = "#809900FF",
  "High_Hgb,High_Ribo,High_MT,Doublet" = "#008099FF",
  "High_Hgb,Low_Complexity" = "#F0E685FF",
  "High_Ribo,Low_nFeature" = "#802268FF",
  "High_UMIs,High_MT,Low_Complexity" = "#837B8DFF",
  "High_Hgb,Low_Complexity,Doublet" = "#E4AF69FF",
  "High_Hgb,High_Ribo,Doublet" = "#AE1F63FF",
  "High_Hgb,High_MT,Low_nFeature,Low_Complexity" = "#99CC00FF",
  "High_Hgb,High_Ribo,High_MT,Low_Complexity,Doublet" = "#33CC00FF",
  "High_Hgb,High_Ribo,High_UMIs,Low_Complexity" = "#0A47FFFF",
  "Low_Complexity,Doublet" = "#990033FF",
  "High_Ribo,Low_Complexity" = "#339900FF",
  "High_Hgb,High_Ribo,Low_nFeature,Low_Complexity" = "#003399FF"
)
