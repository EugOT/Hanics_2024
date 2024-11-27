#!/usr/bin/env Rscript
# code/analysis/ref_dev_cortex_umap_optimisation.R

Sys.setenv(RETICULATE_PYTHON = "/opt/python/3.12/bin/python")
# Load tidyverse infrastructure packages
suppressPackageStartupMessages({
  library(here)
  library(tidyverse)
  library(magrittr)
  library(stringr)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(future)
  library(doFuture)
  library(parallelly)
  library(furrr)
  library(future.apply)
  library(reticulate)
})
reticulate::use_condaenv("/opt/python/3.12/bin/python")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
})

options(future.globals.maxSize = 2e11)
options(Seurat.object.assay.version = "v5")

# Set paths
project <- "hanics2024-cortex-tdtomato"
library(here)
library(workflowr)
source(here("code/preprocessing/constants.R"))
source(here(analysis_src, "scDEED.R"))


# set seed
reseed <- 42
set.seed(seed = reseed)
plan(multicore)


# Create a vector with the stage of development for each object
stage_info <- c("E11.5", "E12.5", "E13.5", "E14.5", "E15.5", "E16", "E18.5", "E18", "P1", "P1", "E10", "E17.5", "P4")

orig_umap <- readr::read_tsv(
  here("data/SCP1290/cluster/cluster_scDevSC.merged.umap.txt"),
  skip = 2,
  col_names = c("cell_name", "UMAP_1", "UMAP_2"),
  col_types = list(col_character(), col_double(), col_double())
)

orig_umap %<>% tibble::column_to_rownames("cell_name")
orig_umap %<>% as.matrix()
orig_tsne <- readr::read_tsv(
  here("data/SCP1290/cluster/cluster_scDevSC.merged.tsne.txt"),
  skip = 2,
  col_names = c("cell_name", "tSNE_1", "tSNE_2"),
  col_types = list(col_character(), col_double(), col_double())
)
orig_tsne %<>% tibble::column_to_rownames("cell_name")
orig_tsne %<>% as.matrix()
orig_metadata <- readr::read_tsv(here(
  "data/SCP1290/metadata/metaData_scDevSC.txt"
))
orig_metadata %<>% dplyr::rename("cell_name" = "NAME")
orig_metadata_types <- orig_metadata[1, ] |> purrr::simplify()
orig_metadata %<>% dplyr::filter(!cell_name == "TYPE")
glimpse(orig_metadata)

change_column_types <- function(df, types) {
  for (col_name in names(types)) {
    col_type <- types[col_name]

    if (col_type == "character") {
      df[[col_name]] <- as.character(df[[col_name]])
    } else if (col_type == "numeric") {
      df[[col_name]] <- as.numeric(df[[col_name]])
    } else if (col_type == "integer") {
      df[[col_name]] <- as.integer(df[[col_name]])
    } else if (col_type == "logical") {
      df[[col_name]] <- as.logical(df[[col_name]])
    } else if (col_type == "factor") {
      df[[col_name]] <- as.factor(df[[col_name]])
    } else if (col_type == "group") {
      df[[col_name]] <- as.factor(df[[col_name]])
    } else {
      warning(paste("Unknown type:", col_type, "for column", col_name))
    }
  }

  return(df)
}

# Apply the function to the metadata
orig_metadata <- change_column_types(orig_metadata, orig_metadata_types)

# Print the modified metadata
glimpse(orig_metadata)

orig_srt <- Read10X(data.dir = here("data/SCP1290/expression/601ae2f4771a5b0d72588bfb"))

# Convert the log1p normalized matrix to a standard matrix if it's not already
normalized_matrix <- as.matrix(orig_srt)

# Reverse the log1p transformation to get the scaled count matrix
count_matrix <- expm1(normalized_matrix)

# Extract scaling factors
scaling_factors <- orig_metadata[orig_metadata$cell_name == colnames(count_matrix), ]$nCount_RNA / 1e4

# Multiply each column by its scaling factor and round the results (it's not necessary but just to be sure)
scaled_count_matrix <- sweep(count_matrix, 2, scaling_factors, FUN = "*")
scaled_count_matrix <- round(scaled_count_matrix)

# Convert the count matrix to a sparse matrix format (dgCMatrix) as needed
count_matrix_sparse <- as(scaled_count_matrix, "dgCMatrix")

# Create a Seurat object using the recovered count matrix
merged_cortex <- CreateSeuratObject(counts = count_matrix_sparse, meta.data = orig_metadata)

merged_cortex$stage <- merged_cortex$orig.ident

rm(
  count_matrix,
  count_matrix_sparse,
  normalized_matrix,
  scaled_count_matrix,
  orig_metadata,
  orig_srt
)

table(merged_cortex$New_cellType)
Idents(merged_cortex) <- "New_cellType"
merged_cortex <- subset(merged_cortex, idents = c("Doublet", "Low quality cells", "Red blood cells"), invert = TRUE)

# Set cell type annotations as identities if available
Idents(merged_cortex) <- "orig.ident"

merged_cortex <- FindVariableFeatures(merged_cortex, nfeatures = 5000, verbose = FALSE)
merged_cortex <- NormalizeData(
  merged_cortex,
  verbose = FALSE
)
merged_cortex <- subset(merged_cortex, downsample = 1000)

# Scale data
merged_cortex <- ScaleData(
  merged_cortex,
  verbose = FALSE
)

# Run PCA
merged_cortex <- RunPCA(merged_cortex, npcs = 20, verbose = FALSE)

invisible(gc())
set.seed(reseed)
if (!file.exists(here(data_dir, glue::glue("{project}-init/{project}-init-umap-search-ref.Rds")))) {

  permuted.cortex <- Permuted(merged_cortex, K = 20)

  invisible(gc())
  set.seed(reseed)

  umap_example <- scDEED(
    input_data = merged_cortex,
    K = 20,
    n_neighbors = seq(from = 15, to = 55, by = 10),
    min.dist = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.8),
    reduction.method = "umap",
    rerun = FALSE,
    permuted = permuted.cortex,
    default_assay = "RNA"
  )

  dir.create(here(data_dir, sprintf("%s-init", project)))
  readr::write_rds(
    x = umap_example,
    file = here(data_dir, glue::glue("{project}-init/{project}-init-umap-search-ref.Rds"))
  )
} else {
  umap_example <-
    read_rds(here(data_dir, glue::glue("{project}-init/{project}-init-umap-search-ref.Rds")))
}
