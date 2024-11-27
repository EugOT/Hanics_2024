#!/usr/bin/env Rscript
# code/analysis/vis.R

set.seed(seed = reseed)

# ggplot2 theme
theme_set(ggmin::theme_powerpoint())

# Function to create feature plots
create_feature_plot <- function(seurat_object, gene_list, file_name, embedding = "wnn.umap", single_pdf = TRUE) {
  # Prepare gene list
  gene_list <- prepare_gene_list(gene_list, seurat_object)

  Iterate_FeaturePlot_scCustom(
    seurat_object = seurat_object,
    reduction = embedding,
    gene_list = gene_list,
    single_pdf = single_pdf,
    raster = length(gene_list) > 12,
    colors_use = viridis(
      n = 30,
      alpha = .55,
      direction = -1,
      option = "E"
    ),
    pt.size = 3,
    alpha_exp = 0.45,
    alpha_na_exp = 0.1,
    file_path = here(plots_dir),
    file_name = str_c(file_name, embedding, sep = "."),
    file_type = ".pdf"
  )
}

# Function to prepare gene list
prepare_gene_list <- function(genes, seurat_object = combined_srt) {
  genes[genes %in% (GetAssayData(seurat_object, slot = "data") %>% as.data.frame() %>% rowSums() %>% .[. > 1] %>% names())]
}

