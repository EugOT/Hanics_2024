#!/usr/bin/env Rscript
# vis.R

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

plot_enriched_motif <- function(holder, cluster) {
  src <- c(
    "###### {{cluster}} {.unnumbered}",
    "",
    "```{r w3-control-{{cluster}}-motifs, fig.height=12, fig.width=12}",
    'pluck(holder, as.character(cluster), 1, "motif_plot")',
    "```",
    "",
    "",
    "```{r w3-control-{{cluster}}-coverage, fig.height=26, fig.width=12}",
    'pluck(holder, as.character(cluster), 1, "top_genes") |>',
    "  walk(\\(gene) CoveragePlot(",
    "    combined_srt,",
    "    region = gene,",
    "    features = gene,",
    '    expression.assay = "RNA3",',
    '    expression.slot = "data",',
    '    group.by = "Diet",',
    '    split.by = "cca_clusters",',
    '    assay = "peaks",',
    "    extend.upstream = 5000,",
    "    extend.downstream = 5000,",
    "    annotation = TRUE,",
    "    peaks = TRUE,",
    "    tile = FALSE,",
    "    links = FALSE,",
    "    heights = c(16, 2, 1),",
    "    widths = c(13, 1)",
    "  ))",
    "```",
    "",
    "",
    "```{r w3-control-{{cluster}}-boxplots, fig.height=8, fig.width=16}",
    'pluck(holder, as.character(cluster), 1, "top_genes") |>',
    "  walk(\\(gene) SCpubr::do_BoxPlot(",
    "    sample = combined_srt,",
    "    feature = gene,",
    '    assay = "RNA3",',
    '    slot = "data",',
    '    group.by = "cca_clusters",',
    "    order = FALSE,",
    '    split.by = "Diet"',
    "  ))",
    "```",
    "",
    "",
    ""
  )
  knitr::knit_expand(text = src)
}

generate_motif_plots <- function(results, cluster, comparison_label) {
  # Generates motif plots, expression plots, and box plots for a given cluster and comparison
  src <- c(
    "###### Cluster {{cluster}} {.unnumbered}",
    "```{r plt-overrepresented-motifs-{{comparison_label}}-{{cluster}}-motifs, fig.height=5, fig.width=7}",
    "sig_motifs <- results[['{{cluster}}']] %>% filter(p_val < 0.05) %>% slice_head(n = 6)",
    "MotifPlot(object = combined_srt, motifs = rownames(sig_motifs), assay = 'peaks')",
    "```",
    "",
    "```{r plt-overrepresented-motifs-{{comparison_label}}-{{cluster}}-expression, fig.height=5, fig.width=7}",
    "FeaturePlot(object = combined_srt, features = rownames(sig_motifs), reduction = 'pacmap.cca', min.cutoff = 'q10', max.cutoff = 'q90', pt.size = 0.1)",
    "```",
    "",
    # Corrected part where comparison_label is split into W3 and control
    "```{r plt-overrepresented-motifs-{{comparison_label}}-{{cluster}}-expression-subset, fig.height=6, fig.width=14}",
    "SCpubr::do_FeaturePlot(sample = combined_srt, features = rownames(sig_motifs), assay = 'RNA3', reduction = 'pacmap.cca', order = TRUE, split.by = 'celltype.stim', idents.keep = c(str_c('{{cluster}}', str_split_1(comparison_label, '-'), sep = '_')), ncol = 2)",
    "```",
    ""
  )

  # Expand the text using knitr::knit_expand
  knitr::knit_expand(text = src)
}

render_motif_plots_for_clusters <- function(results, clusters, comparison_label) {
  # Loop through clusters and generate the plots for each
  src_list <- clusters |> map(\(cluster) generate_motif_plots(results, cluster, comparison_label))

  # Render the child document using knit_child
  out <- knitr::knit_child(text = unlist(src_list), options = list(cache = FALSE))
  return(out)
}
