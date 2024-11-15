#!/usr/bin/env Rscript
# analysis.R

preprocess_raw_srt <- function(project, combined_srt) {
  DefaultAssay(combined_srt) <- "RNA"
  combined_srt$Run <- combined_srt$Run.y
  Idents(combined_srt) <- "Run"

  sex_genes <-
    sex_genes %>% .[. %in% rownames(combined_srt)]
  stress_genes <-
    stress_genes %>% .[. %in% rownames(combined_srt)]

  combined_srt <-
    Store_Palette_Seurat(
      seurat_object = combined_srt,
      palette = rev(brewer.pal(n = 11, name = "RdYlGn")),
      palette_name = "mdat_Colour_Pal"
    )
  combined_srt <-
    Store_Palette_Seurat(
      seurat_object = combined_srt,
      palette = rev(brewer.pal(n = 11, name = "Spectral")),
      palette_name = "expr_Colour_Pal"
    )
  combined_srt <-
    Store_Palette_Seurat(
      seurat_object = combined_srt,
      palette = qc_palette,
      palette_name = "qc_Colour_Pal"
    )


  combined_srt <-
    Add_Mito_Ribo_Seurat(combined_srt, species = "mouse")
  combined_srt[["percent_hb"]] <-
    PercentageFeatureSet(combined_srt, pattern = "^Hb[^(p)]")
  combined_srt <-
    Add_Cell_Complexity_Seurat(combined_srt)
  combined_srt[["var_regex"]] <-
    PercentageFeatureSet(combined_srt, pattern = var_regex)

  DefaultAssay(combined_srt) <- "ATAC"
  combined_srt <- NucleosomeSignal(combined_srt)
  combined_srt <- TSSEnrichment(combined_srt)
  DefaultAssay(combined_srt) <- "RNA"

  combined_srt$QC <-
    ifelse(test = combined_srt$hto_demux == "Doublet",
      yes = "Doublet",
      no = "Pass"
    )

  combined_srt$QC <- ifelse(
    test = combined_srt$hto_demux == "Negative",
    yes = "Negative",
    no = combined_srt$QC
  )

  # Apply QC thresholds to derive categories
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$log10GenesPerUMI < high_cutoff_complexity &
        combined_srt@meta.data$QC == "Pass",
      "Low_Complexity",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$log10GenesPerUMI < high_cutoff_complexity &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "Low_Complexity",
      paste("Low_Complexity", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nFeature_RNA < low_cutoff_gene &
        combined_srt@meta.data$QC == "Pass",
      "Low_nFeature",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nFeature_RNA < low_cutoff_gene &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "Low_nFeature",
      paste("Low_nFeature", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$percent.mt > high_cutoff_pc_mt &
        combined_srt@meta.data$QC == "Pass",
      "High_MT",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$percent.mt > high_cutoff_pc_mt &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "High_MT",
      paste("High_MT", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nCount_RNA > high_cutoff_umis &
        combined_srt@meta.data$QC == "Pass",
      "High_UMIs",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nCount_RNA > high_cutoff_umis &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "High_UMIs",
      paste("High_UMIs", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nCount_ATAC < low_cutoff_fragments &
        combined_srt@meta.data$QC == "Pass",
      "Low_nCount_ATAC",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nCount_ATAC < low_cutoff_fragments &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "Low_nCount_ATAC",
      paste("Low_nCount_ATAC", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nCount_ATAC > high_cutoff_fragments &
        combined_srt@meta.data$QC == "Pass",
      "High_nCount_ATAC",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nCount_ATAC > high_cutoff_fragments &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "High_nCount_ATAC",
      paste("High_nCount_ATAC", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$TSS.enrichment < low_cutoff_TSS_enrichment &
        combined_srt@meta.data$QC == "Pass",
      "Low_TSS_enrichment",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$TSS.enrichment < low_cutoff_TSS_enrichment &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "Low_TSS_enrichment",
      paste("Low_TSS_enrichment", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nucleosome_signal > high_cutoff_nucleosome_signal &
        combined_srt@meta.data$QC == "Pass",
      "High_nucleosome_signal",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$nucleosome_signal > high_cutoff_nucleosome_signal &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "High_nucleosome_signal",
      paste("High_nucleosome_signal", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$percent_ribo > high_cutoff_pc_ribo &
        combined_srt@meta.data$QC == "Pass",
      "High_Ribo",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$percent_ribo > high_cutoff_pc_ribo &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "High_Ribo",
      paste("High_Ribo", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$percent_hb > high_cutoff_pc_hb &
        combined_srt@meta.data$QC == "Pass",
      "High_Hgb",
      combined_srt@meta.data$QC
    )
  combined_srt$QC <-
    ifelse(
      combined_srt@meta.data$percent_hb > high_cutoff_pc_hb &
        combined_srt@meta.data$QC != "Pass" &
        combined_srt@meta.data$QC != "High_Hgb",
      paste("High_Hgb", combined_srt@meta.data$QC, sep = ","),
      combined_srt@meta.data$QC
    )

  return(combined_srt)
}

preprocess_init_filt_srt <- function(project, combined_srt) {
  combined_srt <- NormalizeData(combined_srt)
  combined_srt <-
    FindVariableFeatures(combined_srt,
      selection.method = "vst",
      nfeatures = 5000
    )
  all_genes <- rownames(combined_srt)
  hvg <- VariableFeatures(combined_srt)
  hvg <-
    hvg[str_detect(
      pattern = var_regex,
      string = hvg,
      negate = TRUE
    )]
  agg_genes <-
    GetAssayData(JoinLayers(combined_srt[["RNA"]]), layer = "counts", assay = "RNA") |> rowSums()
  all_genes <-
    all_genes[agg_genes > 0.001 * ncol(combined_srt)]
  all_genes <-
    all_genes[str_detect(
      pattern = var_regex,
      string = all_genes,
      negate = TRUE
    )]

  keep_genes <-
    c(gene_int, hvg) %>%
    unique() %>%
    .[!. %in% housekeeping_mouse] %>%
    .[!. %in% sex_genes] %>%
    .[!. %in% stress_genes]
  glimpse(keep_genes)

  out_of_hvg <- keep_genes[!keep_genes %in% hvg]
  kable_material(
    kable(out_of_hvg, "html"),
    bootstrap_options = c(
      "bordered",
      "condensed",
      "responsive",
      "striped"
    ),
    position = "left",
    font_size = 14
  )

  hvg <- hvg[hvg %in% keep_genes]

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)
  combined_srt[["RNA"]] <- JoinLayers(combined_srt[["RNA"]])
  combined_srt <- ScaleData(
    combined_srt,
    features = hvg,
    vars.to.regress = c("log10GenesPerUMI")
  )

  npcs <- 50
  combined_srt <- RunPCA(
    combined_srt,
    features = hvg,
    npcs = npcs,
    seed.use = reseed,
    verbose = TRUE
  )

  common_row_names <-
    rownames(GetAssayData(JoinLayers(combined_srt[["RNA"]]), layer = "scale.data"))

  c(
    npr,
    np,
    gene_int,
    irs_genes,
    neurotrans,
    glut,
    gaba,
    dopam,
    ach,
    mcr_genes,
    genes.embed,
    genes.manual
  ) %<-% furrr::future_map(
    list(
      npr,
      np,
      gene_int,
      irs_genes,
      neurotrans,
      glut,
      gaba,
      dopam,
      ach,
      mcr_genes,
      genes.embed,
      genes.manual
    ),
    ~ .x[.x %in% common_row_names]
  )


  selected_pcs <-
    seq_along(combined_srt[["pca"]]@stdev)[
      combined_srt[["pca"]]@stdev >
        quantile(combined_srt[["pca"]]@stdev, .10)
    ]

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)
  combined_srt[["RNA"]]$scale.data <- as(object = combined_srt[["RNA"]]$scale.data, Class = "dgCMatrix")

  Idents(combined_srt) <- combined_srt$Run.x
  if (!file.exists(here(data_dir, glue::glue("{project}-init/{project}-init-umap-search.Rds")))) {
    umap_example <- scDEED_simplified(
      input_data = combined_srt,
      num_pc = length(selected_pcs),
      n_neighbors = seq(from = 5, to = 55, by = 10),
      min.dist = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.8),
      use_method = "umap",
      visualization = FALSE
    )

    dir.create(here(data_dir, sprintf("%s-init", project)))
    readr::write_rds(
      x = umap_example,
      file = here(data_dir, glue::glue("{project}-init/{project}-init-umap-search.Rds"))
    )
  } else {
    umap_example <-
      read_rds(here(data_dir, glue::glue("{project}-init/{project}-init-umap-search.Rds")))
  }


  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)


  if (!file.exists(here(data_dir, glue::glue("{project}-init/{project}-init-tsne-search.Rds")))) {
    tsne_example <- scDEED_simplified(
      combined_srt,
      num_pc = length(selected_pcs),
      use_method = "tsne",
      visualization = FALSE
    )

    dir.create(here(data_dir, sprintf("%s-init", project)))
    readr::write_rds(
      x = tsne_example,
      file = here(data_dir, glue::glue("{project}-init/{project}-init-tsne-search.Rds"))
    )
  } else {
    tsne_example <-
      read_rds(here(data_dir, glue::glue("{project}-init/{project}-init-tsne-search.Rds")))
  }

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  Idents(combined_srt) <- combined_srt$Run

  combined_srt <-
    combined_srt |>
    FindNeighbors(
      dims = selected_pcs,
      k.param = umap_example$`best pair of n.neighbors and min.dist` |> pull(n.neighbors),
      annoy.metric = "euclidean",
      n.trees = 100,
      verbose = FALSE
    ) |>
    RunUMAP(
      dims = selected_pcs,
      reduction.name = "umap.rna",
      reduction.key = "rnaUMAP_",
      return.model = FALSE,
      umap.method = "uwot",
      n.epochs = 1000L,
      n.neighbors = umap_example$`best pair of n.neighbors and min.dist` |> pull(n.neighbors),
      min.dist = umap_example$`best pair of n.neighbors and min.dist` |> pull(min.dist),
      seed.use = reseed,
      verbose = FALSE
    )

  combined_srt <-
    RunTSNE(
      combined_srt,
      reduction = "pca",
      dims = selected_pcs,
      seed.use = reseed,
      reduction.name = "tsne.rna",
      reduction.key = "rnatSNE_",
      perplexity = tsne_example$`best perplexity`
    )

  pacmap <- reticulate::import("pacmap")
  # Initialize PaCMAP instance
  reducer <- pacmap$PaCMAP(
    n_components = 2L,
    MN_ratio = 0.5,
    FP_ratio = 2.0,
    apply_pca = FALSE
  )

  # Perform dimensionality Reduction
  pacmap_embedding <-
    reducer$fit_transform(Embeddings(combined_srt[["pca"]])[, selected_pcs])

  colnames(pacmap_embedding) <- paste0("rnaPaCMAP_", 1:2)
  rownames(pacmap_embedding) <- colnames(combined_srt)
  # We will now store this as a custom dimensional reduction called 'pacmap'
  combined_srt[["pacmap.rna"]] <-
    CreateDimReducObject(
      embeddings = pacmap_embedding,
      key = "rnaPaCMAP_",
      assay = DefaultAssay(combined_srt)
    )

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  metadata <- combined_srt@meta.data
  rownames(metadata) <- colnames(combined_srt)
  ref_labels <- metadata$QC

  resolutions <-
    modularity_event_sampling(
      A = combined_srt@graphs$RNA_snn,
      n.res = 4,
      gamma.min = 0.5,
      gamma.max = 3.000001
    ) # sample based on the similarity matrix

  # clustering using Suerat
  combined_srt <- combined_srt %>%
    FindClusters(
      algorithm = "leiden",
      partition.type = "ModularityVertexPartition",
      method = "igraph",
      n.iter = -1,
      resolution = resolutions,
      random.seed = reseed,
      verbose = FALSE
    )

  out <- mrtree(
    combined_srt,
    prefix = "RNA_snn_res.",
    n.cores = n_cores,
    consensus = FALSE,
    sample.weighted = TRUE,
    augment.path = FALSE,
    verbose = FALSE
  )
  # weight per sample is encoraged if the classes are imbalanced

  ks_flat <- apply(
    out$labelmat.flat,
    2,
    FUN = function(x) {
      length(unique(x))
    }
  )
  ks_mrtree <- apply(
    out$labelmat.mrtree,
    2,
    FUN = function(x) {
      length(unique(x))
    }
  )
  amri_flat <-
    sapply(seq_len(ncol(out$labelmat.flat)), function(i) {
      AMRI(out$labelmat.flat[, i], ref_labels)$amri
    })
  amri_flat <-
    aggregate(amri_flat, by = list(k = ks_flat), FUN = mean)
  amri_recon <-
    sapply(seq_len(ncol(out$labelmat.mrtree)), function(i) {
      AMRI(out$labelmat.mrtree[, i], ref_labels)$amri
    })

  df <- rbind(
    data.frame(
      k = amri_flat$k,
      amri = amri_flat$x,
      method = "Seurat flat"
    ),
    data.frame(k = ks_mrtree, amri = amri_recon, method = "MRtree")
  )

  stab_out <- stability_plot(out)

  kable_material(
    kable(
      stab_out$df,
      "html"
    ),
    bootstrap_options = c(
      "bordered",
      "condensed",
      "responsive",
      "striped"
    ),
    position = "left",
    font_size = 14
  )

  res_k <- select_resolution(stab_out$df)

  kable_material(
    kable(
      table(out$labelmat.mrtree[, which.min(abs(as.integer(str_remove(
        dimnames(out$labelmat.mrtree)[[2]], "K"
      )) - res_k))]),
      "html"
    ),
    bootstrap_options = c(
      "bordered",
      "condensed",
      "responsive",
      "striped"
    ),
    position = "left",
    font_size = 14
  )

  combined_srt$k_tree <-
    out$labelmat.mrtree[, which.min(abs(as.integer(str_remove(
      dimnames(out$labelmat.mrtree)[[2]], "K"
    )) - res_k))] %>%
    as.numeric() %>%
    as.factor()

  return(combined_srt)
}

preprocess_post_qc_srt <- function(project, combined_srt) {
  combined_srt <- subset(combined_srt, subset = QC == "Pass")
  DefaultAssay(combined_srt) <- "RNA"
  combined_srt$comb_clstr1 <- Idents(combined_srt)
  s_genes <-
    gorth(
      cc.genes.updated.2019$s.genes,
      source_organism = "hsapiens",
      target_organism = "mmusculus"
    )$ortholog_name
  g2m_genes <-
    gorth(
      cc.genes.updated.2019$g2m.genes,
      source_organism = "hsapiens",
      target_organism = "mmusculus"
    )$ortholog_name
  combined_srt <-
    CellCycleScoring(combined_srt,
      s.features = s_genes,
      g2m.features = g2m_genes
    )

  return(combined_srt)
}

preprocess_srt_norm <- function(project, combined_srt, metadata) {
  combined_srt <- subset(combined_srt, cells = metadata$cell_name)
  DefaultAssay(combined_srt) <- "RNA"

  metadata <-
    metadata |>
    dplyr::select(orig.ident:QC, k_tree:ATAC.weight) |>
    as.data.frame()
  rownames(metadata) <- metadata$cell_name
  metadata <- metadata[colnames(combined_srt), ]
  combined_srt@meta.data <- metadata

  combined_srt <-
    Store_Palette_Seurat(
      seurat_object = combined_srt,
      palette = rev(brewer.pal(n = 11, name = "RdYlGn")),
      palette_name = "mdat_Colour_Pal"
    )
  combined_srt <-
    Store_Palette_Seurat(
      seurat_object = combined_srt,
      palette = rev(brewer.pal(n = 11, name = "Spectral")),
      palette_name = "expr_Colour_Pal"
    )
  combined_srt <-
    Store_Palette_Seurat(
      seurat_object = combined_srt,
      palette = qc_palette,
      palette_name = "qc_Colour_Pal"
    )

  Idents(combined_srt) <- "Run"

  # combined_srt <- Azimuth::RunAzimuth(combined_srt, reference = "mousecortexref")

  return(combined_srt)
}

preprocess_sct_srt <- function(project, combined_srt) {
  # normalize and run dimensionality reduction on control dataset
  npcs <- 50
  metadata <- combined_srt@meta.data
  rownames(metadata) <- colnames(combined_srt)
  combined_srt <-
    SCTransform(
      combined_srt,
      vst.flavor = "v2",
      ncells = ncol(combined_srt),
      variable.features.n = 4500,
      vars.to.regress = c(
        "log10GenesPerUMI",
        "S.Score", "G2M.Score"
      ),
      return.only.var.genes = TRUE,
      seed.use = reseed,
      verbose = FALSE
    )
  hvg <- VariableFeatures(combined_srt)
  hvg <-
    hvg[str_detect(
      pattern = var_regex,
      string = hvg,
      negate = TRUE
    )]

  keep_genes <-
    c(gene_int, hvg) %>%
    unique() %>%
    .[!. %in% housekeeping_mouse] %>%
    .[!. %in% sex_genes] %>%
    .[!. %in% stress_genes]
  glimpse(keep_genes)

  out_of_hvg <- keep_genes[!keep_genes %in% hvg]
  kable_material(
    kable(out_of_hvg, "html"),
    bootstrap_options = c(
      "bordered",
      "condensed",
      "responsive",
      "striped"
    ),
    position = "left",
    font_size = 14
  )

  hvg <- hvg[hvg %in% keep_genes]

  combined_srt <- combined_srt %>%
    RunPCA(
      features = hvg,
      npcs = npcs,
      seed.use = reseed,
      verbose = FALSE
    )

  common_row_names <-
    rownames(combined_srt[["SCT"]]$scale.data)

  c(
    npr,
    np,
    gene_int,
    irs_genes,
    neurotrans,
    glut,
    gaba,
    dopam,
    ach,
    mcr_genes,
    genes.embed,
    genes.manual
  ) %<-% furrr::future_map(
    list(
      npr,
      np,
      gene_int,
      irs_genes,
      neurotrans,
      glut,
      gaba,
      dopam,
      ach,
      mcr_genes,
      genes.embed,
      genes.manual
    ),
    ~ .x[.x %in% common_row_names]
  )

  combined_srt[["SCT"]]$scale.data <- as(object = combined_srt[["SCT"]]$scale.data, Class = "matrix")

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)


  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-sct/{project}-sct-umap-search.Rds")
  ))) {
    umap_example <- scDEED_simplified(
      input_data = combined_srt,
      num_pc = npcs,
      n_neighbors = c(
        seq(
          from = 5, to = 10, by = 1
        ),
        seq(
          from = 15, to = 35, by = 5
        ),
        50, 80
      ),
      min.dist = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.8),
      use_method = "umap",
      visualization = FALSE
    )

    dir.create(here(data_dir, sprintf("%s-sct", project)))
    readr::write_rds(
      x = umap_example,
      file = here(
        data_dir,
        glue::glue("{project}-sct/{project}-sct-umap-search.Rds")
      )
    )
  } else {
    umap_example <-
      read_rds(here(
        data_dir,
        glue::glue("{project}-sct/{project}-sct-umap-search.Rds")
      ))
  }


  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)


  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-sct/{project}-sct-tsne-search.Rds")
  ))) {
    tsne_example <- scDEED_simplified(
      input_data = combined_srt,
      num_pc = npcs,
      use_method = "tsne",
      visualization = FALSE
    )

    dir.create(here(data_dir, sprintf("%s-sct", project)))
    readr::write_rds(
      x = tsne_example,
      file = here(
        data_dir,
        glue::glue("{project}-sct/{project}-sct-tsne-search.Rds")
      )
    )
  } else {
    tsne_example <-
      read_rds(here(
        data_dir,
        glue::glue("{project}-sct/{project}-sct-tsne-search.Rds")
      ))
  }

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  combined_srt <-
    combined_srt |>
    FindNeighbors(
      dims = seq_along(combined_srt[["pca"]]@stdev),
      k.param = umap_example$`best pair of n.neighbors and min.dist` |> pull(n.neighbors),
      annoy.metric = "euclidean",
      n.trees = 100,
      verbose = FALSE
    ) |>
    RunUMAP(
      dims = seq_along(combined_srt[["pca"]]@stdev),
      reduction.name = "umap.rna",
      reduction.key = "rnaUMAP_",
      return.model = FALSE,
      umap.method = "uwot",
      n.epochs = 1000L,
      n.neighbors = umap_example$`best pair of n.neighbors and min.dist` |> pull(n.neighbors),
      min.dist = umap_example$`best pair of n.neighbors and min.dist` |> pull(min.dist),
      seed.use = reseed,
      verbose = FALSE
    )

  combined_srt <-
    RunTSNE(
      combined_srt,
      reduction = "pca",
      dims = seq_along(combined_srt[["pca"]]@stdev),
      seed.use = reseed,
      reduction.name = "tsne.rna",
      reduction.key = "rnatSNE_",
      perplexity = tsne_example$`best perplexity`
    )

  pacmap <- reticulate::import("pacmap")
  # Initialize PaCMAP instance
  reducer <- pacmap$PaCMAP(
    n_components = 2L,
    MN_ratio = 0.5,
    FP_ratio = 2.0,
    apply_pca = FALSE
  )

  # Perform dimensionality Reduction
  pacmap_embedding <-
    reducer$fit_transform(Embeddings(combined_srt[["pca"]])[, seq_along(combined_srt[["pca"]]@stdev)], init = "pca")

  colnames(pacmap_embedding) <- paste0("rnaPaCMAP_", 1:2)
  rownames(pacmap_embedding) <- colnames(combined_srt)
  # We will now store this as a custom dimensional reduction called 'pacmap'
  combined_srt[["pacmap.rna"]] <-
    CreateDimReducObject(
      embeddings = pacmap_embedding,
      key = "rnaPaCMAP_",
      assay = DefaultAssay(combined_srt)
    )

  DefaultAssay(combined_srt) <- "ATAC"
  combined_srt <- RunTFIDF(combined_srt)
  combined_srt <- FindTopFeatures(combined_srt, min.cutoff = "q0")
  combined_srt <- RunSVD(combined_srt)
  combined_srt <-
    RunUMAP(
      combined_srt,
      reduction = "lsi",
      dims = 2:50,
      reduction.name = "umap.atac",
      reduction.key = "atacUMAP_",
      return.model = FALSE,
      umap.method = "uwot",
      n.epochs = 1000L,
      n.neighbors = umap_example$`best pair of n.neighbors and min.dist` |> pull(n.neighbors),
      min.dist = umap_example$`best pair of n.neighbors and min.dist` |> pull(min.dist),
      seed.use = reseed,
      verbose = FALSE
    )

  combined_srt <-
    RunTSNE(
      combined_srt,
      reduction = "lsi",
      dims = 2:50,
      seed.use = reseed,
      reduction.name = "tsne.atac",
      reduction.key = "atactSNE_",
      perplexity = tsne_example$`best perplexity`
    )

  # Perform dimensionality Reduction
  pacmap_embedding <-
    reducer$fit_transform(Embeddings(combined_srt[["lsi"]])[, 2:50], init = "pca")

  colnames(pacmap_embedding) <- paste0("atacPaCMAP_", 1:2)
  rownames(pacmap_embedding) <- colnames(combined_srt)
  # We will now store this as a custom dimensional reduction called 'pacmap'
  combined_srt[["pacmap.atac"]] <-
    CreateDimReducObject(
      embeddings = pacmap_embedding,
      key = "atacPaCMAP_",
      assay = DefaultAssay(combined_srt)
    )



  combined_srt <-
    FindMultiModalNeighbors(
      combined_srt,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:npcs, 2:50),
      k.nn = umap_example$`best pair of n.neighbors and min.dist` |> pull(n.neighbors)
    )
  combined_srt <-
    RunUMAP(
      combined_srt,
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_",
      return.model = FALSE,
      umap.method = "uwot",
      n.epochs = 1000L,
      n.neighbors = umap_example$`best pair of n.neighbors and min.dist` |> pull(n.neighbors),
      min.dist = umap_example$`best pair of n.neighbors and min.dist` |> pull(min.dist),
      seed.use = reseed,
      verbose = FALSE
    )

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  metadata <- combined_srt@meta.data
  rownames(metadata) <- colnames(combined_srt)
  ref_labels <- metadata$k_tree

  resolutions <-
    modularity_event_sampling(
      A = combined_srt@graphs$wsnn,
      n.res = 20,
      gamma.min = 0.8,
      gamma.max = 4.000001
    ) # sample based on the similarity matrix

  # clustering using Suerat
  combined_srt <- combined_srt %>%
    FindClusters(
      algorithm = "leiden",
      partition.type = "ModularityVertexPartition",
      method = "igraph",
      n.iter = -1,
      resolution = resolutions,
      random.seed = reseed,
      verbose = FALSE,
      graph.name = "wsnn"
    )

  out <- mrtree(
    combined_srt,
    prefix = "wsnn_res.",
    n.cores = n_cores,
    consensus = FALSE,
    sample.weighted = TRUE,
    augment.path = FALSE,
    verbose = FALSE
  )
  # weight per sample is encoraged if the classes are imbalanced

  # Adjusted Multiresolution Rand Index (AMRI)
  ks_flat <- apply(
    out$labelmat.flat,
    2,
    FUN = function(x) {
      length(unique(x))
    }
  )
  ks_mrtree <- apply(
    out$labelmat.mrtree,
    2,
    FUN = function(x) {
      length(unique(x))
    }
  )
  amri_flat <-
    sapply(seq_len(ncol(out$labelmat.flat)), function(i) {
      AMRI(out$labelmat.flat[, i], ref_labels)$amri
    })
  amri_flat <-
    aggregate(amri_flat, by = list(k = ks_flat), FUN = mean)
  amri_recon <-
    sapply(seq_len(ncol(out$labelmat.mrtree)), function(i) {
      AMRI(out$labelmat.mrtree[, i], ref_labels)$amri
    })

  df <- rbind(
    data.frame(
      k = amri_flat$k,
      amri = amri_flat$x,
      method = "Seurat flat"
    ),
    data.frame(k = ks_mrtree, amri = amri_recon, method = "MRtree")
  )

  stab_out <- stability_plot(out)

  kable_material(
    kable(
      stab_out$df,
      "html"
    ),
    bootstrap_options = c(
      "bordered",
      "condensed",
      "responsive",
      "striped"
    ),
    position = "left",
    font_size = 14
  )

  res_k <- select_resolution(stab_out$df)

  kable_material(
    kable(
      table(out$labelmat.mrtree[, which.min(abs(as.integer(str_remove(
        dimnames(out$labelmat.mrtree)[[2]], "K"
      )) - res_k))]),
      "html"
    ),
    bootstrap_options = c(
      "bordered",
      "condensed",
      "responsive",
      "striped"
    ),
    position = "left",
    font_size = 14
  )

  combined_srt$k_tree <-
    out$labelmat.mrtree[, which.min(abs(as.integer(str_remove(
      dimnames(out$labelmat.mrtree)[[2]], "K"
    )) - res_k))] %>%
    as.numeric() %>%
    as.factor()

  Idents(combined_srt) <- "k_tree"

  # TODO: first re-call peaks with MACS2 for particular clusters and then
  # do motifs assignment
  # Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
  # probably it should be another function and notebook
  #
  # DefaultAssay(combined_srt) <- "ATAC"
  #
  # # call peaks using MACS2
  # peaks <- CallPeaks(combined_srt, group.by = "k_tree")
  #
  # # remove peaks on nonstandard chromosomes and in genomic blacklist regions
  # peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
  # peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_mm10_unified, invert = TRUE)
  # # quantify counts in each peak
  # macs2_counts <- FeatureMatrix(
  #   fragments = Fragments(combined_srt),
  #   features = peaks,
  #   cells = colnames(combined_srt)
  # )
  # # create a new assay using the MACS2 peak set and add it to the Seurat object
  # combined_srt[["peaks"]] <- CreateChromatinAssay(
  #   counts = macs2_counts,
  #   fragments = fragpath,
  #   annotation = annotation
  # )
  #
  # pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 10090, all_versions = FALSE))
  # motif.matrix <- CreateMotifMatrix(features = granges(combined_srt), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
  # motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
  # combined_srt <- SetAssayData(combined_srt, assay = 'ATAC', slot = 'motifs', new.data = motif.object)
  #
  # # Note that this step can take 30-60 minutes
  # combined_srt <- RunChromVAR(
  #   object = combined_srt,
  #   genome = BSgenome.Mmusculus.UCSC.mm10
  # )

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  combined_srt <-
    PrepSCTFindMarkers(combined_srt, assay = "SCT")

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  markers_logreg <-
    FindAllMarkers(
      combined_srt,
      assay = "SCT",
      verbose = FALSE,
      random.seed = reseed,
      only.pos = TRUE,
      min.pct = 0.2,
      base = 10,
      logfc.threshold = 0.2,
      densify = TRUE,
      test.use = "LR"
    )

  write_csv(
    markers_logreg,
    here(
      tables_dir,
      sprintf(
        "%s_all_mrk-logreg_sct-combined-whole_dataset.csv",
        project
      )
    )
  )

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  markers_logreg %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log10FC) %>%
    kable("html") %>%
    kable_material(
      bootstrap_options = c(
        "bordered",
        "condensed",
        "responsive",
        "striped"
      ),
      position = "left",
      font_size = 14
    )

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  markers_MAST <-
    FindAllMarkers(
      combined_srt,
      assay = "SCT",
      verbose = FALSE,
      random.seed = reseed,
      only.pos = TRUE,
      min.pct = 0.1,
      base = 10,
      logfc.threshold = 0.2,
      test.use = "MAST"
    )
  write_csv(
    markers_MAST,
    here(
      tables_dir,
      sprintf(
        "%s_all_mrk-MAST_sct-combined-whole_dataset.csv",
        project
      )
    )
  )
  markers_MAST %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log10FC) %>%
    kable("html") %>%
    kable_material(
      bootstrap_options = c(
        "bordered",
        "condensed",
        "responsive",
        "striped"
      ),
      position = "left",
      font_size = 14
    )

  return(combined_srt)
}

integrate_srt <- function(project, combined_srt, unintegrated_srt) {
  # normalize and run dimensionality reduction on control dataset
  npcs <- 50
  DefaultAssay(combined_srt) <- "RNA"
  combined_srt[["RNA"]]@layers <-
    combined_srt[["RNA"]]@layers[combined_srt[["RNA"]]@layers |>
                                   names() |>
                                   str_detect(pattern = "Gene Expression")]
  metadata <- combined_srt@meta.data
  rownames(metadata) <- colnames(combined_srt)

  # DefaultAssay(unintegrated_srt) <- "SCT"
  # unintegrated_srt[["SCT"]] <-
  #   split(unintegrated_srt[["SCT"]], f = unintegrated_srt$Run)
  #
  # all_genes <- rownames(unintegrated_srt)
  # hvg <- VariableFeatures(unintegrated_srt)
  # hvg <-
  #   hvg[str_detect(pattern = var_regex,
  #                  string = hvg,
  #                  negate = TRUE)]
  # agg_genes <-
  #   JoinLayers(unintegrated_srt[["SCT"]])$counts |> rowSums()
  # all_genes <-
  #   all_genes[agg_genes > 0.001 * ncol(unintegrated_srt)]
  # all_genes <-
  #   all_genes[str_detect(pattern = var_regex,
  #                        string = all_genes,
  #                        negate = TRUE)]
  #
  # keep_genes <-
  #   c(gene_int, hvg) %>%
  #   unique() %>%
  #   .[. %in% all_genes] %>%
  #   .[!. %in% housekeeping_mouse] %>%
  #   .[!. %in% sex_genes] %>%
  #   .[!. %in% stress_genes]

  # combined_srt <-
  #   SCTransform(
  #     combined_srt,
  #     residual.features = keep_genes,
  #     vars.to.regress = c("log10GenesPerUMI"),
  #     return.only.var.genes = FALSE,
  #     do.scale = TRUE,
  #     do.center = TRUE,
  #     seed.use = reseed
  #     )
  combined_srt <- NormalizeData(combined_srt)
  combined_srt <-
    FindVariableFeatures(combined_srt,
                         selection.method = "vst",
                         nfeatures = 5000
    )

  invisible(gc())
  set.seed(seed = reseed)

  all_genes <- rownames(combined_srt)
  hvg <- VariableFeatures(combined_srt)
  hvg <-
    hvg[str_detect(pattern = var_regex,
                   string = hvg,
                   negate = TRUE)]
  agg_genes <-
    JoinLayers(combined_srt[["RNA"]])$counts |> rowSums()
  all_genes <-
    all_genes[agg_genes > 0.001 * ncol(combined_srt)]
  all_genes <-
    all_genes[str_detect(pattern = var_regex,
                         string = all_genes,
                         negate = TRUE)]

  keep_genes <-
    c(gene_int, hvg) %>%
    unique() %>%
    .[. %in% all_genes] %>%
    .[!. %in% housekeeping_mouse] %>%
    .[!. %in% sex_genes] %>%
    .[!. %in% stress_genes]
  glimpse(keep_genes)

  out_of_hvg <- keep_genes[!keep_genes %in% hvg]

  hvg <- hvg[hvg %in% keep_genes]

  combined_srt <- ScaleData(
    combined_srt,
    features = keep_genes,
    vars.to.regress = c("log10GenesPerUMI")
  )

  combined_srt <- RunPCA(
    combined_srt,
    features = keep_genes,
    npcs = npcs,
    seed.use = reseed,
    verbose = TRUE
  )

  # combined_srt[["SCT"]] <- split(combined_srt[["SCT"]], f = combined_srt$Run)

  # BP matrix doesn't work yet!
  # unintegrated_srt <- Azimuth::RunAzimuth(unintegrated_srt, reference = "mousecortexref")
  # # Warning: Overwriting miscellanous data for model
  # # Warning: Adding a dimensional reduction (refUMAP) without the associated assay being present
  # # Warning: Adding a dimensional reduction (refUMAP) without the associated assay being present
  # # detected inputs from MOUSE with id type Gene.name
  # # reference rownames detected MOUSE with id type Gene.name
  # # Normalizing query using reference SCT model
  # # Warning: 105 features of the features specified were not present in both the reference query assays.
  # # Continuing with remaining 2895 features.
  # # Projecting cell embeddings
  # # Finding query neighbors
  # # Error in as(object = y[[i]], Class = "CsparseMatrix") :
  # #   no method or default for coercing "RenameDims" to "CsparseMatrix"

  combined_srt <- IntegrateLayers(
    object = combined_srt, method = CCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.cca",
    features = keep_genes, normalization.method = "LogNormalize",
    verbose = FALSE
  )
  combined_srt <- IntegrateLayers(
    object = combined_srt, method = RPCAIntegration,
    orig.reduction = "pca", new.reduction = "integrated.rpca",
    features = keep_genes, normalization.method = "LogNormalize",
    verbose = FALSE
  )
  combined_srt <- IntegrateLayers(
    object = combined_srt, method = HarmonyIntegration,
    orig.reduction = "pca", new.reduction = "harmony",
    features = keep_genes, normalization.method = "LogNormalize",
    verbose = FALSE
  )
  combined_srt <- IntegrateLayers(
    object = combined_srt, method = FastMNNIntegration, k = 25,
    new.reduction = "integrated.mnn", batch = combined_srt$Run, verbose = FALSE
  )

  # combined_srt <- IntegrateLayers(
  #   object = combined_srt, method = scVIIntegration,
  #   new.reduction = "integrated.scvi", normalization.method = "LogNormalize",
  #   conda_env = "/opt/python/3.8.8", verbose = T
  # )

  pacmap <- reticulate::import("pacmap")
  # Initialize PaCMAP instance
  reducer <- pacmap$PaCMAP(
    n_components = 2L,
    MN_ratio = 0.5,
    FP_ratio = 2.0,
    apply_pca = FALSE
  )

  # Perform dimensionality Reduction
  pacmap_embedding <-
    reducer$fit_transform(Embeddings(combined_srt[["integrated.cca"]])[, 1:30], init = "pca")

  colnames(pacmap_embedding) <- paste0("ccaPaCMAP_", 1:2)
  rownames(pacmap_embedding) <- colnames(combined_srt)
  # We will now store this as a custom dimensional reduction called 'pacmap'
  combined_srt[["pacmap.cca"]] <-
    CreateDimReducObject(
      embeddings = pacmap_embedding,
      key = "ccaPaCMAP_",
      assay = DefaultAssay(combined_srt)
    )

  # Initialize PaCMAP instance
  reducer <- pacmap$PaCMAP(
    n_components = 2L,
    MN_ratio = 0.5,
    FP_ratio = 2.0,
    apply_pca = FALSE
  )

  # Perform dimensionality Reduction
  pacmap_embedding <-
    reducer$fit_transform(Embeddings(combined_srt[["integrated.rpca"]])[, 1:30], init = "pca")

  colnames(pacmap_embedding) <- paste0("rpcaPaCMAP_", 1:2)
  rownames(pacmap_embedding) <- colnames(combined_srt)
  # We will now store this as a custom dimensional reduction called 'pacmap'
  combined_srt[["pacmap.rpca"]] <-
    CreateDimReducObject(
      embeddings = pacmap_embedding,
      key = "rpcaPaCMAP_",
      assay = DefaultAssay(combined_srt)
    )

  # Initialize PaCMAP instance
  reducer <- pacmap$PaCMAP(
    n_components = 2L,
    MN_ratio = 0.5,
    FP_ratio = 2.0,
    apply_pca = FALSE
  )

  # Perform dimensionality Reduction
  pacmap_embedding <-
    reducer$fit_transform(Embeddings(combined_srt[["harmony"]])[, 1:30], init = "pca")

  colnames(pacmap_embedding) <- paste0("harmonyPaCMAP_", 1:2)
  rownames(pacmap_embedding) <- colnames(combined_srt)
  # We will now store this as a custom dimensional reduction called 'pacmap'
  combined_srt[["pacmap.harmony"]] <-
    CreateDimReducObject(
      embeddings = pacmap_embedding,
      key = "harmonyPaCMAP_",
      assay = DefaultAssay(combined_srt)
    )

  # Initialize PaCMAP instance
  reducer <- pacmap$PaCMAP(
    n_components = 2L,
    MN_ratio = 0.5,
    FP_ratio = 2.0,
    apply_pca = FALSE
  )

  # Perform dimensionality Reduction
  pacmap_embedding <-
    reducer$fit_transform(Embeddings(combined_srt[["integrated.mnn"]])[, 1:30], init = "pca")

  colnames(pacmap_embedding) <- paste0("mnnPaCMAP_", 1:2)
  rownames(pacmap_embedding) <- colnames(combined_srt)
  # We will now store this as a custom dimensional reduction called 'pacmap'
  combined_srt[["pacmap.mnn"]] <-
    CreateDimReducObject(
      embeddings = pacmap_embedding,
      key = "mnnPaCMAP_",
      assay = DefaultAssay(combined_srt)
    )

  # # Initialize PaCMAP instance
  # reducer <- pacmap$PaCMAP(
  #   n_components = 2L,
  #   MN_ratio = 0.5,
  #   FP_ratio = 2.0,
  #   apply_pca = FALSE
  # )
  #
  # # Perform dimensionality Reduction
  # pacmap_embedding <-
  #   reducer$fit_transform(Embeddings(combined_srt[["integrated.scvi"]])[, 1:30], init = "pca")
  #
  # colnames(pacmap_embedding) <- paste0("scviPaCMAP_", 1:2)
  # rownames(pacmap_embedding) <- colnames(combined_srt)
  # # We will now store this as a custom dimensional reduction called 'pacmap'
  # combined_srt[["pacmap.scvi"]] <-
  #   CreateDimReducObject(
  #     embeddings = pacmap_embedding,
  #     key = "scviPaCMAP_",
  #     assay = DefaultAssay(combined_srt)
  #   )

  combined_srt <-
    FindNeighbors(
      combined_srt,
      reduction = "integrated.cca",
      graph.name = "integrated.cca.snn",
      dims = 1:30
    )
  combined_srt <- combined_srt %>%
    FindClusters(
      algorithm = "leiden",
      cluster.name = "cca_clusters",
      partition.type = "ModularityVertexPartition",
      method = "igraph",
      n.iter = -1,
      resolution = 2,
      random.seed = reseed,
      verbose = FALSE,
      graph.name = "integrated.cca.snn"
    )

  combined_srt <-
    FindNeighbors(
      combined_srt,
      reduction = "integrated.rpca",
      graph.name = "integrated.rpca.snn",
      dims = 1:30
    )
  combined_srt <- combined_srt %>%
    FindClusters(
      algorithm = "leiden",
      cluster.name = "rpca_clusters",
      partition.type = "ModularityVertexPartition",
      method = "igraph",
      n.iter = -1,
      resolution = 2,
      random.seed = reseed,
      verbose = FALSE,
      graph.name = "integrated.rpca.snn"
    )

  combined_srt <-
    FindNeighbors(
      combined_srt,
      reduction = "harmony",
      graph.name = "harmony.snn",
      dims = 1:30
    )
  combined_srt <- combined_srt %>%
    FindClusters(
      algorithm = "leiden",
      cluster.name = "harmony_clusters",
      partition.type = "ModularityVertexPartition",
      method = "igraph",
      n.iter = -1,
      resolution = 2,
      random.seed = reseed,
      verbose = FALSE,
      graph.name = "harmony.snn"
    )

  combined_srt <-
    FindNeighbors(
      combined_srt,
      reduction = "integrated.mnn",
      graph.name = "integrated.mnn.snn",
      dims = 1:30
    )
  combined_srt <- combined_srt %>%
    FindClusters(
      algorithm = "leiden",
      cluster.name = "mnn_clusters",
      partition.type = "ModularityVertexPartition",
      method = "igraph",
      n.iter = -1,
      resolution = 2,
      random.seed = reseed,
      verbose = FALSE,
      graph.name = "integrated.mnn.snn"
    )

  # combined_srt <-
  #   FindNeighbors(
  #     combined_srt,
  #     reduction = "integrated.scvi",
  #     graph.name = "integrated.scvi.snn",
  #     dims = 1:30
  #   )
  # combined_srt <- combined_srt %>%
  #   FindClusters(
  #     algorithm = "leiden",
  #     cluster.name = "scvi_clusters",
  #     partition.type = "ModularityVertexPartition",
  #     method = "igraph",
  #     n.iter = -1,
  #     resolution = 2,
  #     random.seed = reseed,
  #     verbose = FALSE,
  #     graph.name = "integrated.scvi.snn"
  #   )

  combined_srt <- JoinLayers(combined_srt)

  return(combined_srt)
}

quantify_de_integrated <- function(clusters) {
  Idents(combined_srt) <- clusters

  markers_logreg <-
    FindAllMarkers(
      combined_srt,
      assay = "RNA",
      verbose = FALSE,
      random.seed = reseed,
      only.pos = TRUE,
      min.pct = 0.2,
      base = 10,
      logfc.threshold = 0.2,
      densify = TRUE,
      test.use = "LR"
    )

  write_csv(
    markers_logreg,
    here(
      tables_dir,
      glue::glue("{project}_all_mrk-logreg-{clusters}-whole_dataset.csv")
    )
  )

  markers_logreg %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log10FC) %>%
    kable("html") %>%
    kable_material(
      bootstrap_options = c(
        "bordered",
        "condensed",
        "responsive",
        "striped"
      ),
      position = "left",
      font_size = 14
    )

  plan(sequential)
  invisible(gc())
  set.seed(seed = reseed)
  plan(multisession, workers = n_cores)

  markers_MAST <-
    FindAllMarkers(
      combined_srt,
      assay = "RNA",
      verbose = FALSE,
      random.seed = reseed,
      only.pos = TRUE,
      min.pct = 0.1,
      base = 10,
      logfc.threshold = 0.2,
      test.use = "MAST"
    )
  write_csv(
    markers_MAST,
    here(
      tables_dir,
      glue::glue("{project}_all_mrk-MAST-{clusters}-whole_dataset.csv")
    )
  )
  markers_MAST %>%
    group_by(cluster) %>%
    slice_max(n = 20, order_by = avg_log10FC) %>%
    kable("html") %>%
    kable_material(
      bootstrap_options = c(
        "bordered",
        "condensed",
        "responsive",
        "striped"
      ),
      position = "left",
      font_size = 14
    )
  return(
    list(
      "LR" = markers_logreg,
      "MAST" = markers_MAST
    )
  )
}
