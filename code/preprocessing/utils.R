#!/usr/bin/env Rscript
# utils.R
# This script contains various utility functions for data preprocessing and analysis.

suppressPackageStartupMessages({
  library(BPCells)
  library(Seurat)
  library(scDEED)
  library(Signac)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(EnsDb.Mmusculus.v79)
  library(chromVAR)
  library(JASPAR2020)
  library(TFBSTools)
  library(motifmatchr)
  library(BSgenome.Mmusculus.UCSC.mm10)
  library(SeuratWrappers)
  library(SeuratDisk)
  library(sctransform)
  library(glmGamPoi)
  library(clustree)
  library(patchwork)
  library(qs)
  library(Scillus)
  library(scCustomize)
  library(Nebulosa)
  library(mrtree)
  library(gprofiler2)
  library(scBubbletree)
})

#' Get download markdown link
#'
#' Convert an output file name and location to a URL that can be used to
#' download the file.
#'
#' @param file name of the output file
#' @param folder name of the directory in the output directory containing the
#' output file
#'
#' @return Markdown URL link to the file
get_download_link <- function(file, folder = NULL) {
  remote <- workflowr::wflow_git_remote(verbose = FALSE)["origin"]

  url <- gsub(":", "/", remote)
  url <- gsub("git@", "http://", url)
  url <- gsub(".git", "", url, fixed = TRUE)
  url <- paste(url, "raw/main/output", sep = "/")

  if (is.null(folder)) {
    url <- paste(url, file, sep = "/")
  } else {
    url <- paste(url, folder, file, sep = "/")
  }

  link <- glue::glue("[{file}]({url})")

  return(link)
}

#' Helper function to write number of cells per file
#'
#' From the SRA id of dataset check number of cells in filtered seurat object
#'
#' @param sra id of dataset
#'
#' @return list with number of cells parameter, value and description
n_cells_per_file <- function(sra, srt = combined_srt) {
  dataset <- list(
    Parameter = sprintf("n_cells-%s", sra),
    Value = sum(srt$orig.ident == sra),
    Description = sprintf("Number of cells in the filtered %s dataset", sra)
  )
  return(dataset)
}

#' Write gene table
#'
#' Save a gene table to a file
#'
#' @param gene_table gene table to save
#' @param path file path to save location
#'
#' @details
#' data.frame objects will be saved as a (zipped) CSV file. List objects will
#' be saved in XLSX format.
write_gene_table <- function(gene_table, path) {
  if (is.data.frame(gene_table)) {
    zip_path <- paste0(path, ".zip")
    if (file.exists(zip_path)) {
      file.remove(zip_path)
    }
    readr::write_csv(gene_table, path, na = "")
    zip(zip_path, path, flags = "-q -j")
    invisible(file.remove(path))
  } else {
    writexl::write_xlsx(gene_table, path)
  }
}


#' Save plot to file
#'
#' Saves a plot to the specified path with customizable dimensions and format.
#'
#' @param name Name of the plot.
#' @param plt The plot object.
#' @param type The type of plot.
#' @param h Height of the plot (default is 12).
#' @param asp Aspect ratio of the plot (default is 1.618).
#' @param path Path to save the plot (default is `plots_dir`).
#' @param format File format of the plot (default is `.pdf`).
save_my_plot <- function(name,
                         plt,
                         type,
                         h = 12,
                         asp = 1.618,
                         path = plots_dir,
                         format = ".pdf") {
  cowplot::save_plot(
    filename = here::here(
      path,
      stringr::str_glue(type,
        as.character(name),
        format,
        .sep = "_"
      )
    ),
    plot = plt,
    base_height = h,
    base_asp = asp,
    limitsize = FALSE
  )
}

#' Save or load result of a function
#'
#' Saves the result of a computational function to a file or loads the result if it exists.
#'
#' @param filename Name of the file to save/load.
#' @param compute_function The function to compute if the file does not exist.
#' @param force_recompute Logical, whether to force recomputation (default is FALSE).
#' @return The result of the computational function.
save_or_load <- function(filename, compute_function, force_recompute = FALSE) {
  if (!dir.exists(filename)) {
    dir.create(here::here(data_dir, project))
  }
  filename <- here::here(data_dir, project, filename)
  if (!file.exists(filename) || force_recompute) {
    result <- compute_function()
    saveRDS(result, file = filename)
  } else {
    result <- readRDS(filename)
  }
  return(result)
}

#' Read in HTO-demultiplex and QC results
#'
#' Reads in HTO-demultiplexing and QC results from a CSV file.
#'
#' @param org The organism identifier.
#' @return A tibble with the HTO-demultiplexing and QC results.
read_hto_demux_qc <- function(org) {
  read_csv(here(
    data_dir,
    "COUNT",
    sprintf("%sRi_transcriptome", org),
    "QC_categories.csv"
  )) %>%
    mutate(
      Run = org,
      cell_name = str_c(org, barcode, sep = "_")
    )
}


#' Read in 10X Multiomics results from CellRanger ARC pipeline
#'
#' Reads in 10X Multiomics results from the CellRanger ARC pipeline and creates a Seurat object.
#'
#' @param sample Sample name.
#' @return A Seurat object with both RNA and ATAC data.
read_arc <- function(sample) {
  options(Seurat.object.assay.version = "v5")

  # the 10x hdf5 file contains both data types.
  # Check if we already ran import on RNA results
  rna_counts <- Read10X(
    data.dir = here(
      data_dir,
      "COUNT",
      # sprintf("%sMi_multiome", sample),
      sprintf("%sRi_transcriptome", sample),
      "filtered_feature_bc_matrix"
    )
  )
  rna_counts
  inputdata.10x <-
    Read10X_h5(here(
      data_dir,
      "COUNT",
      sprintf("%sMi_multiome", sample),
      "filtered_feature_bc_matrix.h5"
    ))
  atac_counts <- inputdata.10x$Peaks

  # Create Seurat object
  srt <- CreateSeuratObject(counts = rna_counts, project = sample)
  srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-")

  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <-
    StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <-
    seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- read_rds("~/src/data/reference_sources/EnsDb.Mmusculus.v79.rds")
  seqinfo_mm10 <- read_rds("~/src/data/reference_sources/seqinfo_ucsc_mm10.rds")

  frag.file <-
    here(
      data_dir,
      "COUNT",
      sprintf("%sMi_multiome", sample),
      "atac_fragments.tsv.gz"
    )
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = seqinfo_mm10,
    fragments = frag.file,
    min.cells = 3,
    annotation = annotations
  )
  keep_cells <-
    intersect(
      colnames(srt),
      colnames(chrom_assay)
    )
  srt <- subset(srt, cells = keep_cells)
  chrom_assay <- subset(chrom_assay, cells = keep_cells)
  srt[["ATAC"]] <- chrom_assay
  srt$cell_name <- str_c(sample, colnames(srt), sep = "_")
  srt <- RenameCells(srt, new.names = as.character(srt$cell_name))
  srt <-
    subset(srt, subset = cell_name %in% qc_categories$cell_name)
  srt@meta.data <-
    srt@meta.data |> left_join(qc_categories, by = "cell_name")
  metadata <- samples_table |>
    dplyr::select(Run, Samples:AnimalID) |>
    dplyr::filter(Run == sample)
  metadata <-
    srt@meta.data |> left_join(y = metadata, by = join_by(hto_demux == SampleAlias))
  rownames(metadata) <- metadata$cell_name
  srt@meta.data <- metadata
  srt <- RenameCells(srt, new.names = as.character(srt$cell_name))
  return(srt)
}

#' Read in 10X Multiomics results from CellRanger ARC pipeline (using BPCells)
#'
#' Reads in 10X Multiomics results from the CellRanger ARC pipeline using BPCells.
#'
#' @param sample Sample name.
#' @return A Seurat object with both RNA and ATAC data.
read_arc_bp <- function(sample) {
  options(Seurat.object.assay.version = "v5")

  # the 10x hdf5 file contains both data types.
  # Check if we already ran import on RNA results
  if (!file.exists(here(
    data_dir,
    "COUNT",
    sprintf("%s-rna-raw", sample)
  ))) {
    rna_counts <- open_matrix_10x_hdf5(
      here(
        data_dir,
        "COUNT",
        # sprintf("%sMi_multiome", sample),
        sprintf("%sRi_transcriptome", sample),
        "filtered_feature_bc_matrix.h5"
      ),
      feature_type = "Gene Expression"
    ) %>%
      write_matrix_dir(here(
        data_dir,
        "COUNT",
        sprintf("%s-rna-raw", sample)
      ))
  } else {
    rna_counts <- open_matrix_dir(here(
      data_dir,
      "COUNT",
      sprintf("%s-rna-raw", sample)
    ))
  }

  # Get all Ensembl mouse gene IDs
  ensembl_ids <- rownames(rna_counts)

  # Fetch MGI gene symbols for Ensembl mouse gene IDs
  gene_symbols <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  gene_symbols <- gene_symbols[!is.na(gene_symbols)]
  gene_names <- gene_symbols[!duplicated(gene_symbols)]
  rna_counts <- rna_counts[c(names(gene_names), "tdTomato"), ]
  rownames(rna_counts) <- c(gene_names, "tdTomato")
  rna_counts
  inputdata.10x <-
    Read10X_h5(here(
      data_dir,
      "COUNT",
      sprintf("%sMi_multiome", sample),
      "filtered_feature_bc_matrix.h5"
    ))
  atac_counts <- inputdata.10x$Peaks

  # Create Seurat object
  srt <- CreateSeuratObject(counts = rna_counts, project = sample)
  srt[["percent.mt"]] <- PercentageFeatureSet(srt, pattern = "^MT-")

  # Now add in the ATAC-seq data
  # we'll only use peaks in standard chromosomes
  grange.counts <-
    StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <-
    seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- read_rds("~/src/data/reference_sources/EnsDb.Mmusculus.v79.rds")
  seqinfo_mm10 <- read_rds("~/src/data/reference_sources/seqinfo_ucsc_mm10.rds")

  frag.file <-
    here(
      data_dir,
      "COUNT",
      sprintf("%sMi_multiome", sample),
      "atac_fragments.tsv.gz"
    )
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = seqinfo_mm10,
    fragments = frag.file,
    min.cells = 3,
    annotation = annotations
  )
  keep_cells <-
    intersect(
      colnames(srt),
      colnames(chrom_assay)
    )
  srt <- subset(srt, cells = keep_cells)
  chrom_assay <- subset(chrom_assay, cells = keep_cells)
  srt[["ATAC"]] <- chrom_assay
  srt$cell_name <- str_c(sample, colnames(srt), sep = "_")
  srt <- RenameCells(srt, new.names = as.character(srt$cell_name))
  srt <-
    subset(srt, subset = cell_name %in% qc_categories$cell_name)
  srt@meta.data <-
    srt@meta.data |> left_join(qc_categories, by = "cell_name")
  metadata <- samples_table |>
    dplyr::select(Run, Samples:AnimalID) |>
    dplyr::filter(Run == sample)
  metadata <-
    srt@meta.data |> left_join(y = metadata, by = join_by(hto_demux == SampleAlias))
  rownames(metadata) <- metadata$cell_name
  srt@meta.data <- metadata
  srt <- RenameCells(srt, new.names = as.character(srt$cell_name))
  return(srt)
}

#' Select MRTree resolution
#'
#' Selects the best MRTree resolution using ARI score for clustering.
#'
#' @param df A data frame with ARI scores and resolutions.
#' @return The selected resolution.
select_resolution <- function(df) {
  #' Get distance between two resolutions with top ARI score
  get_top_res_diff <- function(dat) {
    tmp.ari <-
      dat |>
      top_n(n = 2, wt = ari) |>
      purrr::pluck(2)
    tmp.res <-
      dat |>
      top_n(n = 2, wt = ari) |>
      purrr::pluck(1)
    tmp.ari <- tmp.ari[1] - tmp.ari[2]
    tmp.res <- tmp.res[1] - tmp.res[2]
    return(c(tmp.ari, tmp.res))
  }

  #' Pick one of two top resolutions with parameters
  pick_res_param <- function(dat,
                             ari.dif,
                             res.dif,
                             ari.thd = .05,
                             res.thd = 0) {
    if (ari.dif < ari.thd & res.dif < res.thd) {
      res <-
        dat |>
        top_n(n = 2, wt = ari) |>
        purrr::pluck(1)
      res <- res[2]
    } else {
      res <-
        dat |>
        top_n(n = 1, wt = ari) |>
        purrr::pluck(1)
    }
    return(res)
  }


  df %<>% as_tibble()
  ein.check <-
    df |>
    top_n(n = 2, wt = ari) |>
    purrr::pluck(2) |>
    purrr::map_lgl(~ .x == 1)
  if (all(ein.check)) {
    df %<>% arrange(-resolution) |> distinct(ari, .keep_all = TRUE)
    c(tmp.ari, tmp.res) %<-% get_top_res_diff(df)
    resK <- pick_res_param(df, ari.dif = tmp.ari, res.dif = tmp.res)
  } else {
    df %<>%
      dplyr::filter(ari != 1) |>
      arrange(-resolution) |>
      distinct(ari, .keep_all = TRUE)
    c(tmp.ari, tmp.res) %<-% get_top_res_diff(df)
    resK <- pick_res_param(df, ari.dif = tmp.ari, res.dif = tmp.res)
  }
  return(resK)
}

#' Compute PC score
#'
#' Computes the p-value score for each principal component.
#'
#' @param object A Seurat object.
#' @param PCs The principal components to assess.
#' @param score.thresh Threshold for p-value significance (default is 1e-05).
#' @return A data frame with the PC scores.
pc_score <-
  function(object = srt,
           PCs = 1:5,
           score.thresh = 1e-05) {
    pAll <- object[["pca"]]@jackstraw$empirical.p.values
    pAll <- pAll[, PCs, drop = FALSE]
    pAll <- as.data.frame(pAll)
    pAll$Contig <- rownames(x = pAll)
    pAll.l <- reshape2::melt(data = pAll, id.vars = "Contig")
    colnames(x = pAll.l) <- c("Contig", "PC", "Value")
    score.df <- NULL
    for (i in PCs) {
      pc.score <-
        suppressWarnings(prop.test(x = c(
          length(x = which(x = pAll[
            ,
            i
          ] <= score.thresh)),
          floor(x = nrow(x = pAll) *
            score.thresh)
        ), n = c(nrow(pAll), nrow(pAll)))$p.val)
      if (length(x = which(x = pAll[, i] <= score.thresh)) == 0) {
        pc.score <- 1
      }
      if (is.null(x = score.df)) {
        score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
      } else {
        score.df <- rbind(score.df, data.frame(PC = paste0(
          "PC",
          i
        ), Score = pc.score))
      }
    }
    return(score.df)
  }


#' Derive MRTree clustering from Seurat object
#'
#' Performs clustering using MRTree based on RNA and ATAC data, selecting the best resolution based on ARI scores.
#'
#' @param srt A Seurat object with RNA and ATAC data.
#' @param n.pcs Number of principal components to use (default is `n_pcs`).
#' @param vseed Random seed for reproducibility (default is `reseed`).
#' @param n.cores Number of cores for parallelization (default is `n_cores`).
#' @return A list containing the Seurat object and marker genes identified by LR and MAST tests.
derive_k_tree <- function(srt,
                          n.pcs = n_pcs,
                          vseed = reseed,
                          n.cores = n_cores) {
  srt <- NormalizeData(srt)
  srt <-
    FindVariableFeatures(srt,
      selection.method = "vst",
      nfeatures = 3000
    )
  all.genes <- rownames(srt)
  hvg <- VariableFeatures(srt)
  var_regex <- "^Hla-|^Ig[hjkl]|^Rna|^mt-|^Rp[sl]|^Hb[^(p)]|^Gm"
  hvg <-
    hvg[str_detect(
      pattern = var_regex,
      string = hvg,
      negate = TRUE
    )]
  srt <-
    ScaleData(
      srt,
      features = all.genes,
      vars.to.regress = c(
        "var_regex", "log10GenesPerUMI",
        "S.Score", "G2M.Score"
      )
    )
  srt <-
    RunPCA(
      srt,
      features = hvg,
      npcs = n.pcs,
      seed.use = vseed,
      verbose = TRUE
    )
  srt <-
    JackStraw(
      object = srt,
      assay = "RNA",
      reduction = "pca",
      dims = n.pcs,
      num.replicate = 100,
      prop.freq = 0.01,
      maxit = 1000
    )
  srt <-
    ScoreJackStraw(srt,
      dims = seq_along(srt[["pca"]]@stdev)
    )
  test_pc <-
    pc_score(
      object = srt,
      PCs = seq_along(srt[["pca"]]@stdev),
      score.thresh = 1e-05
    )
  selected_pcs <-
    seq_along(srt[["pca"]]@stdev)[test_pc$Score <= 1e-03 &
      srt[["pca"]]@stdev > quantile(srt[["pca"]]@stdev, .25)]
  srt <-
    srt |>
    FindNeighbors(
      dims = selected_pcs,
      k.param = 15,
      annoy.metric = "euclidean",
      n.trees = 100,
      verbose = FALSE
    ) |>
    RunUMAP(
      dims = selected_pcs,
      reduction.name = "umap",
      reduction.key = "UMAP_",
      return.model = FALSE,
      umap.method = "umap-learn",
      densmap = TRUE,
      dens.lambda = 1L,
      dens.frac = 0.3,
      n.epochs = 1000L,
      n.neighbors = 15L,
      min.dist = 0.01,
      spread = 2L,
      metric = "correlation",
      init = "pca",
      seed.use = vseed,
      verbose = FALSE
    )
  resolutions <-
    modularity_event_sampling(
      A = srt@graphs$RNA_snn,
      n.res = 30,
      gamma.min = 0.2,
      gamma.max = 2.50001
    ) # sample based on the similarity matrix
  srt <- FindClusters(
    srt,
    algorithm = 4,
    method = "igraph",
    resolution = resolutions,
    random.seed = vseed,
    verbose = FALSE
  )
  out <- mrtree(
    srt,
    prefix = "RNA_snn_res.",
    n.cores = n.cores,
    consensus = FALSE,
    sample.weighted = TRUE,
    augment.path = FALSE,
    verbose = FALSE
  )
  # Adjusted Multiresolution Rand Index (AMRI)
  ks.flat <- apply(
    out$labelmat.flat,
    2,
    FUN = function(x) {
      length(unique(x))
    }
  )
  ks.mrtree <- apply(
    out$labelmat.mrtree,
    2,
    FUN = function(x) {
      length(unique(x))
    }
  )
  amri.flat <-
    sapply(seq_len(ncol(out$labelmat.flat)), function(i) {
      AMRI(out$labelmat.flat[, i], srt$seurat_clusters)$amri
    })
  amri.flat <-
    aggregate(amri.flat, by = list(k = ks.flat), FUN = mean)
  amri.recon <-
    sapply(seq_len(ncol(out$labelmat.mrtree)), function(i) {
      AMRI(out$labelmat.mrtree[, i], srt$seurat_clusters)$amri
    })
  df <- rbind(
    data.frame(
      k = amri.flat$k,
      amri = amri.flat$x,
      method = "Seurat flat"
    ),
    data.frame(k = ks.mrtree, amri = amri.recon, method = "MRtree")
  )
  stab.out <- stability_plot(out)
  resK <- SelectResolution(stab.out$df)
  srt$k_tree <- out$labelmat.mrtree[, which.min(abs(as.integer(str_remove(
    dimnames(out$labelmat.mrtree)[[2]], "K"
  )) - resK))] %>%
    as.numeric() %>%
    as.factor()
  Idents(srt) <- "k_tree"
  if (length(unique(srt$k_tree)) > 1) {
    srt.markers.lr <-
      FindAllMarkers(
        srt,
        assay = "RNA",
        verbose = FALSE,
        random.seed = reseed,
        latent.vars = c(
          "var_regex", "log10GenesPerUMI",
          "S.Score", "G2M.Score"
        ),
        only.pos = TRUE,
        min.pct = 0.1,
        base = 10,
        logfc.threshold = 0.2,
        test.use = "LR"
      )
    if (length(unique(srt.markers.lr$cluster)) > 1) {
      write_csv(
        srt.markers.lr,
        here(
          tables_dir,
          sprintf(
            "%s_all-mrk_logreg.csv",
            unique(srt$orig.ident)
          )
        )
      )
    }
    srt.markers.mast <-
      FindAllMarkers(
        srt,
        assay = "RNA",
        verbose = FALSE,
        random.seed = reseed,
        latent.vars = c(
          "var_regex", "log10GenesPerUMI",
          "S.Score", "G2M.Score"
        ),
        only.pos = TRUE,
        min.pct = 0.1,
        base = 10,
        logfc.threshold = 0.2,
        test.use = "MAST"
      )
    if (length(unique(srt.markers.lr$cluster)) > 1) {
      write_csv(
        srt.markers.lr,
        here(
          tables_dir,
          sprintf(
            "%s_all-mrk_mast.csv",
            unique(srt$orig.ident)
          )
        )
      )
    }
  }
  return(list(srt, srt.markers.lr, srt.markers.mast))
}


#' Modify original scDEED Permuted function to make it parallel
#'
#' This function permutes the Seurat object data in parallel for subsequent DEED analysis.
#'
#' @param pbmc A Seurat object.
#' @return A permuted version of the Seurat object.
permuted_parallel <- function(pbmc) {
  plan("sequential")
  invisible(gc())
  plan(multicore, workers = n_cores)
  pbmc.permuted <- pbmc
  if (Seurat::DefaultAssay(pbmc) == "RNA") {
    X <- pbmc[["RNA"]]$scale.data
    X_permuted <- pbmc.permuted[["RNA"]]$scale.data
  } else if (Seurat::DefaultAssay(pbmc) == "SCT") {
    X <- pbmc[["SCT"]]$scale.data
    X_permuted <- pbmc.permuted[["SCT"]]$scale.data
  }

  set.seed(reseed)

  # Use future_map to parallelize the loop
  X_permuted <-
    furrr::future_pmap(list(1:dim(X)[1], list(X), list(X_permuted)), function(i, X, X_permuted) {
      row <- pracma::randperm(dim(X)[2])
      X_permuted[i, ] <- X[i, row]
      return(X_permuted[i, ])
    })
  X_permuted <- do.call(rbind, X_permuted)
  rownames(X_permuted) <- rownames(X)
  colnames(X_permuted) <- colnames(X)

  pbmc.permuted[[Seurat::DefaultAssay(pbmc)]]$scale.data <- X_permuted

  plan("sequential")
  invisible(gc())
  plan(multicore, workers = n_cores)

  return(pbmc.permuted)
}


#' Modify original scDEED Cell.Similarity.tSNE function for parallel computation
#'
#' Computes cell similarity based on t-SNE distances in parallel.
#'
#' @param Euc_distances Euclidean distances matrix.
#' @param Euc_distances_permuted Permuted Euclidean distances matrix.
#' @param tSNE_distances t-SNE distances matrix.
#' @param tSNE_distances_permuted Permuted t-SNE distances matrix.
#' @param percent Proportion of neighbors to consider for similarity calculation.
#' @return A list containing similarity results for t-SNE.
cell_similarity_tsne <- function(Euc_distances,
                                 Euc_distances_permuted,
                                 tSNE_distances,
                                 tSNE_distances_permuted,
                                 percent) {
  numberselected <- floor((dim(Euc_distances)[2]) * percent)

  rho_tSNE <- future_map_dbl(1:(dim(Euc_distances)[2]), ~ cor(
    tSNE_distances[.x, order(Euc_distances[.x, ])][2:(numberselected + 1)],
    sort(tSNE_distances[.x, ])[2:(numberselected + 1)]
  ))
  print("tSNE done")
  rho_tSNE_permuted <-
    future_map_dbl(1:(dim(Euc_distances)[2]), ~ cor(
      tSNE_distances_permuted[.x, order(Euc_distances_permuted[.x, ])][2:(numberselected + 1)],
      sort(tSNE_distances_permuted[.x, ])[2:(numberselected + 1)]
    ))
  print("permuted tSNE done")
  similarity_results_tSNE <-
    list("rho_tSNE" = rho_tSNE, "rho_tSNE_permuted" = rho_tSNE_permuted)
  return(similarity_results_tSNE)
}


#' Modify original scDEED Cell.Similarity.UMAP function for parallel computation
#'
#' Computes cell similarity based on UMAP distances in parallel.
#'
#' @param Euc_distances Euclidean distances matrix.
#' @param Euc_distances_permuted Permuted Euclidean distances matrix.
#' @param UMAP_distances UMAP distances matrix.
#' @param UMAP_distances_permuted Permuted UMAP distances matrix.
#' @param percent Proportion of neighbors to consider for similarity calculation.
#' @return A list containing similarity results for UMAP.
cell_similarity_umap <- function(Euc_distances,
                                 Euc_distances_permuted,
                                 UMAP_distances,
                                 UMAP_distances_permuted,
                                 percent) {
  numberselected <- floor((dim(Euc_distances)[2]) * percent)

  rho_UMAP <- future_map_dbl(1:(dim(Euc_distances)[2]), ~ cor(
    UMAP_distances[.x, order(Euc_distances[.x, ])][2:(numberselected + 1)],
    sort(UMAP_distances[.x, ])[2:(numberselected + 1)]
  ))
  print("UMAP done")

  rho_UMAP_permuted <-
    future_map_dbl(1:(dim(Euc_distances)[2]), ~ cor(
      UMAP_distances_permuted[.x, order(Euc_distances_permuted[.x, ])][2:(numberselected + 1)],
      sort(UMAP_distances_permuted[.x, ])[2:(numberselected + 1)]
    ))
  print("permuted UMAP done")

  similarity_results_UMAP <-
    list("rho_UMAP" = rho_UMAP, "rho_UMAP_permuted" = rho_UMAP_permuted)

  return(similarity_results_UMAP)
}


#' Modify original scDEED Cell.Classify.tSNE function for parallel computation
#'
#' Classifies cells based on their t-SNE similarity scores, identifying good, intermediate, and bad cells.
#'
#' @param rho_tSNE The similarity scores for t-SNE.
#' @param rho_tSNE_permuted The permuted similarity scores for t-SNE.
#' @return A list containing indices of bad, intermediate, and good cells.
cell_classify_tsne <- function(rho_tSNE, rho_tSNE_permuted) {
  y <- seq(1, length(rho_tSNE), length = length(rho_tSNE))

  rho_tSNE_upper <- quantile(rho_tSNE_permuted, 0.95)
  rho_tSNE_lower <- quantile(rho_tSNE_permuted, 0.05)

  tSNE_badindex <- which(rho_tSNE < rho_tSNE_lower)
  tSNE_intindex <-
    which(rho_tSNE < rho_tSNE_upper & rho_tSNE >= rho_tSNE_lower)
  tSNE_goodindex <- setdiff(y, union(tSNE_badindex, tSNE_intindex))

  ClassifiedCells_tSNE.results <-
    list(
      "tSNE_badindex" = tSNE_badindex,
      "tSNE_intindex" = tSNE_intindex,
      "tSNE_goodindex" = tSNE_goodindex
    )
  return(ClassifiedCells_tSNE.results)
}


#' Simplified scDEED function for pre-processed Seurat object
#'
#' This function performs DEED analysis on a pre-processed Seurat object, optimizing UMAP or t-SNE parameters.
#'
#' @param input_data A Seurat object.
#' @param num_pc Number of principal components to use.
#' @param n_neighbors A vector of neighbor values for UMAP.
#' @param min.dist A vector of minimum distance values for UMAP.
#' @param similarity_percent Proportion of neighbors to consider for similarity calculation.
#' @param visualization Logical, whether to visualize the results (default is FALSE).
#' @param use_method The dimensionality reduction method to use ("umap" or "tsne").
#' @param perplexity A vector of perplexity values for t-SNE.
#' @param perplexity_score Target perplexity value for t-SNE.
#' @param optimize_neib Logical, whether to optimize neighbors (default is TRUE).
#' @param optimize_min Logical, whether to optimize minimum distance (default is TRUE).
#' @return A list of results from the DEED analysis, including dubious and trustworthy cells.
scDEED_simplified <-
  function(input_data,
           num_pc,
           n_neighbors = c(
             seq(from = 5, to = 10, by = 1),
             seq(from = 15, to = 50, by = 5),
             80
           ),
           min.dist = c(
             seq(from = 0.01, to = 0.1, by = 0.025),
             seq(from = 0.2, to = 0.8, by = 0.1)
           ),
           similarity_percent = 0.5,
           visualization = FALSE,
           use_method = "umap",
           perplexity = c(seq(from = 20, to = 410, by = 20)),
           perplexity_score = 30,
           optimize_neib = TRUE,
           optimize_min = TRUE) {
    # Check if the input data is a Seurat object
    if (!"Seurat" %in% class(input_data)) {
      stop("Input data must be a Seurat object")
    }

    nsamples <- ncol(input_data)
    if (nsamples > 30000) {
      # Downsample the number of cells per identity class
      input_data <-
        subset(
          x = input_data,
          downsample = 1000
        )
    }

    if (use_method == "umap") {
      if (is.null(input_data@reductions$umap)) {
        input_data <- Seurat::RunUMAP(input_data, dims = 1:num_pc)
      }
      input_data.permuted <- permuted_parallel(input_data)
      input_data.permuted <-
        Seurat::RunPCA(
          input_data.permuted,
          npcs = num_pc,
          features = Seurat::VariableFeatures(object = input_data.permuted)
        )

      results.PCA <-
        Distances.PCA.UMAP.big(input_data, input_data.permuted, K = num_pc)


      # set up the default parameters for runumap
      default_neib <- 30L
      default_min <- 0.3

      # set up the parameters used in the final runumap.
      final_neib <- default_neib
      final_min <- default_min

      if (optimize_neib == TRUE && optimize_min == TRUE) {
        plan("sequential")
        invisible(gc())
        plan(multicore, workers = n_cores)
        # Create the grid for the parameters
        all_pairs <- expand.grid(n_neighbors, min.dist)

        # Run umap_optimize in parallel for each pair of parameters
        all_dub <- future_map2(
          .x = all_pairs$Var1,
          .y = all_pairs$Var2,
          .f = ~ umap_optimize(
            input_data = input_data,
            input_data.permuted = input_data.permuted,
            reduction.method = "pca",
            K = num_pc,
            n = .x,
            # n_neighbors
            m = .y,
            # min.dist
            results.PCA = results.PCA,
            similarity_percent = similarity_percent
          )
        ) |> simplify()
        dubious_number_UMAP <- cbind(all_pairs, all_dub)

        dub_para <- colnames(dubious_number_UMAP) <-
          c(
            "n.neighbors",
            "min.dist",
            "number of dubious cells"
          )

        # Find the best parameters
        best_para <-
          dubious_number_UMAP[which.min(all_dub), c(1, 2)]

        # In case of multiple best parameters, choose the one with the smallest sum of parameters
        if (nrow(best_para) > 1) {
          best_para <- best_para[which.min(rowSums(best_para)), ]
        }
        colnames(best_para) <- c("n.neighbors", "min.dist")

        # Get the final parameters
        final_neib <- best_para$n.neighbors
        final_min <- best_para$min.dist
      }

      if (optimize_neib == TRUE && optimize_min == FALSE) {
        dubious_number_UMAP_neib <-
          foreach::`%do%`(
            foreach::foreach(n = n_neighbors, .combine = "c"),
            umap_optimize(
              input_data = input_data,
              input_data.permuted = input_data.permuted,
              reduction.method = "pca",
              K = num_pc,
              n = n,
              m = default_min,
              results.PCA = results.PCA,
              similarity_percent = similarity_percent
            )
          )



        best_para_neib <-
          n_neighbors[which(dubious_number_UMAP_neib == min(dubious_number_UMAP_neib))]


        if (length(best_para_neib) != 0) {
          best_para_neib <- min(best_para_neib)
        }
        dub_neighbor <-
          data.frame(
            "n.neighbors" = n_neighbors,
            "number of dubious cells" = dubious_number_UMAP_neib
          )
        colnames(dub_neighbor) <-
          c("n.neighbors", "number of dubious cells")
        final_neib <- best_para_neib
      }


      if (optimize_min == TRUE && optimize_neib == FALSE) {
        dubious_number_UMAP_min <-
          foreach::`%do%`(
            foreach::foreach(m = min.dist, .combine = "c"),
            umap_optimize(
              input_data = input_data,
              input_data.permuted = input_data.permuted,
              reduction.method = "pca",
              K = num_pc,
              n = default_neib,
              m = m,
              results.PCA = results.PCA,
              similarity_percent = similarity_percent
            )
          )



        best_para_min <-
          min.dist[which(dubious_number_UMAP_min == min(dubious_number_UMAP_min))]
        if (length(best_para_min) != 0) {
          best_para_min <- min(best_para_min)
        }

        dub_min_dist <-
          data.frame(
            "min.dist" = min.dist,
            "number of dubious cells" = dubious_number_UMAP_min
          )
        colnames(dub_min_dist) <-
          c("min.dist", "number of dubious cells")
        final_min <- best_para_min
      }



      res <-
        ChoosenNeighbors(input_data,
          input_data.permuted,
          "pca",
          num_pc,
          n = final_neib,
          m = final_min
        )
      similarity_score_UMAP <-
        cell_similarity_umap(
          results.PCA$PCA_distances,
          results.PCA$PCA_distances_permuted,
          res$UMAP_distances,
          res$UMAP_distances_permuted,
          similarity_percent
        )

      ClassifiedCells_UMAP <-
        Cell.Classify.UMAP(
          similarity_score_UMAP$rho_UMAP,
          similarity_score_UMAP$rho_UMAP_permuted
        )
      cell_list <- rownames(res$object@meta.data)
      if (visualization == TRUE) {
        bad_graph <-
          Seurat::DimPlot(
            res$object,
            reduction = "umap",
            cells.highlight = list(`Dubious cells` = ClassifiedCells_UMAP$UMAP_badindex),
            cols.highlight = "red"
          )
        levels(bad_graph$data$highlight)[match("Unselected", levels(bad_graph$data$highlight))] <-
          "Other Cells"
        trust_graph <-
          Seurat::DimPlot(
            res$object,
            reduction = "umap",
            cells.highlight = list(`Trustworthy cells` = ClassifiedCells_UMAP$UMAP_goodindex),
            cols.highlight = "blue"
          )
        levels(trust_graph$data$highlight)[match("Unselected", levels(trust_graph$data$highlight))] <-
          "Other Cells"


        if (optimize_neib == TRUE && optimize_min == FALSE) {
          highlight_neib <-
            subset(dub_neighbor, n.neighbors == best_para_neib)
          input_data_dubious_plot_neib <-
            ggplot2::ggplot(
              data = dub_neighbor,
              ggplot2::aes(
                x = n.neighbors,
                y = `number of dubious cells`,
                group = 1
              )
            ) +
            ggplot2::geom_point(size = 5) +
            ggplot2::geom_point(
              data = highlight_neib,
              ggplot2::aes(x = n.neighbors, y = `number of dubious cells`),
              color = "cyan",
              size = 5
            ) +
            ggplot2::geom_vline(
              xintercept = highlight_neib$n.neighbors,
              linetype = "dotted"
            ) +
            ggplot2::annotate(
              geom = "text",
              x = highlight_neib$n.neighbors,
              y = max(dub_neighbor$`number of dubious cells`),
              label = "optimized",
              color = "cyan",
              size = 4,
              vjust = "inward",
              hjust = "inward"
            ) +
            ggplot2::labs(x = "n.neighbors", y = "# of dubious cell embeddings") +
            ggplot2::theme_bw() +
            ggplot2::theme(
              text = ggplot2::element_text(size = 20),
              panel.border = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.text = ggplot2::element_text(size = 20),
              axis.line = ggplot2::element_line(colour = "black")
            )
          output <-
            list(
              dub_neighbor,
              best_para_neib,
              cell_list[ClassifiedCells_UMAP$UMAP_badindex],
              cell_list[ClassifiedCells_UMAP$UMAP_goodindex],
              bad_graph,
              trust_graph,
              input_data_dubious_plot_neib
            )

          names(output) <-
            c(
              "number of dubious cells corresponding to n.neighbors list",
              "best n.neighbors",
              "list of dubious cells corresponding to best n.neighbors",
              "list of trustworthy cells corresponding to best n.neighbors",
              "UMAP plot with dubious cells - best n.neighbors",
              "UMAP plot with trustworthy cells - best n.neighbors",
              "plot. # of dubious embeddings vs parameters"
            )
        } else if (optimize_neib == FALSE && optimize_min == TRUE) {
          highlight_min <- subset(dub_min_dist, min.dist == best_para_min)
          input_data_dubious_plot_min <-
            ggplot2::ggplot(
              data = dub_min_dist,
              ggplot2::aes(
                x = min.dist,
                y = `number of dubious cells`,
                group = 1
              )
            ) +
            ggplot2::geom_point(size = 5) +
            ggplot2::geom_point(
              data = highlight_min,
              ggplot2::aes(x = min.dist, y = `number of dubious cells`),
              color = "cyan",
              size = 5
            ) +
            ggplot2::geom_vline(
              xintercept = highlight_min$min.dist, linetype =
                "dotted"
            ) +
            ggplot2::annotate(
              geom = "text",
              x = highlight_min$min.dist,
              y = max(dub_min_dist$`number of dubious cells`),
              label = "optimized",
              color = "cyan",
              size = 4,
              vjust = "inward",
              hjust = "inward"
            ) +
            ggplot2::labs(x = "min.dist", y = "# of dubious cell embeddings") +
            ggplot2::theme_bw() +
            ggplot2::theme(
              text = ggplot2::element_text(size = 20),
              panel.border = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.text = ggplot2::element_text(size = 20),
              axis.line = ggplot2::element_line(colour = "black")
            )
          output <-
            list(
              dub_min_dist,
              best_para_min,
              cell_list[ClassifiedCells_UMAP$UMAP_badindex],
              cell_list[ClassifiedCells_UMAP$UMAP_goodindex],
              bad_graph,
              trust_graph,
              input_data_dubious_plot_min
            )

          names(output) <-
            c(
              "number of dubious cells corresponding to min.dist list",
              "best min.dist",
              "list of dubious cells corresponding to best min.dist",
              "list of trustworthy cells corresponding to best min.dist",
              "UMAP plot with dubious cells - best min.dist",
              "UMAP plot with trustworthy cells - best min.dist",
              "plot. # of dubious embeddings vs parameters"
            )
        } else if (optimize_neib == TRUE && optimize_min == TRUE) {
          highlight <-
            subset(
              dub_para,
              n.neighbors == best_para[1, 1] &
                min.dist == best_para[1, 2]
            )
          input_data_dubious_plot <-
            ggplot2::ggplot(
              data = dub_para,
              ggplot2::aes(
                x = factor(n.neighbors),
                y = factor(min.dist),
                fill = `number of dubious cells`
              )
            ) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(low = "blue", high = "red") + # , breaks=c(min(dub_para$`number of dubious cells`), max(dub_para$`number of dubious cells`)),labels=c(min(dub_para$`number of dubious cells`) ,max(dub_para$`number of dubious cells`)), limits=c(min(dub_para$`number of dubious cells`),max(dub_para$`number of dubious cells`))) +

            ggplot2::geom_vline(xintercept = factor(highlight$n.neighbors)) +
            ggplot2::annotate(
              geom = "text",
              x = factor(highlight$n.neighbors),
              y = factor(max(dub_para[, 2])),
              label = "optimized",
              color = "cyan",
              size = 4,
              vjust = "inward",
              hjust = "inward"
            ) +
            ggplot2::geom_hline(yintercept = factor(highlight$min.dist)) +
            ggplot2::annotate(
              geom = "text",
              x = factor(max(dub_para[, 1])),
              y = factor(highlight$min.dist),
              label = "optimized",
              color = "cyan",
              size = 4,
              vjust = "inward",
              hjust = "inward"
            ) +
            ggplot2::labs(x = "n.neighbors", y = "min.dist") +
            ggplot2::theme_bw() +
            ggplot2::theme(
              text = ggplot2::element_text(size = 20),
              panel.border = ggplot2::element_blank(),
              panel.grid.major = ggplot2::element_blank(),
              panel.grid.minor = ggplot2::element_blank(),
              axis.text = ggplot2::element_text(size = 5),
              axis.line = ggplot2::element_line(colour = "black")
            )

          output <- list(
            dub_para,
            best_para,
            cell_list[ClassifiedCells_UMAP$UMAP_badindex],
            cell_list[ClassifiedCells_UMAP$UMAP_goodindex],
            bad_graph,
            trust_graph,
            input_data_dubious_plot
          )

          names(output) <-
            c(
              "number of dubious cells corresponding to pair of n.neighbors and min.dist list",
              "best pair of n.neighbors and min.dist",
              "list of dubious cells corresponding to best pair of n.neighbors and min.dist",
              "list of trustworthy cells corresponding to best pair of n.neighbors and min.dist",
              "UMAP plot with dubious cells - best pair of n.neighbors and min.dist",
              "UMAP plot with trustworthy cells - best pair of n.neighbors and min.dist",
              "plot. # of dubious embeddings vs pair of n.neighbors and min.dist"
            )
        } else {
          output <-
            list(
              cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex],
              bad_graph, trust_graph
            )
          names(output) <- c(
            "list of dubious cells",
            "list of trustworthy cells",
            "UMAP plot with dubious cells",
            "UMAP plot with trustworthy cells"
          )
        }


        return(output)
      } else {
        if (optimize_neib == TRUE && optimize_min == FALSE) {
          output <-
            list(dub_neighbor, best_para_neib, cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex])
          names(output) <-
            c(
              "number of dubious cells corresponding to n.neighbors list",
              "best n.neighbors",
              "list of dubious cells corresponding to best n.neighbors",
              "list of trustworthy cells corresponding to best n.neighbors"
            )
        } else if (optimize_neib == FALSE && optimize_min == TRUE) {
          output <-
            list(dub_min_dist, best_para_min, cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex])
          names(output) <-
            c(
              "number of dubious cells corresponding to min.dist list",
              "best min.dist",
              "list of dubious cells corresponding to best min.dist",
              "list of trustworthy cells corresponding to best min.dist"
            )
        } else if (optimize_neib == TRUE && optimize_min == TRUE) {
          output <-
            list(dub_para, best_para, cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex])
          names(output) <-
            c(
              "number of dubious cells corresponding to pair of n.neighbors and min.dist list",
              "best pair of n.neighbors and min.dist",
              "list of dubious cells corresponding to best pair of n.neighber and min.dist",
              "list of trustworthy cells corresponding to best pair of n.neighber and min.dist"
            )
        } else {
          output <-
            list(cell_list[ClassifiedCells_UMAP$UMAP_badindex], cell_list[ClassifiedCells_UMAP$UMAP_goodindex])
          names(output) <- c(
            "list of dubious cells corresponding to best n.neighber and min.dist",
            "list of trustworthy cells corresponding to best n.neighber and min.dist"
          )
        }
        return(output)
      }
    }




    if (use_method == "tsne") {
      if (nsamples < 91) {
        stop("sample size too small")
      }

      if (is.null(input_data@reductions$tsne)) {
        input_data <- Seurat::RunTSNE(input_data)
      }

      input_data.permuted <- permuted_parallel(input_data)
      results.PCA <-
        Distances.PCA.big(input_data,
          input_data.permuted,
          K = num_pc,
          perplexity_score = perplexity_score
        )

      input_data.permuted <-
        Seurat::RunPCA(
          input_data.permuted,
          npcs = num_pc,
          features = Seurat::VariableFeatures(object = input_data.permuted)
        )
      perplexity <- sort(perplexity)
      perplexity <-
        perplexity[perplexity <= floor((nsamples - 1) / 3)]




      plan("sequential")
      invisible(gc())
      plan(multicore, workers = n_cores)

      dubious_number_tSNE <-
        future_map_dbl(
          perplexity,
          ~ tsne_optimize(
            input_data,
            input_data.permuted,
            num_pc,
            perplexity = .x,
            results.PCA,
            similarity_percent
          )
        )


      best_para <-
        perplexity[which(dubious_number_tSNE == min(dubious_number_tSNE))]
      if (length(best_para) != 0) {
        best_para <- min(best_para)
      }
      dub_perplex <-
        data.frame(
          "perplexity" = perplexity,
          "number of dubious cells" = dubious_number_tSNE
        )
      colnames(dub_perplex) <-
        c("perplexity", "number of dubious cells")
      res <-
        ChoosePerplexity(input_data, input_data.permuted, num_pc, best_para)
      similarity_score_tSNE <-
        cell_similarity_tsne(
          results.PCA$PCA_distances,
          results.PCA$PCA_distances_permuted,
          res$tSNE_distances,
          res$tSNE_distances_permuted,
          similarity_percent
        )

      ClassifiedCells_tSNE <-
        cell_classify_tsne(
          similarity_score_tSNE$rho_tSNE,
          similarity_score_tSNE$rho_tSNE_permuted
        )
      cell_list <- rownames(res$object@meta.data)
      if (visualization == TRUE) {
        bad_graph <-
          Seurat::DimPlot(
            res$object,
            reduction = "tsne",
            cells.highlight = list(`Dubious cells` = ClassifiedCells_tSNE$tSNE_badindex),
            cols.highlight = "red"
          )
        levels(bad_graph$data$highlight)[match("Unselected", levels(bad_graph$data$highlight))] <-
          "Other Cells"
        trust_graph <-
          Seurat::DimPlot(
            res$object,
            reduction = "tsne",
            cells.highlight = list(`Trustworthy cells` = ClassifiedCells_tSNE$tSNE_goodindex),
            cols.highlight = "blue"
          )
        levels(trust_graph$data$highlight)[match("Unselected", levels(trust_graph$data$highlight))] <-
          "Other Cells"

        # input_data_dubious <- data.frame(perplexity, dubious_number_tSNE)
        highlight <- subset(dub_perplex, perplexity == best_para)
        input_data_dubious_plot <-
          ggplot2::ggplot(
            data = dub_perplex,
            ggplot2::aes(x = perplexity, y = `number of dubious cells`, group = 1)
          ) +
          ggplot2::geom_point(size = 5) +
          ggplot2::geom_point(
            data = highlight,
            ggplot2::aes(x = perplexity, y = `number of dubious cells`),
            color = "cyan",
            size = 5
          ) +
          ggplot2::geom_vline(
            xintercept = highlight$perplexity, linetype =
              "dotted"
          ) +
          ggplot2::annotate(
            geom = "text",
            x = highlight$perplexity,
            y = max(dub_perplex$`number of dubious cells`),
            label = "optimized",
            color = "cyan",
            size = 4,
            vjust = "inward",
            hjust = "inward"
          ) +
          ggplot2::labs(x = "perplexity", y = "# of dubious cell embeddings") +
          ggplot2::theme_bw() +
          ggplot2::theme(
            text = ggplot2::element_text(size = 20),
            panel.border = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            axis.text = ggplot2::element_text(size = 20),
            axis.line = ggplot2::element_line(colour = "black")
          )

        output <-
          list(
            dub_perplex,
            best_para,
            cell_list[ClassifiedCells_tSNE$tSNE_badindex],
            cell_list[ClassifiedCells_tSNE$tSNE_goodindex],
            bad_graph,
            trust_graph,
            input_data_dubious_plot
          )
        names(output) <-
          c(
            "number of dubious cells corresponding to perplexity list",
            "best perplexity",
            "list of dubious cell corresponding to best perplexity",
            "list of trustworthy cell corresponding to best perplexity",
            "tSNE plot with dubious cells - best perplexity",
            "tSNE plot with trustworthy cells - best perplexity",
            "plot. # of dubious embeddings vs parameters"
          )
        return(output)
      } else {
        output <-
          list(dub_perplex, best_para, cell_list[ClassifiedCells_tSNE$tSNE_badindex], cell_list[ClassifiedCells_tSNE$tSNE_goodindex])
        names(output) <-
          c(
            "number of dubious cells corresponding to perplexity list",
            "best perplexity",
            "list of dubious cell corresponding to best perplexity",
            "list of trustworthy cell corresponding to best perplexity"
          )
        return(output)
      }
    }
  }

#' Dual Filter Genes
#'
#' Filters genes based on expression and chromatin assays from a Seurat object.
#'
#' @param genes Vector of gene names to be filtered.
#' @param seurat_object A Seurat object containing the assays.
#' @param expression_assay Name of the expression assay in the Seurat object.
#' @param chromatine_assay Name of the chromatin assay in the Seurat object.
#'
#' @return A vector of filtered gene names.
#'
#' @examples
#' filtered_genes <- dual_filter_genes(genes = c("Gene1", "Gene2"), seurat_object = combined_srt)
dual_filter_genes <- function(
    genes,
    seurat_object = combined_srt,
    expression_assay = "RNA",
    chromatine_assay = "ATAC") {
  common_row_names <-
    rownames(seurat_object[[expression_assay]]$scale.data)

  genes <- genes[genes %in% common_row_names]
  genes <-
    genes |>
    map(~ try(
      Signac:::FindRegion(
        object = seurat_object,
        region = .x,
        assay = chromatine_assay
      ),
      silent = TRUE
    )) |>
    set_names(genes) |>
    map(class) |>
    simplify() %>%
    .[. == "GRanges"] |>
    names()

  return(genes)
}

#' Process markers between two conditions
#'
#' Processes marker genes for a given cluster and two diet conditions.
#'
#' @param combined_srt A Seurat object.
#' @param cluster The cluster being analyzed.
#' @param diet1 The first diet condition.
#' @param diet2 The second diet condition.
#' @param p_value_threshold The p-value threshold for filtering markers (default is 0.005).
#' @return A vector of top marker genes.
process_markers <- function(combined_srt, cluster, diet1, diet2, p_value_threshold = 0.005) {
  ident_1 <- paste0(cluster, "_", diet1)
  ident_2 <- paste0(cluster, "_", diet2)

  markers <- FindMarkers(
    combined_srt,
    ident.1 = ident_1,
    ident.2 = ident_2,
    only.pos = TRUE,
    test.use = "LR",
    min.pct = 0.05,
    latent.vars = "nCount_peaks",
    verbose = FALSE
  )

  top_markers <- markers %>%
    dplyr::filter(p_val < p_value_threshold) %>%
    rownames()

  return(top_markers)
}

#' Get genes overlapping specific peak regions
#'
#' Identifies genes that overlap with specified peak regions in the chromatin assay.
#'
#' @param combined_srt A Seurat object.
#' @param regions_to_plot A vector of peak regions to plot.
#' @param assay_name Name of the assay (default is "RNA3").
#' @return A vector of genes that overlap with the specified peak regions.
get_genes_in_peaks <- function(combined_srt, regions_to_plot, assay_name = "RNA3") {
  peak.regions <-
    combined_srt[["peaks"]]@ranges |>
    as.data.frame() |>
    mutate(n = row_number(), range = str_c(seqnames, start, end, sep = "-")) |>
    dplyr::filter(range %in% regions_to_plot) %>%
    pull(n, range) %>%
    .[regions_to_plot] %>%
    combined_srt[["peaks"]]@ranges[., ]
  peak.granges <- granges(peak.regions)

  # Now, find the overlap between peak regions and gene annotations
  overlaps <- findOverlaps(gene.annotations, peak.granges)

  # Extract the gene names that overlap with the peaks
  genes.in.peaks <- unique(mcols(gene.annotations)[queryHits(overlaps), "gene_name"])

  # Print or save the list of genes
  print(genes.in.peaks)

  genes.in.peaks <-
    genes.in.peaks[genes.in.peaks %in% rownames(combined_srt[[assay_name]]@data)]

  return(genes.in.peaks)
}

#' Plot coverage and box plots for genes
#'
#' Generates coverage and box plots for a list of genes from a Seurat object.
#'
#' @param combined_srt A Seurat object.
#' @param genes A vector of gene names to plot.
#' @param assay_name The assay to use for expression data (default is "RNA3").
#' @return A list containing coverage and box plots for the specified genes.
plot_coverage_and_box <- function(combined_srt, genes, assay_name = "RNA3") {
  coverage_plots <- lapply(genes, function(gene) {
    CoveragePlot(
      combined_srt,
      region = gene,
      features = gene,
      expression.assay = assay_name,
      expression.slot = "data",
      group.by = "DietGroup",
      split.by = "cca_clusters",
      assay = "peaks",
      extend.upstream = 5000,
      extend.downstream = 5000,
      annotation = TRUE,
      peaks = TRUE,
      tile = TRUE,
      links = TRUE,
      heights = c(16, 6, 1, 2, 3),
      widths = c(13, 1)
    )
  }) |> set_names(genes)

  box_plots <- lapply(genes, function(gene) {
    SCpubr::do_BoxPlot(
      sample = combined_srt,
      feature = gene,
      assay = assay_name,
      slot = "data",
      group.by = "cca_clusters",
      order = FALSE,
      split.by = "Diet"
    )
  }) |> set_names(genes)

  return(list(coverage_plots = coverage_plots, box_plots = box_plots))
  # return(list(coverage_plots = coverage_plots))
}

#' Find and plot motifs for top marker regions
#'
#' Identifies enriched motifs for a set of top marker regions and generates motif plots.
#'
#' @param combined_srt A Seurat object.
#' @param top_markers A vector of top marker regions.
#' @return A list containing enriched motifs and a motif plot.
find_and_plot_motifs <- function(combined_srt, top_markers) {
  enriched_motifs <- FindMotifs(
    object = combined_srt,
    features = top_markers
  )

  motif_plot <- MotifPlot(
    object = combined_srt,
    motifs = rownames(enriched_motifs)[1:25]
  )

  return(list(enriched_motifs = enriched_motifs, motif_plot = motif_plot))
}

#' Perform analysis for a given cluster and diet conditions
#'
#' This function processes marker genes, retrieves associated genes in peaks, plots coverage and box plots,
#' and performs motif analysis for a given cluster and diet conditions.
#'
#' @param combined_srt A Seurat object.
#' @param cluster The cluster being analyzed.
#' @param diet_conditions A vector of diet conditions to compare.
#' @param force_recompute Logical, whether to force recomputation of results (default is FALSE).
#' @return A list containing results of the analysis for the cluster and diet conditions.
perform_analysis_for_cluster <-
  function(combined_srt,
           cluster,
           diet_conditions,
           force_recompute = FALSE) {
    results <- list()
    for (i in seq_along(diet_conditions)) {
      for (j in i:length(diet_conditions)) {
        if (i != j) {
          results_filename <-
            paste0(
              "results_",
              cluster,
              "_",
              paste(diet_conditions, collapse = "_"),
              ".rds"
            )
          top_regions <- save_or_load(
            results_filename,
            compute_function = function() {
              process_markers(
                combined_srt = combined_srt,
                cluster = cluster,
                diet1 = diet_conditions[i],
                diet2 = diet_conditions[j]
              )
            },
            force_recompute = force_recompute
          )

          top_genes <- get_genes_in_peaks(
            combined_srt = combined_srt,
            regions_to_plot = top_regions
          )

          all_motifs <- find_and_plot_motifs(
            combined_srt = combined_srt,
            top_markers = top_regions
          )

          results[[paste0(
            diet_conditions[i],
            "_",
            diet_conditions[j],
            "_",
            cluster
          )]] <- list(
            top_genes = top_genes,
            enriched_motifs = all_motifs$enriched_motifs,
            motif_plot = all_motifs$motif_plot
          )
        }
      }
    }

    return(results)
  }

#' Convert motif to transcription factor gene
#'
#' This function retrieves the transcription factor (TF) gene names corresponding to a given motif
#' from a Seurat object's chromatin assay.
#'
#' @param tfmotif A motif or set of motifs for which to find the corresponding transcription factor genes.
#' @param srt A Seurat object containing the chromatin assay (default is `combined_srt`).
#' @return A vector of transcription factor gene names corresponding to the input motif(s).
#'
#' @details
#' This function extracts the transcription factor gene names associated with a given motif in the
#' chromatin assay of a Seurat object. It filters out names that contain punctuation and converts
#' the remaining names to title case.
#'
#' @examples
#' # Retrieve transcription factor genes for a given motif
#' tf_genes <- motif_to_tfgene(tfmotif = "MA0139.1", srt = combined_srt)
motif_to_tfgene <- function(tfmotif, srt = combined_srt) {
  tf_genes <- srt@assays$peaks@motifs@motif.names[tfmotif] |>
    map(~pluck(.x, 1)) %>%
    .[str_detect(., pattern = "[:punct:]", negate = T)] |>
    str_to_sentence()
  return(tf_genes)
}

#' Perform motif overrepresentation analysis for a cluster
#'
#' Runs motif overrepresentation analysis for a given cluster and comparison of two conditions.
#'
#' @param cluster_number The cluster number.
#' @param condition1 The first condition for comparison.
#' @param condition2 The second condition for comparison.
#' @param comparison_label A label for the comparison.
#' @return A data frame of overrepresented motifs for the cluster.
perform_overrepresentation_analysis_for_cluster <- function(cluster_number, condition1, condition2, comparison_label) {
  # Run FindMarkers for the given cluster and comparison
  result <- FindMarkers(combined_srt,
    ident.1 = paste0(cluster_number, "_", condition1),
    ident.2 = paste0(cluster_number, "_", condition2),
    only.pos = TRUE, mean.fxn = rowMeans, fc.name = "avg_diff", verbose = FALSE
  )

  # Check if there are significant results
  if (!is.null(result) && nrow(result) > 0) {
    result <- result %>%
      rownames_to_column(var = "rname") %>%
      mutate(gene = as.character(combined_srt@assays$peaks@motifs@motif.names[rname]))
    sig_motifs <- result %>%
      filter(p_val < 0.05) %>%
      slice_head(n = 6)
    if (nrow(sig_motifs) > 0) {
      sig_motifs$cluster <- cluster_number
      sig_motifs$comparison <- comparison_label
      # Store results in a list
      motifs_by_comparison[str_c(cluster_number, comparison_label)] <<- sig_motifs
      return(result)
    } else {
      return(NULL)
    }
  } else {
    return(NULL)
  }
}

#' Run motif overrepresentation analysis for multiple clusters
#'
#' Runs motif overrepresentation analysis for a list of clusters and two conditions.
#'
#' @param clusters A vector of cluster numbers.
#' @param comparison1 The first condition for comparison.
#' @param comparison2 The second condition for comparison.
#' @param comparison_label A label for the comparison.
#' @return A list of results from the overrepresentation analysis for each cluster.
run_overrepresentation_analysis_for_clusters <- function(clusters, comparison1, comparison2, comparison_label) {
  results <- lapply(clusters, function(cluster) {
    perform_overrepresentation_analysis_for_cluster(cluster, comparison1, comparison2, comparison_label)
  })

  return(results)
}

