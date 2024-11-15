#!/usr/bin/env Rscript
# data.R

source(here(src_dir, "genes.R"))

# samples_table <-
#   readr::read_tsv(here(data_dir, "cellranger_metrics_summary.tsv")) %>% arrange(Run)
# srr_set <- unique(samples_table$Run)
#
# qc_categories <-
#   purrr::reduce(
#     srr_set %>% map(~ read_hto_demux_qc(.x)),
#     bind_rows
#   ) %>%
#   as.data.frame()
# rownames(qc_categories) <- qc_categories$cell_name

load_raw_srt <- function(project) {
  if (!file.exists(
    here(
      data_dir,
      glue::glue("{project}-raw/{project}-raw.Rds")
    )
  )) {
    options(Seurat.object.assay.version = "v5")

    srt_list <- samples_table |>
      dplyr::filter(Region == project) |>
      dplyr::distinct(Run) |>
      dplyr::pull(Run) |>
      map(~ read_arc_bp(.x)) |>
      setNames(
        samples_table |>
          dplyr::filter(Region == project) |>
          dplyr::distinct(Run) |>
          dplyr::pull(Run)
      )

    combined_srt <-
      merge(srt_list[[1]], y = srt_list[2:length(srt_list)])

    dir.create(here(data_dir, sprintf("%s-raw", project)))
    SeuratObject::SaveSeuratRds(
      combined_srt,
      file = here(
        data_dir,
        glue::glue("{project}-raw/{project}-raw.Rds")
      )
    )

    plan(sequential)
    invisible(gc())
    options(future.globals.maxSize = 999999 * 1024^2)
    set.seed(seed = reseed)
    plan(multisession, workers = n_cores)
  }

  combined_srt <-
    read_rds(here(
      data_dir,
      glue::glue("{project}-raw/{project}-raw.Rds")
    ))

  return(combined_srt)
}

load_raw_srt_norm <- function(project) {
  options(Seurat.object.assay.version = "v5")

  srt_list <- samples_table |>
    dplyr::filter(Region == project) |>
    dplyr::distinct(Run) |>
    dplyr::pull(Run) |>
    map(~ read_arc(.x)) |>
    setNames(
      samples_table |>
        dplyr::filter(Region == project) |>
        dplyr::distinct(Run) |>
        dplyr::pull(Run)
    )

  combined_srt <-
    merge(srt_list[[1]], y = srt_list[2:length(srt_list)])

  return(combined_srt)
}

load_init_filt_srt <- function(project) {
  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-init-filt/{project}-init-filt.Rds")
  ))) {
    combined_srt <- load_raw_srt(project)

    plan(sequential)
    invisible(gc())
    options(future.globals.maxSize = 999999 * 1024^2)
    set.seed(seed = reseed)
    plan(multisession, workers = n_cores)

    combined_srt <- preprocess_raw_srt(project, combined_srt)

    dir.create(here(data_dir, sprintf("%s-init-filt", project)))
    SeuratObject::SaveSeuratRds(
      combined_srt,
      file = here(
        data_dir,
        glue::glue("{project}-init-filt/{project}-init-filt.Rds")
      )
    )
  } else {
    combined_srt <-
      read_rds(here(
        data_dir,
        glue::glue("{project}-init-filt/{project}-init-filt.Rds")
      ))
  }

  return(combined_srt)
}

load_init_filt_srt_norm <- function(project, metadata) {
  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-norm/{project}-norm.Rds")
  ))) {
    combined_srt <- load_raw_srt_norm(project)

    plan(sequential)
    invisible(gc())
    options(future.globals.maxSize = 999999 * 1024^2)
    set.seed(seed = reseed)
    plan(multisession, workers = n_cores)

    combined_srt <- preprocess_srt_norm(project, combined_srt, metadata)

    dir.create(here(data_dir, sprintf("%s-norm", project)))
    SeuratObject::SaveSeuratRds(
      combined_srt,
      file = here(
        data_dir,
        glue::glue("{project}-norm/{project}-norm.Rds")
      )
    )
  } else {
    combined_srt <-
      SeuratObject::LoadSeuratRds(here(
        data_dir,
        glue::glue("{project}-norm/{project}-norm.Rds")
      ))
  }

  return(combined_srt)
}

load_init_qc_srt <- function(project) {
  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-init/{project}-init.Rds")
  ))) {
    combined_srt <- load_init_filt_srt(project)

    plan(sequential)
    invisible(gc())
    options(future.globals.maxSize = 999999 * 1024^2)
    set.seed(seed = reseed)
    plan(multisession, workers = n_cores)

    combined_srt <- preprocess_init_filt_srt(project, combined_srt)
    dir.create(here(data_dir, sprintf("%s-init", project)))
    SeuratObject::SaveSeuratRds(
      combined_srt,
      file = here(
        data_dir,
        glue::glue("{project}-init/{project}-init.Rds")
      )
    )
  } else {
    combined_srt <-
      read_rds(here(
        data_dir,
        glue::glue("{project}-init/{project}-init.Rds")
      ))
  }

  return(combined_srt)
}

load_post_qc_srt <- function(project) {
  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-postqc/{project}-postqc.Rds")
  ))) {
    combined_srt <- load_init_qc_srt(project)

    plan(sequential)
    invisible(gc())
    options(future.globals.maxSize = 999999 * 1024^2)
    set.seed(seed = reseed)
    plan(multisession, workers = n_cores)

    combined_srt <- preprocess_post_qc_srt(project, combined_srt)

    dir.create(here(data_dir, sprintf("%s-postqc", project)))
    SeuratObject::SaveSeuratRds(
      combined_srt,
      file = here(
        data_dir,
        glue::glue("{project}-postqc/{project}-postqc.Rds")
      )
    )
  } else {
    combined_srt <-
      read_rds(here(
        data_dir,
        glue::glue("{project}-postqc/{project}-postqc.Rds")
      ))
  }

  return(combined_srt)
}

load_sct_srt <- function(project) {
  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-sct/{project}-sct.Rds")
  ))) {
    combined_srt <- load_post_qc_srt(project)

    plan(sequential)
    invisible(gc())
    options(future.globals.maxSize = 999999 * 1024^2)
    set.seed(seed = reseed)
    plan(multisession, workers = n_cores)

    combined_srt <- preprocess_sct_srt(project, combined_srt)

    dir.create(here(data_dir, sprintf("%s-sct", project)))
    SeuratObject::SaveSeuratRds(
      object = combined_srt,
      file = here(
        data_dir,
        glue::glue("{project}-sct/{project}-sct.Rds")
      )
    )
  } else {
    combined_srt <-
      SeuratObject::LoadSeuratRds(here(
        data_dir,
        glue::glue("{project}-sct/{project}-sct.Rds")
      ))
  }

  return(combined_srt)
}

load_integrated_transcriptome_srt <- function(project) {
  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-integrated/{project}-integrated.Rds")
  ))) {
    unintegrated_srt <- load_sct_srt(project)
    combined_srt <- load_init_filt_srt_norm(project, unintegrated_srt@meta.data)

    # plan(sequential)
    invisible(gc())
    options(future.globals.maxSize = 999999 * 1024^2)
    set.seed(seed = reseed)
    # plan(multisession, workers = n_cores)

    combined_srt <- integrate_srt(project, combined_srt, unintegrated_srt)

    dir.create(here(data_dir, sprintf("%s-integrated", project)))
    SeuratObject::SaveSeuratRds(
      combined_srt,
      file = here(
        data_dir,
        glue::glue("{project}-integrated/{project}-integrated.Rds")
      )
    )
  } else {
    combined_srt <-
      SeuratObject::LoadSeuratRds(here(data_dir, glue::glue("{project}-integrated/{project}-integrated.Rds")))
  }

  return(combined_srt)
}

load_sct_markers <- function(project) {
  markers_logreg <- readr::read_csv(here(
    tables_dir,
    sprintf(
      "%s_all_mrk-logreg_sct-combined-whole_dataset.csv",
      project
    )
  ))

  markers_MAST <- readr::read_csv(here(
    tables_dir,
    sprintf(
      "%s_all_mrk-MAST_sct-combined-whole_dataset.csv",
      project
    )
  ))


  return(list(markers_logreg, markers_MAST))
}

load_integrated_markers <- function(project) {
  if (!file.exists(here(
    data_dir,
    glue::glue("{project}-integrated/{project}-integrated-markers.Rds")
  ))) {
    combined_srt <- load_integrated_transcriptome_srt(project)

    plan(sequential)
    invisible(gc())
    options(future.globals.maxSize = 999999 * 1024^2)
    set.seed(seed = reseed)
    plan(multisession, workers = n_cores)

    markers <- c(
      "cca_clusters", "rpca_clusters", "harmony_clusters", "mnn_clusters"
    ) |>
      map(quantify_de_integrated) |>
      set_names(c(
        "cca_clusters", "rpca_clusters", "harmony_clusters", "mnn_clusters"
      ))

    dir.create(here(data_dir, sprintf("%s-integrated", project)))
    write_rds(
      markers,
      file = here(
        data_dir,
        glue::glue("{project}-integrated/{project}-integrated-markers.Rds")
      )
    )
  } else {
    markers <-
      read_rds(here(data_dir, glue::glue("{project}-integrated/{project}-integrated-markers.Rds")))
  }

  return(markers)
}
