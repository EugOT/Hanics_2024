import os
import random

import matplotlib
import numpy as np
import pacmap
import pandas as pd
import scanpy as sc

PLOTS_DIR = os.path.join(snakemake.params["plots"])

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["figure.max_open_warning"] = 300

reseed = 42
random.seed(reseed)
np.random.seed(reseed)
n_cores = snakemake.threads

sc.settings.verbosity = 2  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.figdir = PLOTS_DIR
sc.settings.set_figure_params(
    dpi=120, dpi_save=600, vector_friendly=True, format="pdf", transparent=True
)
sc.settings.autoshow = False
sc.settings.autosave = True
sc.logging.print_versions()


def cellbender_to_scrublet(data_path, scr_out_path, h5ad_out_path, expected_dblt):
    # load the data
    adata = sc.read_h5ad(data_path)
    sc.external.pp.scrublet(
        adata, expected_doublet_rate=expected_dblt, sim_doublet_ratio=3.0
    )
    embedding = pacmap.PaCMAP(
        n_components=2, MN_ratio=0.5, FP_ratio=2.0, apply_pca=False
    )
    adata.obsm["X_pacmap"] = embedding.fit_transform(
        adata.obsm["X_cellbender"], init="pca"
    )
    sc.pl.embedding(
        adata,
        basis="X_pacmap",
        color="doublet_score",
        title="PaCMAP: Doublets score derived using Scrublet in {sample}".format(
            sample=adata.uns["name"]
        ),
        save="_doublet-score_{sample}.pdf".format(sample=adata.uns["name"]),
    )
    sc.tl.tsne(
        adata,
        use_rep="X_cellbender",
        n_jobs=n_cores,
        random_state=reseed,
        perplexity=30,
        metric="euclidean",
    )
    sc.pl.tsne(
        adata,
        color="doublet_score",
        title="tSNE: Doublets score derived using Scrublet in {sample}".format(
            sample=adata.uns["name"]
        ),
        save="_doublet-score_{sample}.pdf".format(sample=adata.uns["name"]),
    )
    # Ensure that 'predicted_doublet' is a boolean column
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)

    # Subset the dataframe
    adata = adata[adata.obs["predicted_doublet"].isin([True, False])]

    # Write to .h5ad file
    adata.write(h5ad_out_path)
    pd.DataFrame(adata.obs).to_csv(
        scr_out_path, sep="\t", header=True
    )  # scrublet_calls.tsv


cellbender_to_scrublet(
    data_path=snakemake.input["h5ad"],
    scr_out_path=snakemake.output["scrublet_calls"],
    h5ad_out_path=snakemake.output["h5ad"],
    expected_dblt=snakemake.params["expected_dblt"],
)
