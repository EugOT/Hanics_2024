import numpy as np
from cellbender.remove_background.downstream import anndata_from_h5


def cellbender_to_h5ad(data_path, cb_z_out_path, h5ad_out_path, sample_name):
    # load the data
    adata = anndata_from_h5(data_path)
    adata.uns["name"] = sample_name
    adata.obsm["X_cellbender"] = adata.obsm["gene_expression_encoding"]
    adata.var_names_make_unique()
    # load the latent representation into a new slot called 'X_cellbender'
    z = []
    z = adata.obsm["X_cellbender"]
    np.savetxt(
        cb_z_out_path, z, delimiter=","
    )  # save the latent representation to a csv file

    # Save results:
    adata.write(h5ad_out_path)


cellbender_to_h5ad(
    data_path=snakemake.input["filt_h5"],
    cb_z_out_path=snakemake.output["dr"],
    h5ad_out_path=snakemake.output["h5ad"],
    sample_name=snakemake.params["sample_run_name"],
)
