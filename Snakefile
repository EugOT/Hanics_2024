"""Hanics et al., 2024"""
import pandas as pd
from os import listdir, rename, getcwd
from os.path import join, basename, dirname, abspath
from pathlib import Path
from snakemake.utils import validate, min_version

##### set minimum snakemake version #####
min_version("7.20.0")

##### load config and sample sheets #####
configfile: "config.yaml"
samples = pd.read_table(config["samples"]).set_index("Run", drop=False)
resolutions = [0.001]
prj  = "hanics2024-cortex-tdtomato"

def plots_doublets_raw(wildcards):
    x = "output/figures/{wildcards.run}_raw/doublets_call".format(wildcards=wildcards)
    return x.replace("\.", "_")


def get_mem_mb(wildcards, attempt):
    return attempt * 500000


##### target rules #####

shell.executable("/usr/bin/bash")

rule all:
    input:
        expand("cellranger_{run}/outs/raw_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("cellranger_{run}/outs/filtered_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("cellbender/{run}_output_filtered.h5",
                run=samples["Run"]),
        expand("scrublet/{run}_initial_annotation.h5ad",
                run=samples["Run"])
        # expand(["output/tables/{prj}-whole_dataset-0.001-graph-regulons.graphml"], prj=prj)

##### load rules #####

CELLRANGER="source /home/etretiakov/src/cellranger-7.1.0/sourceme.bash && cellranger "

rule cellranger_count:
    input:
        sample=directory("fastq"),
        idx=directory("mm10_tdTomato")
    output:
        raw="cellranger_{run}/outs/raw_feature_bc_matrix.h5",
        filtered="cellranger_{run}/outs/filtered_feature_bc_matrix.h5",
        summary="cellranger_{run}/outs/web_summary.html",
        bam="cellranger_{run}/outs/possorted_genome_bam.bam",
    params:
        ids="cellranger_{run}",
        sample="{run}"
    threads: 32
    resources:
        mem_mb=64000
    shell:
        ("{CELLRANGER} count --include-introns true \
            --id={params.sample} \
            --sample={params.sample} \
            --transcriptome={input.idx} \
            --fastqs={input.sample} \
            --jobmode=local \
            --localcores={threads} ")

rule cellbender:
    input:
        "cellranger_{run}/outs/raw_feature_bc_matrix.h5"
    output:
        expand(["cellbender/{{run}}_output.h5", "cellbender/{{run}}_output_filtered.h5"], res=resolutions)
    params:
        ndroplets=lambda wildcards: samples["NTotalDropletsIncluded"][wildcards.run],
        ncells=lambda wildcards: samples["NTotalCells"][wildcards.run],
        h5="cellbender/{run}_output.h5"
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender:latest"
    threads: 4
    resources:
        nvidia_gpu=1,
        mem_mb=10000
    shell:
        ("cellbender remove-background \
            --input {input} \
            --output {params.h5} \
            --cuda \
            --expected-cells {params.ncells} \
            --total-droplets-included {params.ndroplets} \
            --fpr 0.001 \
            --epochs 150")

rule save_h5ad:
    input:
        filt_h5="cellbender/{run}_output_filtered.h5"
    output:
        dr="cellbender/{run}_latent_gene_expression.csv",
        h5ad="cellbender/{run}_filtered.h5ad"
    params:
        sample_run_name="{run}"
    container:
        "docker://us.gcr.io/broad-dsde-methods/cellbender:0.3.0"
    threads: 8
    resources:
        mem_mb=20000
    script:
        "code/cb-z.py"

rule doublets_call:
    input:
        h5ad="cellbender/{run}_filtered.h5ad"
    output:
        scrublet_calls="scrublet/{run}_scrublet_calls.tsv",
        h5ad="scrublet/{run}_initial_annotation.h5ad"
    params:
        expected_dblt=lambda wildcards: samples["NExpectedDoubletRate"][wildcards.run],
        plots=plots_doublets_raw
    container:
        "docker://etretiakov/scrna-seq:jammy-2024.10.11-v0.0.2"
    threads: 8
    resources:
        mem_mb=20000
    script:
        "code/scrublet.py"

# rule exploratory_data_analysis_0_001:
#     input:
#         rmd="analysis/01A-eda.Rmd"
#     output:
#         "data/{prj}-whole_dataset-fpr_0.001-clusters.h5ad",
#         "data/{prj}-whole_dataset-fpr_0.001-clusters.h5Seurat"
#     container:
#         "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
#     threads: 32
#     resources:
#         mem_mb=get_mem_mb
#     shell:
#         ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")

# rule export_for_regulons:
#     input:
#         h5ad="data/{prj}-whole_dataset-fpr_0.001-clusters.h5ad"
#     output:
#         expr_mtx="data/{prj}-whole_dataset-0.001-expr-mtx.csv",
#         expr_data="data/{prj}-whole_dataset-0.001-expr-mtx.tsv"
#     params:
#         prj=prj
#     container:
#         "docker://etretiakov/workbench-session-complete:jammy-2023.04.30-custom-11.7"
#     threads: 32
#     resources:
#         mem_mb=30000
#     shell:
#         ("quarto render analysis/02-export.qmd --to html")


# rule grn:
#     input:
#         expr_mtx="data/{prj}-whole_dataset-0.001-expr-mtx.csv"
#     output:
#         modules="data/{prj}-whole_dataset-0.001.adjacencies.tsv"
#     params:
#         prj=prj
#     container:
#         "docker://aertslab/pyscenic:0.12.1"
#     threads: 32
#     resources:
#         mem_mb=120000
#     shell:
#         ("pyscenic grn {input.expr_mtx} /data/references/allTFs_hg38.txt -o {output.modules} --num_workers {threads}")


# rule motifs:
#     input:
#         expr_mtx="data/{prj}-whole_dataset-0.001-expr-mtx.csv",
#         modules="data/{prj}-whole_dataset-0.001.adjacencies.tsv"
#     output:
#         motifs="data/{prj}-whole_dataset-0.001.motifs.csv"
#     params:
#         prj=prj
#     container:
#         "docker://aertslab/pyscenic:0.12.1"
#     threads: 64
#     resources:
#         mem_mb=200000
#     shell:
#         ("pyscenic ctx {input.modules} /data/references/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /data/references/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather /data/references/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /data/references/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.scores.feather --annotations_fname /data/references/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl --expression_mtx_fname {input.expr_mtx} --output {output.motifs} --num_workers {threads}")


# rule aucell:
#     input:
#         h5ad="data/{prj}-whole_dataset-fpr_0.001-clusters.h5ad",
#         motifs="data/{prj}-whole_dataset-0.001.motifs.csv"
#     output:
#         regulons="data/{prj}-whole_dataset-0.001.regulons.dat",
#         auc_mtx="data/{prj}-whole_dataset-0.001-auc-mtx.csv",
#         h5ad_scenic="data/{prj}-whole_dataset-0.001-scenic_plus.h5ad",
#         h5ad_regulons="data/{prj}-whole_dataset-0.001-regulons.h5ad",
#         gephi="output/tables/{prj}-whole_dataset-0.001-graph-regulons.graphml"
#     params:
#         prj=prj
#     container:
#         "docker://aertslab/pyscenic_scanpy:0.12.1_1.9.1"
#     threads: 32
#     resources:
#         mem_mb=200000
#     script:
#         "code/aucell.py"
