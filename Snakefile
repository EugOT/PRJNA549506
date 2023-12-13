"""Mu et al., 2021"""
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
resolutions = [0.001, 0.01, 0.05, 0.1]

def plots_doublets_raw(wildcards):
    x = "output/figures/{wildcards.run}_raw/doublets_call_FPR_{wildcards.res}".format(wildcards=wildcards)
    return x.replace("\.", "_")


def get_mem_mb(wildcards, attempt):
    return attempt * 500000


##### target rules #####

shell.executable("/bin/bash")

rule all:
    input:
        expand("cellbender/{run}/{run}_output_FPR_{res}_filtered.h5",
                run=samples["Run"], res=resolutions),
        expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("cellranger/{run}/outs/filtered_feature_bc_matrix.h5",
                run=samples["Run"]),
        expand("scrublet/{run}/{run}_initial_annotation_FPR_{res}.h5ad",
                run=samples["Run"], res=resolutions),
        expand(["output/figures/combined-top5_logreg-umap-whole_dataset-fpr_{res}.pdf",
            "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_{res}.pdf",
            "output/tables/01A-eda-whole_dataset-fpr_{res}/parameters.json",
            "output/tables/01A-eda-whole_dataset-fpr_{res}/mu2021_ME_all_mrk-MAST_sct-combined-whole_dataset-fpr_{res}.csv",
            "output/tables/01A-eda-whole_dataset-fpr_{res}/mu2021_ME_all_mrk-logreg_sct-combined-whole_dataset-fpr_{res}.csv"], res=resolutions)

##### load rules #####

CELLRANGER="source /home/etretiakov/src/cellranger-7.1.0/sourceme.bash && cellranger "

rule cellranger_count:
    input:
        sample=directory("fastq"),
        idx=directory("../mm10_optimized")
    output:
        raw="cellranger/{run}/outs/raw_feature_bc_matrix.h5",
        filtered="cellranger/{run}/outs/filtered_feature_bc_matrix.h5",
        summary="cellranger/{run}/outs/web_summary.html",
        bam="cellranger/{run}/outs/possorted_genome_bam.bam",
    params:
        ids="cellranger/{run}",
        sample="{run}"
    threads: 32
    resources:
        mem_mb=64000
    shell:
        ("{run} count --include-introns true \
            --id={params.ids} \
            --sample={params.sample} \
            --transcriptome={input.idx} \
            --fastqs={input.sample} \
            --jobmode=local \
            --localcores={threads} ")

rule cellbender:
    input:
        "cellranger/{run}/outs/raw_feature_bc_matrix.h5"
    output:
        expand(["cellbender/{{run}}/{{run}}_output_FPR_{res}.h5", "cellbender/{{run}}/{{run}}_output_FPR_{res}_filtered.h5"], res=resolutions)
    params:
        ndroplets=lambda wildcards: samples["NTotalDropletsIncluded"][wildcards.run],
        ncells=lambda wildcards: samples["NTotalCells"][wildcards.run],
        h5="cellbender/{run}/{run}_output.h5"
    container:
        "docker://etretiakov/cellbender:v0.0.1"
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
            --fpr 0.001 0.01 0.05 0.1 \
            --epochs 150")

rule doublets_call:
    input:
        filt_h5="cellbender/{run}/{run}_output_FPR_{res}_filtered.h5"
    output:
        scrublet_calls="scrublet/{run}/{run}_scrublet_calls_FPR_{res}.tsv",
        dr="cellbender/{run}/{run}_latent_gene_expression_FPR_{res}.csv",
        h5ad="scrublet/{run}/{run}_initial_annotation_FPR_{res}.h5ad"
    params:
        expected_dblt=lambda wildcards: samples["NExpectedDoubletRate"][wildcards.run],
        sample_run_name="{run}",
        plots=plots_doublets_raw
    container:
        "docker://etretiakov/scrna-seq:jammy-2022.12.09-v0.0.1"
    threads: 8
    script:
        "../scrublet_cb-z.py"


rule exploratory_data_analysis_0_001:
    input:
        rmd="analysis/01A-eda-whole_dataset-fpr_0.001.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        cellbender=expand("cellbender/{run}/{run}_output_FPR_0.001_filtered.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.001.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.001.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.001.pdf",
        "output/tables/01A-eda-whole_dataset-fpr_0.001/mu2021_ME_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.001.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.001/mu2021_ME_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.001.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.001/parameters.json"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 32
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


rule exploratory_data_analysis_0_01:
    input:
        rmd="analysis/01A-eda-whole_dataset-fpr_0.01.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        cellbender=expand("cellbender/{run}/{run}_output_FPR_0.01_filtered.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.01.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.01.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.01.pdf",
        "output/tables/01A-eda-whole_dataset-fpr_0.01/mu2021_ME_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.01.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.01/mu2021_ME_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.01.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.01/parameters.json"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 32
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


rule exploratory_data_analysis_0_05:
    input:
        rmd="analysis/01A-eda-whole_dataset-fpr_0.05.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        cellbender=expand("cellbender/{run}/{run}_output_FPR_0.05_filtered.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.05.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.05.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.05.pdf",
        "output/tables/01A-eda-whole_dataset-fpr_0.05/mu2021_ME_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.05.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.05/mu2021_ME_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.05.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.05/parameters.json"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 32
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")


rule exploratory_data_analysis_0_1:
    input:
        rmd="analysis/01A-eda-whole_dataset-fpr_0.1.Rmd",
        raw=expand("cellranger/{run}/outs/raw_feature_bc_matrix.h5", run=samples["Run"]),
        cellbender=expand("cellbender/{run}/{run}_output_FPR_0.1_filtered.h5", run=samples["Run"]),
        scrublet_calls=expand("scrublet/{run}/{run}_scrublet_calls_FPR_0.1.tsv", run=samples["Run"])
    output:
        "output/figures/combined-top5_logreg-umap-whole_dataset-fpr_0.1.pdf",
        "output/figures/combined-top5_MAST-umap-whole_dataset-fpr_0.1.pdf",
        "output/tables/01A-eda-whole_dataset-fpr_0.1/mu2021_ME_all_mrk-MAST_sct-combined-whole_dataset-fpr_0.1.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.1/mu2021_ME_all_mrk-logreg_sct-combined-whole_dataset-fpr_0.1.csv",
        "output/tables/01A-eda-whole_dataset-fpr_0.1/parameters.json"
    container:
        "docker://etretiakov/workbench-session-complete:jammy-2022.12.09-custom-11.2"
    threads: 32
    shell:
        ("R -e 'workflowr::wflow_build(\"{input.rmd}\", verbose = TRUE, log_dir = here::here(\"logs_workflowr\"))'")
