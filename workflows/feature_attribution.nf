#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.adata = null
params.checkpoint_dir = null
params.cluster_source = "leiden"
params.cluster_labels = null
params.cluster_column = null
params.cov_columns = null
params.output_dir = "${projectDir}/results/attribution"
params.attribution_method = "both"
params.baseline = "both"
params.use_absolute = false
params.rank_aggregation = "none"
params.output_ranked = false
params.top_n_genes = 50
params.group_column = null

process feature_attribution {
    label 'gpu'
    label 'high_memory'

    input:
    path adata
    path checkpoint_dir
    path labels_file

    output:
    path "output/*.csv", emit: attributions

    script:
    def cov_arg = params.cov_columns ? "--cov_columns ${params.cov_columns}" : ""
    def grp_arg = params.group_column ? "--group_column ${params.group_column}" : ""
    def abs_flag = params.use_absolute ? "--use_absolute" : ""
    def rank_flag = params.output_ranked ? "--output_ranked" : ""
    def cluster_source_arg = "--cluster_source ${params.cluster_source}"
    def cluster_labels_arg = (params.cluster_source == "leiden" && params.cluster_labels) ? "--cluster_labels ${labels_file}" : ""
    def cluster_col_arg = (params.cluster_source == "obs_column" && params.cluster_column) ? "--cluster_column ${params.cluster_column}" : ""
    """
    mkdir -p output
    python ${projectDir}/scripts/run_feature_attribution.py \\
        --adata ${adata} \\
        --checkpoint_dir ${checkpoint_dir} \\
        ${cluster_source_arg} \\
        ${cluster_labels_arg} \\
        ${cluster_col_arg} \\
        --output_dir output \\
        --attribution_method ${params.attribution_method} \\
        --baseline ${params.baseline} \\
        --rank_aggregation ${params.rank_aggregation} \\
        --top_n_genes ${params.top_n_genes} \\
        ${cov_arg} \\
        ${grp_arg} \\
        ${abs_flag} \\
        ${rank_flag}
    """
}

workflow {
    if (!params.adata || !params.checkpoint_dir) {
        exit 1, "ERROR: --adata, --checkpoint_dir required"
    }
    if (params.cluster_source == "leiden" && !params.cluster_labels) {
        exit 1, "ERROR: --cluster_labels required when cluster_source=leiden"
    }
    if (params.cluster_source == "obs_column" && !params.cluster_column) {
        exit 1, "ERROR: --cluster_column required when cluster_source=obs_column"
    }
    ch_adata = file(params.adata)
    ch_ckpt = file(params.checkpoint_dir)
    ch_labels = params.cluster_source == "leiden" ? file(params.cluster_labels) : file(params.adata)
    feature_attribution(ch_adata, ch_ckpt, ch_labels)
}
