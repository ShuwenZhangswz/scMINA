#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.input_mode = "h5ad"
params.input_h5ad = null
params.input_rna = null
params.input_atac = null
params.input_meta = null
params.input_index_dir = null
params.cov = null
params.output_dir = "${projectDir}/results/scpair_train"
params.save_adata = false
params.export_embeddings = true
params.seed = 0
params.hidden_layer = "800,30"
params.max_epochs = 1000

process scpair_train {
    label 'gpu'
    label 'high_memory'

    output:
    path "output", emit: output_dir
    path "output/checkpoints/*", emit: checkpoints
    path "output/embeddings_*.csv", emit: embeddings
    path "output/train_metadata.csv", emit: train_metadata
    path "output/cov_columns.json", emit: cov_columns
    path "output/adata_paired.h5ad", optional: true, emit: adata

    script:
    def mode = params.input_mode
    def out = "output"
    def extra = params.cov ? " --cov '${params.cov}'" : ""
    if (params.save_adata) extra += " --save_adata"
    if (!params.export_embeddings) extra += " --no_export_embeddings"
    def input_args = ""
    if (mode == 'h5ad' && params.input_h5ad) {
        input_args = " --input_h5ad ${params.input_h5ad}"
    } else if (mode == 'separate') {
        if (params.input_rna) input_args += " --input_rna ${params.input_rna}"
        if (params.input_atac) input_args += " --input_atac ${params.input_atac}"
        if (params.input_meta) input_args += " --input_meta ${params.input_meta}"
        if (params.input_index_dir) input_args += " --input_index_dir ${params.input_index_dir}"
    }
    """
    mkdir -p ${out} ${out}/checkpoints
    python ${projectDir}/scripts/run_scpair_train.py \\
        --input_mode ${mode} \\
        --output_dir ${out} \\
        --seed ${params.seed} \\
        --hidden_layer ${params.hidden_layer} \\
        --max_epochs ${params.max_epochs} \\
        ${input_args} \\
        ${extra}
    """
}

workflow {
    if (params.input_mode == 'h5ad' && !params.input_h5ad) {
        exit 1, "ERROR: --input_h5ad required when input_mode=h5ad"
    }
    if (params.input_mode == 'separate') {
        if (!params.input_rna || !params.input_atac || !params.input_meta || !params.input_index_dir) {
            exit 1, "ERROR: separate mode requires --input_rna, --input_atac, --input_meta, --input_index_dir"
        }
    }
    scpair_train()
}
