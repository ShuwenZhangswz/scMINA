#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.adata = null
params.checkpoint_dir = null
params.output_dir = "${projectDir}/results/embeddings"
params.cov = null
params.hidden_layer = "800,30"

process scpair_inference {
    label 'gpu'

    input:
    path adata
    path checkpoint_dir

    output:
    path "output/embeddings_*.csv", emit: embeddings

    script:
    """
    mkdir -p output
    python ${projectDir}/scripts/extract_embeddings.py \\
        --adata ${adata} \\
        --checkpoint_dir ${checkpoint_dir} \\
        --output_dir output \\
        ${params.cov ? "--cov '${params.cov}'" : ""} \\
        --hidden_layer ${params.hidden_layer}
    """
}

workflow {
    if (!params.adata || !params.checkpoint_dir) {
        exit 1, "ERROR: --adata and --checkpoint_dir required"
    }
    ch_adata = file(params.adata)
    ch_ckpt = file(params.checkpoint_dir)
    scpair_inference(ch_adata, ch_ckpt)
}
