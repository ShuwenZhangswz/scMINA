#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.embeddings_dir = null
params.output_dir = "${projectDir}/results/clustering"
params.resolutions = "0.5,0.9,1.2,1.5"
params.embedding_mode = "all"
params.save_concat = false

process scpair_cluster {
    label 'high_memory'

    input:
    path embeddings_dir

    output:
    path "output/cluster_labels_res*.csv", emit: cluster_labels
    path "output/concat_embeddings.csv", optional: true, emit: concat_emb

    script:
    def save_flag = params.save_concat ? "--save_concat" : ""
    """
    mkdir -p output
    python ${projectDir}/scripts/run_clustering.py \\
        --embeddings_dir ${embeddings_dir} \\
        --output_dir output \\
        --resolutions ${params.resolutions} \\
        --embedding_mode ${params.embedding_mode} \\
        ${save_flag}
    """
}

workflow {
    if (!params.embeddings_dir) {
        exit 1, "ERROR: --embeddings_dir required"
    }
    ch_emb = file(params.embeddings_dir)
    scpair_cluster(ch_emb)
}
