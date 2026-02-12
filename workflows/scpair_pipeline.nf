#!/usr/bin/env nextflow
/*
 * scMINA scPair Feature Attribution Pipeline
 * train -> embeddings -> scpair_cluster
 * Feature attribution: run feature_attribution.nf separately after choosing resolution
 *
 * --run full: train then scpair_cluster (default)
 * --run clustering_only: scpair_cluster only, requires --embeddings_dir
 */
nextflow.enable.dsl = 2

params.run = "full"
params.input_mode = "h5ad"
params.input_h5ad = null
params.input_rna = null
params.input_atac = null
params.input_meta = null
params.input_index_dir = null
params.cov = null
params.output_dir = "${projectDir}/results"
params.save_adata = true
params.seed = 0
params.hidden_layer = "800,30"
params.max_epochs = 1000
params.resolutions = "0.5,0.9,1.2,1.5"
params.save_concat = false
params.embeddings_dir = null

include { scpair_train } from './scpair_train.nf'
include { scpair_cluster } from './scpair_cluster.nf'

workflow {
    if (params.run == "full") {
        if (params.input_mode == 'h5ad' && !params.input_h5ad) {
            exit 1, "ERROR: --input_h5ad required when run=full and input_mode=h5ad"
        }
        if (params.input_mode == 'separate' && (!params.input_rna || !params.input_atac || !params.input_meta || !params.input_index_dir)) {
            exit 1, "ERROR: separate mode requires --input_rna, --input_atac, --input_meta, --input_index_dir"
        }
        scpair_train()
        scpair_cluster(scpair_train.out.output_dir)
    } else if (params.run == "clustering_only") {
        if (!params.embeddings_dir) exit 1, "ERROR: --embeddings_dir required when run=clustering_only"
        scpair_cluster(file(params.embeddings_dir))
    } else {
        exit 1, "ERROR: run must be 'full' or 'clustering_only'"
    }
}
