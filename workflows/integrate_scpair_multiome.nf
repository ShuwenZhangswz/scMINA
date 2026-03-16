#!/usr/bin/env nextflow

/*
 * scMINA: Integrate scPair embeddings into a multiome Seurat object
 *
 * This workflow assumes that scPair training / embedding export has already
 * been performed in the Python environment. It takes:
 *  - a Seurat multiome object (RDS)
 *  - a scPair embeddings CSV
 *  - a metadata CSV
 * and runs the R helper defined in scripts/integrate_scPair_multiome.R.
 *
 * NOTE: The R script is expected to support a simple CLI interface, e.g.
 *   Rscript integrate_scPair_multiome.R --seurat_obj_path ... --scpair_csv ...
 */

nextflow.enable.dsl = 2

// Parameters
params.seurat_obj_path = null
params.scpair_csv      = null
params.metadata_csv    = null

params.dims_use        = null          // e.g. "1:60" (handled in R)
params.resolution      = 0.9
params.prefix          = "Sample"
params.outdir          = "${projectDir}/results/integrate_scpair_multiome"
params.split_by        = "donor_id"    // set to NULL in R to disable
params.generate_plots  = true

process integrate_scpair_multiome {
    label 'high_memory'

    input:
    path seurat_obj
    path scpair_csv
    path metadata_csv

    output:
    path "*_scPair_final_res*.rds", emit: final_seurat, optional: true
    path "*_scPair_raw.rds",        emit: raw_seurat,   optional: true
    path "*_markers_res*.csv",     emit: markers,      optional: true
    path "*_UMAP*.pdf",            emit: umap_plots,   optional: true
    path "*_heatmap_top5_res*.png", emit: heatmaps,    optional: true

    script:
    def extra_args = []
    if (params.dims_use)       extra_args << "--dims_use '${params.dims_use}'"
    if (params.split_by)       extra_args << "--split_by '${params.split_by}'"
    if (!params.generate_plots) extra_args << "--no_plots"
    def extra_str = extra_args.join(' ')

    """
    mkdir -p "${params.outdir}"

    Rscript "${projectDir}/scripts/integrate_scPair_multiome.R" \\
        --seurat_obj_path "${seurat_obj}" \\
        --scpair_csv "${scpair_csv}" \\
        --metadata_csv "${metadata_csv}" \\
        --resolution ${params.resolution} \\
        --prefix "${params.prefix}" \\
        --outdir "${params.outdir}" \\
        ${extra_str}
    """
}

workflow {
    if (!params.seurat_obj_path) {
        exit 1, "ERROR: --seurat_obj_path is required"
    }
    if (!params.scpair_csv) {
        exit 1, "ERROR: --scpair_csv is required"
    }
    if (!params.metadata_csv) {
        exit 1, "ERROR: --metadata_csv is required"
    }

    def ch_seurat = file(params.seurat_obj_path)
    def ch_scpair = file(params.scpair_csv)
    def ch_meta   = file(params.metadata_csv)

    integrate_scpair_multiome(ch_seurat, ch_scpair, ch_meta)
}

