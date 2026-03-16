#!/usr/bin/env nextflow

/*
 * scMINA FigR pipeline
 *
 * Combined Nextflow workflow that:
 *   1) Prepares FigR inputs from multiome matrices and a Seurat object
 *      with scPair embeddings (prep_FigR_inputs.R)
 *   2) Runs FigR GRN analysis using those inputs (run_FigR_analysis.R)
 *
 * Both steps are implemented in R under scMINA/scripts/.
 * The R scripts are expected to support a simple CLI interface.
 */

nextflow.enable.dsl = 2

// ---------------------------- Parameters ----------------------------

// Preprocessing inputs
params.atac_mtx          = null   // ATAC matrix (.mtx)
params.rna_mtx           = null   // normalized RNA matrix (.mtx)
params.metadata_csv      = null   // metadata CSV (barcodes as rownames)
params.genes_csv         = null   // genes CSV
params.peaks_csv         = null   // peaks CSV
params.seurat_scpair_rds = null   // Seurat object with scPair + UMAP

// Preprocessing options
params.k           = 30
params.prefix      = "FigR"
params.prep_outdir = "${projectDir}/results/FigR_preprocessing"

// Analysis options
params.genome           = "hg38"
params.nCores_corr      = 64
params.nCores_smooth    = 24
params.nCores_grn       = 8
params.dorc_cutoff      = 5
params.dorc_labelTop    = 20
params.dorcK            = 2
params.analysis_prefix  = null       // defaults to params.prefix if null
params.analysis_outdir  = "${projectDir}/results/FigR_results"

// ---------------------------- Processes ----------------------------

process prep_figr_inputs {
    label 'high_memory'

    input:
    path atac_mtx
    path rna_mtx
    path metadata_csv
    path genes_csv
    path peaks_csv
    path seurat_scpair_rds

    output:
    path "*_ATACse_FigR.rds", emit: atac_se
    path "*_RNAmat_FigR.rds", emit: rna_mat
    path "*_cellkNN_FigR.rds", emit: cell_knn
    path "*_UMAP*.pdf",       emit: prep_umap_plots, optional: true

    script:
    """
    mkdir -p "${params.prep_outdir}"

    Rscript "${projectDir}/scripts/prep_FigR_inputs.R" \\
        --atac_mtx "${atac_mtx}" \\
        --rna_mtx "${rna_mtx}" \\
        --metadata_csv "${metadata_csv}" \\
        --genes_csv "${genes_csv}" \\
        --peaks_csv "${peaks_csv}" \\
        --seurat_scpair_rds "${seurat_scpair_rds}" \\
        --k ${params.k} \\
        --prefix "${params.prefix}" \\
        --outdir "${params.prep_outdir}"
    """
}

process run_figr_analysis {
    label 'high_memory'

    input:
    path atac_se_rds
    path rna_mat_rds
    path cellknn_rds

    output:
    path "*_cisCorr.rds",        emit: cisCorr,      optional: true
    path "*_numDorcs.rds",       emit: numDorcs,     optional: true
    path "*_FigR_GRN.rds",       emit: figR_GRN,     optional: true
    path "*_DORC_plot.pdf",      emit: dorc_plot,    optional: true
    path "*_TF_DORC_scatter.pdf", emit: tf_dorc_scatter, optional: true
    path "*_TF_drivers.pdf",     emit: tf_drivers,   optional: true

    script:
    def prefix = params.analysis_prefix ?: params.prefix
    """
    mkdir -p "${params.analysis_outdir}"

    Rscript "${projectDir}/scripts/run_FigR_analysis.R" \\
        --atac_se_rds "${atac_se_rds}" \\
        --rna_mat_rds "${rna_mat_rds}" \\
        --cellknn_rds "${cellknn_rds}" \\
        --genome "${params.genome}" \\
        --nCores_corr ${params.nCores_corr} \\
        --nCores_smooth ${params.nCores_smooth} \\
        --nCores_grn ${params.nCores_grn} \\
        --dorc_cutoff ${params.dorc_cutoff} \\
        --dorc_labelTop ${params.dorc_labelTop} \\
        --dorcK ${params.dorcK} \\
        --prefix "${prefix}" \\
        --outdir "${params.analysis_outdir}"
    """
}

// ---------------------------- Workflow ----------------------------

workflow {
    // Validate required preprocessing inputs
    if (!params.atac_mtx)          exit 1, "ERROR: --atac_mtx is required"
    if (!params.rna_mtx)           exit 1, "ERROR: --rna_mtx is required"
    if (!params.metadata_csv)      exit 1, "ERROR: --metadata_csv is required"
    if (!params.genes_csv)         exit 1, "ERROR: --genes_csv is required"
    if (!params.peaks_csv)         exit 1, "ERROR: --peaks_csv is required"
    if (!params.seurat_scpair_rds) exit 1, "ERROR: --seurat_scpair_rds is required"

    def ch_atac    = file(params.atac_mtx)
    def ch_rna     = file(params.rna_mtx)
    def ch_meta    = file(params.metadata_csv)
    def ch_genes   = file(params.genes_csv)
    def ch_peaks   = file(params.peaks_csv)
    def ch_seurat  = file(params.seurat_scpair_rds)

    // Step 1: Preprocessing
    prep_figr_inputs(ch_atac, ch_rna, ch_meta, ch_genes, ch_peaks, ch_seurat)

    // Step 2: FigR analysis, using outputs from preprocessing
    run_figr_analysis(
        prep_figr_inputs.out.atac_se,
        prep_figr_inputs.out.rna_mat,
        prep_figr_inputs.out.cell_knn
    )
}

