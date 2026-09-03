#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.seurat_obj = null
params.deg_rds = null
params.enrichment = null
params.figr_rds = null
params.genes = null
params.plots = "coverage,volcano,enrichment,network,figr_heatmap"
params.prefix = "scMINA"
params.seed = 0
params.padj_cutoff = 0.05
params.log2fc_cutoff = 1.0
params.top_terms_per_celltype = 2
params.enrichment_direction = "all"
params.network_score_cutoff = 1.5
params.heatmap_score_cutoff = 0.8
params.output_dir = "${projectDir}/results/extra_visualizations"

process extra_visualizations {
    label 'high_memory'

    input:
    path seurat_obj
    path deg_rds
    path enrichment
    path figr_rds

    output:
    path "output/*", emit: figures

    script:
    """
    mkdir -p output

    Rscript "${projectDir}/scripts/run_extra_visualizations.R" \
        --seurat_obj "${seurat_obj}" \
        --deg_rds "${deg_rds}" \
        --enrichment "${enrichment}" \
        --figr_rds "${figr_rds}" \
        --genes "${params.genes}" \
        --plots "${params.plots}" \
        --prefix "${params.prefix}" \
        --seed ${params.seed} \
        --padj_cutoff ${params.padj_cutoff} \
        --log2fc_cutoff ${params.log2fc_cutoff} \
        --top_terms_per_celltype ${params.top_terms_per_celltype} \
        --enrichment_direction "${params.enrichment_direction}" \
        --network_score_cutoff ${params.network_score_cutoff} \
        --heatmap_score_cutoff ${params.heatmap_score_cutoff} \
        --output_dir output
    """
}

workflow {
    if (!params.seurat_obj) exit 1, "ERROR: --seurat_obj is required"
    if (!params.deg_rds) exit 1, "ERROR: --deg_rds is required"
    if (!params.enrichment) exit 1, "ERROR: --enrichment is required"
    if (!params.figr_rds) exit 1, "ERROR: --figr_rds is required"
    if (!params.genes) exit 1, "ERROR: --genes is required (comma-separated)"

    extra_visualizations(
        file(params.seurat_obj),
        file(params.deg_rds),
        file(params.enrichment),
        file(params.figr_rds)
    )
}
