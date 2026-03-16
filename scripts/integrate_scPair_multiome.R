#' Integrate scPair Embeddings into a Multiome Seurat Object
#'
#' This function loads a Seurat multiome object, integrates scPair embeddings,
#' performs neighbor finding, clustering, UMAP visualization, marker detection,
#' and saves cleaned outputs. Designed for GitHub workflows and reproducible pipelines.
#'
#' @param seurat_obj_path Path to input Seurat RDS object.
#' @param scpair_csv Path to scPair embedding .csv file (cells × dimensions).
#' @param metadata_csv Path to metadata containing cell barcodes.
#' @param dims_use Integer vector of embedding dimensions to use (default: 1:60).
#' @param resolution Clustering resolution (default: 0.9).
#' @param prefix Prefix for output file names (default: "Sample").
#' @param outdir Output directory (default: "./scpair_results").
#' @param split_by Metadata column name for sample-level scaling or UMAP splitting
#'                 (default: "donor_id"; set to NULL to disable split scaling).
#' @param generate_plots Whether to save UMAP/heatmap plots (default: TRUE).
#'
#' @return Processed Seurat object with scPair embedding and clustering.
#' @export
#'
#' @examples
#' integrate_scPair_multiome(
#'   seurat_obj_path="oligo_multiome.rds",
#'   scpair_csv="oligo_emb_scpair.csv",
#'   metadata_csv="oligo_metadata.csv",
#'   resolution=0.9,
#'   prefix="OL",
#'   split_by="donor_id"
#' )
#'

integrate_scPair_multiome <- function(
    seurat_obj_path,
    scpair_csv,
    metadata_csv,
    dims_use = 1:60,
    resolution = 0.9,
    prefix = "Sample",
    outdir = "./scpair_results",
    split_by = "donor_id",
    generate_plots = TRUE
){
  suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
  })

  #------------------------------------------------------------
  # 0. Prepare output directory
  #------------------------------------------------------------
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  message("Output directory: ", outdir)

  #------------------------------------------------------------
  # 1. Load Seurat object + scPair embeddings
  #------------------------------------------------------------
  message("Loading Seurat object and scPair embeddings...")
  obj <- readRDS(seurat_obj_path)

  scPair_res <- read.csv(scpair_csv, row.names = 1)
  metadata <- read.csv(metadata_csv, row.names = 1)

  # Validate cell matching
  if (!all(rownames(metadata) %in% rownames(scPair_res))) {
    warning("Not all metadata barcodes appear in scPair matrix.")
  }

  rownames(scPair_res) <- rownames(metadata)
  scPair_res <- scPair_res[colnames(obj), ]

  # Standardize embedding column names (dims)
  colnames(scPair_res) <- seq_len(length(dims_use))

  #------------------------------------------------------------
  # 2. Add scPair embedding as DimReduc object
  #------------------------------------------------------------
  message("Adding scPair embedding as DimReduc...")
  embedding <- CreateDimReducObject(
    embeddings = as.matrix(scPair_res),
    assay = "RNA",
    key = "scPair_"
  )
  obj[["scPair"]] <- embedding

  #------------------------------------------------------------
  # 3. Graph construction + clustering
  #------------------------------------------------------------
  message("Running FindNeighbors + FindClusters using scPair...")
  obj <- FindNeighbors(obj, reduction = "scPair", dims = dims_use)
  obj <- FindClusters(obj, reduction = "scPair", dims = dims_use, resolution = resolution)

  #------------------------------------------------------------
  # 4. Run UMAP
  #------------------------------------------------------------
  umap_name <- paste0("scPair_umap_res", resolution)

  obj <- RunUMAP(
    obj,
    reduction = "scPair",
    dims = dims_use,
    reduction.name = umap_name
  )

  saveRDS(obj, file.path(outdir, paste0(prefix, "_scPair_raw.rds")))

  #------------------------------------------------------------
  # 5. Normalize + scale RNA assay
  #------------------------------------------------------------
  message("Normalizing + scaling RNA assay...")

  obj <- NormalizeData(obj)

  if (!is.null(split_by) && split_by %in% colnames(obj@meta.data)) {
    message("Scaling RNA by split.by = ", split_by)
    obj <- ScaleData(
      obj,
      assay = "RNA",
      features = rownames(obj),
      split.by = split_by
    )
  } else {
    message("split.by is NULL or not found in metadata → scaling WITHOUT split")
    obj <- ScaleData(
      obj,
      assay = "RNA",
      features = rownames(obj)
    )
  }

  #------------------------------------------------------------
  # 6. Set identities
  #------------------------------------------------------------
  res_col <- paste0("RNA_snn_res.", resolution)
  Idents(obj) <- res_col

  #------------------------------------------------------------
  # 7. Generate UMAP plots (optional)
  #------------------------------------------------------------
  if (generate_plots) {
    message("Saving UMAP plots...")

    # Basic UMAP
    p1 <- DimPlot(obj, reduction = umap_name, label = TRUE)
    ggsave(
      filename = file.path(outdir, paste0(prefix, "_UMAP_res", resolution, ".pdf")),
      p1, width = 5, height = 5
    )

    # UMAP split by metadata if split_by exists
    if (!is.null(split_by) && split_by %in% colnames(obj@meta.data)) {
      message("Saving UMAP split.by = ", split_by)
      p2 <- DimPlot(obj,
                    reduction = umap_name,
                    label = TRUE,
                    group.by = res_col,
                    split.by = split_by)

      ggsave(
        filename = file.path(outdir, paste0(prefix, "_UMAP_res", resolution, "_by_", split_by, ".pdf")),
        p2, width = 12, height = 5
      )
    }
  }

  #------------------------------------------------------------
  # 8. Marker discovery
  #------------------------------------------------------------
  message("Computing markers...")
  markers <- FindAllMarkers(
    obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.5
  )

  write.csv(markers,
            file.path(outdir, paste0(prefix, "_markers_res", resolution, ".csv")))

  #------------------------------------------------------------
  # 9. Select top markers + heatmap
  #------------------------------------------------------------
  message("Selecting top markers per cluster...")

  top5 <- markers %>%
    group_by(cluster) %>%
    filter(avg_log2FC > 1) %>%
    slice_head(n = 5) %>%
    ungroup()

  if (generate_plots) {
    message("Generating heatmap...")
    h <- DoHeatmap(obj,
                   features = unique(top5$gene),
                   assay = "RNA",
                   slot = "scale.data") + NoLegend()

    ggsave(
      filename = file.path(outdir, paste0(prefix, "_heatmap_top5_res", resolution, ".png")),
      h, width = 15, height = 15, dpi = 300
    )
  }

  #------------------------------------------------------------
  # 10. Save final Seurat object
  #------------------------------------------------------------
  final_path <- file.path(outdir, paste0(prefix, "_scPair_final_res", resolution, ".rds"))
  saveRDS(obj, final_path)

  message("DONE! All results saved in: ", outdir)
  message("Final Seurat object: ", final_path)

  return(obj)
}
