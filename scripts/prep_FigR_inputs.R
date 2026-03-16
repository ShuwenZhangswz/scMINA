#' Prepare Preprocessing Files for FigR Analysis
#'
#' This function performs preprocessing required for FigR analysis using multiome
#' scATAC + scRNA matrices, metadata, peak annotation, and a Seurat object containing
#' scPair embeddings. It generates:
#'   - A SummarizedExperiment object for ATAC counts
#'   - Normalized RNA expression matrix
#'   - kNN index matrix from scPair embeddings
#' These are saved as RDS files for downstream FigR steps.
#'
#' @param atac_mtx Path to ATAC matrix (.mtx file).
#' @param rna_mtx Path to RNA matrix (.mtx file), normalized matrix required.
#' @param metadata_csv Metadata CSV with cell annotations.
#' @param genes_csv Gene names CSV (row names correspond to RNAmat rows).
#' @param peaks_csv Peaks annotation CSV (row names correspond to ATACmat rows).
#' @param seurat_scpair_rds Path to Seurat object that contains scPair and UMAP embeddings.
#' @param k Number of neighbors for kNN graph (default = 30).
#' @param prefix Prefix for output files (default = "FigR").
#' @param outdir Output directory (default = "./FigR_preprocessing").
#'
#' @return A list containing ATAC.se, RNAmat, and cellkNN.
#' @export
#'
#' @examples
#' prep_FigR_inputs(
#'   atac_mtx = "oligo_ATACmat.mtx",
#'   rna_mtx = "oligo_RNAmat.mtx",
#'   metadata_csv = "oligo_metadata.csv",
#'   genes_csv = "oligo_genes.csv",
#'   peaks_csv = "oligo_peaks.csv",
#'   seurat_scpair_rds = "OL_scpair.rds",
#'   prefix = "oligo"
#' )
#'

prep_FigR_inputs <- function(
    atac_mtx,
    rna_mtx,
    metadata_csv,
    genes_csv,
    peaks_csv,
    seurat_scpair_rds,
    k = 30,
    prefix = "FigR",
    outdir = "./FigR_preprocessing"
){
  suppressPackageStartupMessages({
    library(dplyr)
    library(FNN)
    library(chromVAR)
    library(Matrix)
    library(SummarizedExperiment)
    library(Seurat)
    library(GenomicRanges)
    library(ggplot2)
  })

  #----------------------------------------------------------
  # Create output directory
  #----------------------------------------------------------
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  #----------------------------------------------------------
  # Load input data
  #----------------------------------------------------------
  ATACmat <- readMM(atac_mtx)
  RNAmat <- readMM(rna_mtx)
  metadata <- read.csv(metadata_csv, row.names = 1)
  genes <- read.csv(genes_csv, row.names = 1)
  peaks <- read.csv(peaks_csv, row.names = 1)

  obj_scpair <- readRDS(seurat_scpair_rds)

  #----------------------------------------------------------
  # Extract embeddings: scPair + UMAP
  #----------------------------------------------------------
  DR <- as.data.frame(Embeddings(obj_scpair[["scPair"]]))
  UMAP <- as.data.frame(Embeddings(obj_scpair[["scPair_umap"]]))

  metadata$UMAP1 <- UMAP[, 1]
  metadata$UMAP2 <- UMAP[, 2]

  #----------------------------------------------------------
  # Prepare RNA matrix
  #----------------------------------------------------------
  rownames(RNAmat) <- rownames(genes)
  colnames(RNAmat) <- rownames(metadata)

  #----------------------------------------------------------
  # Format peak annotation into GRanges
  #----------------------------------------------------------
  peak_split <- strsplit(rownames(peaks), "-")
  peak_df <- data.frame(do.call(rbind, peak_split))
  colnames(peak_df) <- c("chrom", "start", "end")
  peak_df$start <- as.numeric(peak_df$start)
  peak_df$end <- as.numeric(peak_df$end)

  peaks_granges <- GRanges(
    seqnames = Rle(peak_df$chrom),
    ranges = IRanges(start = peak_df$start, end = peak_df$end)
  )

  #----------------------------------------------------------
  # Create SummarizedExperiment for ATAC
  #----------------------------------------------------------
  ATAC.se <- SummarizedExperiment(
    assays = list(counts = ATACmat),
    rowRanges = peaks_granges,
    colData = metadata
  )

  #----------------------------------------------------------
  # Compute kNN graph (required by FigR)
  #----------------------------------------------------------
  cellkNN <- get.knn(DR, k = k)$nn.index
  rownames(cellkNN) <- rownames(DR)

  #----------------------------------------------------------
  # Save diagnostic plots (UMAP)
  #----------------------------------------------------------
  p1 <- DimPlot(obj_scpair, reduction = "scPair_umap", label = TRUE, pt.size = 0.001)
  ggsave(
    filename = file.path(outdir, paste0(prefix, "_UMAP.pdf")),
    plot = p1,
    width = 10,
    height = 10,
    dpi = 200,
    limitsize = FALSE
  )

  # UMAP colored by cluster
  df_meta <- as.data.frame(metadata)

  ggplot(df_meta, aes(UMAP1, UMAP2, color = clusters)) +
    geom_point(size = 0.5) +
    theme_classic() +
    guides(colour = guide_legend(override.aes = list(size = 2))) +
    ggsave(file.path(outdir, paste0(prefix, "_UMAP_clusters.pdf")),
           width = 10, height = 10, dpi = 200)

  # UMAP colored by disease
  if ("disease" %in% colnames(metadata)) {
    ggplot(df_meta, aes(UMAP1, UMAP2, color = disease)) +
      geom_point(size = 0.5) +
      theme_classic() +
      guides(colour = guide_legend(override.aes = list(size = 2))) +
      ggsave(file.path(outdir, paste0(prefix, "_UMAP_disease.pdf")),
             width = 10, height = 10, dpi = 200)
  }

  #----------------------------------------------------------
  # Save outputs
  #----------------------------------------------------------
  saveRDS(ATAC.se, file.path(outdir, paste0(prefix, "_ATACse_FigR.rds")))
  saveRDS(RNAmat, file.path(outdir, paste0(prefix, "_RNAmat_FigR.rds")))
  saveRDS(cellkNN, file.path(outdir, paste0(prefix, "_cellkNN_FigR.rds")))

  return(list(
    ATAC_se = ATAC.se,
    RNA_mat = RNAmat,
    cell_kNN = cellkNN
  ))
}
