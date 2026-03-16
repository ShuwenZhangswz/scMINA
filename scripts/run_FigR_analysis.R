#' Run FigR Gene Regulatory Network Analysis
#'
#' This function performs the full FigR analysis workflow, including:
#'   - Peak–gene association testing (runGenePeakcorr)
#'   - DORC identification and visualization
#'   - DORC scoring and smoothing
#'   - RNA smoothing using kNN graph
#'   - GRN inference using FigR::runFigRGRN
#'   - TF-DORC visualization and driver ranking
#'
#' It assumes that preprocessing was performed using prep_FigR_inputs().
#'
#' @param atac_se_rds Path to SummarizedExperiment object produced during preprocessing.
#' @param rna_mat_rds Path to normalized RNA matrix (RDS).
#' @param cellknn_rds Path to cell kNN matrix produced during preprocessing.
#' @param genome Genome assembly for FigR ("hg38", "hg19", or "mm10"). Default: "hg38".
#' @param nCores_corr Number of cores for runGenePeakcorr.
#' @param nCores_smooth Number of cores for smoothing steps.
#' @param nCores_grn Number of cores for GRN inference.
#' @param dorc_cutoff Minimum number of significant peaks to call a DORC (default = 5).
#' @param dorc_labelTop Number of DORCs to label in plots (default = 20).
#' @param dorcK Proportion for DORC selection in GRN (default = 2).
#' @param prefix Output filename prefix (default = "FigR").
#' @param outdir Output directory (default = "./FigR_results").
#'
#' @return A list containing cisCorr, cisCorr.filt, dorcGenes, dorcMat.s, RNAmat.s, and figR.d.
#' @export
#'
#' @examples
#' run_FigR_analysis(
#'   atac_se_rds = "oligo_ATACse_FigR.rds",
#'   rna_mat_rds = "oligo_RNAmat_FigR.rds",
#'   cellknn_rds = "oligo_cellkNN_FigR.rds",
#'   genome = "hg38",
#'   prefix = "oligo"
#' )
#'

run_FigR_analysis <- function(
    atac_se_rds,
    rna_mat_rds,
    cellknn_rds,
    genome = "hg38",
    nCores_corr = 64,
    nCores_smooth = 24,
    nCores_grn = 8,
    dorc_cutoff = 5,
    dorc_labelTop = 20,
    dorcK = 2,
    prefix = "FigR",
    outdir = "./FigR_results"
){
  suppressPackageStartupMessages({
    library(dplyr)
    library(Matrix)
    library(FigR)
    library(BuenColors)
    library(ggplot2)
  })

  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  #------------------------------------------------------
  # Load Inputs
  #------------------------------------------------------
  ATAC.se <- readRDS(atac_se_rds)
  RNAmat  <- readRDS(rna_mat_rds)
  cellkNN <- readRDS(cellknn_rds)

  #------------------------------------------------------
  # Peak–Gene Correlation Testing
  #------------------------------------------------------
  cisCorr <- FigR::runGenePeakcorr(
    ATAC.se = ATAC.se,
    RNAmat = RNAmat,
    genome = genome,
    nCores = nCores_corr,
    p.cut = NULL,
    n_bg = 100
  )
  saveRDS(cisCorr, file.path(outdir, paste0(prefix, "_cisCorr.rds")))

  cisCorr.filt <- cisCorr %>% dplyr::filter(pvalZ <= 0.05)

  #------------------------------------------------------
  # DORC Identification + Visualization
  #------------------------------------------------------
  pdf(file.path(outdir, paste0(prefix, "_DORC_plot.pdf")), width = 10, height = 5)
  dorcGenes <- dorcJPlot(
    dorcTab = cisCorr.filt,
    cutoff = dorc_cutoff,
    labelTop = dorc_labelTop,
    returnGeneList = TRUE,
    force = 2
  )
  dev.off()

  numDorcs <- cisCorr.filt %>%
    group_by(Gene) %>%
    tally() %>%
    arrange(desc(n))
  saveRDS(numDorcs, file.path(outdir, paste0(prefix, "_numDorcs.rds")))

  #------------------------------------------------------
  # Compute DORC Scores
  #------------------------------------------------------
  dorcMat <- getDORCScores(
    ATAC.se = ATAC.se,
    dorcTab = cisCorr.filt,
    geneList = dorcGenes,
    nCores = 1
  )

  #------------------------------------------------------
  # Smooth DORC Scores Using kNN
  #------------------------------------------------------
  dorcMat.s <- smoothScoresNN(
    NNmat = cellkNN[, 1:30],
    mat = dorcMat,
    nCores = 4
  )

  #------------------------------------------------------
  # Smooth RNA Matrix Using kNN
  #------------------------------------------------------
  RNAmat.s <- smoothScoresNN(
    NNmat = cellkNN[, 1:30],
    mat = RNAmat,
    nCores = nCores_smooth
  )

  #------------------------------------------------------
  # FigR Gene Regulatory Network Inference
  #------------------------------------------------------
  figR.d <- runFigRGRN(
    ATAC.se = ATAC.se,
    dorcTab = cisCorr.filt,
    genome = genome,
    dorcMat = dorcMat.s,
    rnaMat = RNAmat.s,
    dorcK = dorcK,
    nCores = nCores_grn
  )

  saveRDS(figR.d, file.path(outdir, paste0(prefix, "_FigR_GRN.rds")))

  #------------------------------------------------------
  # Visualization: TF–DORC Scatter Plot
  #------------------------------------------------------
  pdf(file.path(outdir, paste0(prefix, "_TF_DORC_scatter.pdf")), width = 4, height = 4)
  figR.d %>%
    ggplot(aes(Corr.log10P, Enrichment.log10P, color = Score)) +
    ggrastr::geom_point_rast(size = 0.01, shape = 16) +
    theme_classic() +
    scale_color_gradientn(
      colours = jdb_palette("solar_extra"),
      limits = c(-3, 3),
      oob = scales::squish,
      breaks = scales::breaks_pretty(n = 3)
    )
  dev.off()

  #------------------------------------------------------
  # Visualization: Driver Ranking
  #------------------------------------------------------
  pdf(file.path(outdir, paste0(prefix, "_TF_drivers.pdf")), width = 6, height = 4)
  rankDrivers(figR.d, rankBy = "meanScore", interactive = FALSE)
  dev.off()

  #------------------------------------------------------
  # Return all key results
  #------------------------------------------------------
  return(list(
    cisCorr = cisCorr,
    cisCorr_filt = cisCorr.filt,
    dorcGenes = dorcGenes,
    numDorcs = numDorcs,
    dorcMat_smooth = dorcMat.s,
    RNAmat_smooth = RNAmat.s,
    figR_GRN = figR.d
  ))
}

if (sys.nframe() == 0) {
  # CLI wrapper for Rscript usage with --key value args
  args <- commandArgs(trailingOnly = TRUE)
  arg_list <- list()
  if (length(args) > 0) {
    i <- 1
    while (i <= length(args)) {
      key <- args[[i]]
      if (startsWith(key, "--")) {
        name <- substring(key, 3)
        value <- if (i + 1 <= length(args) && !startsWith(args[[i + 1]], "--")) args[[i + 1]] else TRUE
        arg_list[[name]] <- value
        i <- i + if (isTRUE(value) && (i + 1 > length(args) || startsWith(args[[i + 1]], "--"))) 1 else 2
      } else {
        i <- i + 1
      }
    }
  }

  required <- c("atac_se_rds", "rna_mat_rds", "cellknn_rds")
  missing <- setdiff(required, names(arg_list))
  if (length(missing) > 0) {
    stop("Missing required arguments: ", paste(missing, collapse = ", "))
  }

  genome        <- if (!is.null(arg_list$genome)) arg_list$genome else "hg38"
  nCores_corr   <- if (!is.null(arg_list$nCores_corr)) as.integer(arg_list$nCores_corr) else 64
  nCores_smooth <- if (!is.null(arg_list$nCores_smooth)) as.integer(arg_list$nCores_smooth) else 24
  nCores_grn    <- if (!is.null(arg_list$nCores_grn)) as.integer(arg_list$nCores_grn) else 8
  dorc_cutoff   <- if (!is.null(arg_list$dorc_cutoff)) as.numeric(arg_list$dorc_cutoff) else 5
  dorc_labelTop <- if (!is.null(arg_list$dorc_labelTop)) as.integer(arg_list$dorc_labelTop) else 20
  dorcK         <- if (!is.null(arg_list$dorcK)) as.numeric(arg_list$dorcK) else 2
  prefix        <- if (!is.null(arg_list$prefix)) arg_list$prefix else "FigR"
  outdir        <- if (!is.null(arg_list$outdir)) arg_list$outdir else "./FigR_results"

  run_FigR_analysis(
    atac_se_rds   = arg_list$atac_se_rds,
    rna_mat_rds   = arg_list$rna_mat_rds,
    cellknn_rds   = arg_list$cellknn_rds,
    genome        = genome,
    nCores_corr   = nCores_corr,
    nCores_smooth = nCores_smooth,
    nCores_grn    = nCores_grn,
    dorc_cutoff   = dorc_cutoff,
    dorc_labelTop = dorc_labelTop,
    dorcK         = dorcK,
    prefix        = prefix,
    outdir        = outdir
  )
}
