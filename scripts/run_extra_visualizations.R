#!/usr/bin/env Rscript

# Generate supplementary scMINA visualizations from saved Seurat, differential
# expression, enrichment, and FigR results.

parse_cli <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Unexpected argument: ", key)
    name <- substring(key, 3L)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      out[[name]] <- TRUE
      i <- i + 1L
    } else {
      out[[name]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

require_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || !nzchar(as.character(value))) {
    stop("Missing required argument --", name)
  }
  value
}

require_file <- function(path, label) {
  if (!file.exists(path)) stop(label, " does not exist: ", path)
}

safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  gsub("^_+|_+$", "", x)
}

read_table_auto <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) return(readRDS(path))
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

save_coverage_plots <- function(obj, genes, outdir, prefix) {
  if (!requireNamespace("Signac", quietly = TRUE)) {
    stop("Coverage plots require the Signac package")
  }
  for (gene in genes) {
    plot <- Signac::CoveragePlot(
      object = obj,
      region = gene,
      peaks = FALSE,
      annotation = TRUE
    )
    ggplot2::ggsave(
      filename = file.path(outdir, paste0(prefix, "_", safe_name(gene), "_ATAC_coverage.pdf")),
      plot = plot,
      device = "pdf",
      width = 5,
      height = 10,
      units = "in",
      limitsize = FALSE
    )
  }
}

save_volcano <- function(genes, outdir, prefix, padj_cutoff, log2fc_cutoff) {
  required <- c("p_val_adj", "avg_log2FC")
  missing <- setdiff(required, colnames(genes))
  if (length(missing)) stop("DEG table is missing columns: ", paste(missing, collapse = ", "))

  genes$Significant <- "Non-significant"
  genes$Significant[genes$p_val_adj < padj_cutoff & genes$avg_log2FC >= log2fc_cutoff] <- "Up"
  genes$Significant[genes$p_val_adj < padj_cutoff & genes$avg_log2FC <= -log2fc_cutoff] <- "Down"
  genes$Significant <- factor(genes$Significant, levels = c("Down", "Non-significant", "Up"))
  genes$minus_log10_padj <- -log10(pmax(genes$p_val_adj, .Machine$double.xmin))

  significant <- genes[genes$Significant != "Non-significant", , drop = FALSE]
  write.csv(significant, file.path(outdir, paste0(prefix, "_DEGs_significant.csv")), row.names = FALSE)

  plot <- ggplot2::ggplot(genes, ggplot2::aes(avg_log2FC, minus_log10_padj, color = Significant)) +
    ggplot2::geom_point(size = 1, alpha = 0.8, na.rm = TRUE) +
    ggplot2::scale_color_manual(values = c(Down = "blue3", `Non-significant` = "grey70", Up = "red3")) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = 4, linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = 4, linewidth = 0.5) +
    ggplot2::labs(x = expression(Log[2] * FC), y = expression(-Log[10] * adjusted~P), title = prefix) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "right")

  ggplot2::ggsave(
    file.path(outdir, paste0(prefix, "_DEG_volcano.pdf")),
    plot = plot, device = "pdf", width = 6, height = 4, units = "in"
  )
}

save_enrichment_heatmap <- function(df, outdir, prefix, top_n, enrichment_direction) {
  required <- c("Term", "celltype", "Odds.Ratio", "Adjusted.P.value")
  missing <- setdiff(required, colnames(df))
  if (length(missing)) stop("Enrichment table is missing columns: ", paste(missing, collapse = ", "))
  if (!"Combined.Score" %in% colnames(df)) df$Combined.Score <- df$Odds.Ratio

  if ("direction" %in% colnames(df) && enrichment_direction != "all") {
    df <- df[tolower(df$direction) == tolower(enrichment_direction), , drop = FALSE]
  }
  if (!nrow(df)) stop("No enrichment rows remain after filtering")

  df <- df[order(df$celltype, -df$Combined.Score, df$Term), , drop = FALSE]
  top <- do.call(rbind, lapply(split(df, df$celltype), function(x) head(x, top_n)))
  selected_terms <- unique(top$Term)
  celltypes <- unique(df$celltype)
  mat <- matrix(0, nrow = length(selected_terms), ncol = length(celltypes),
                dimnames = list(selected_terms, celltypes))
  sig <- matrix(1, nrow = length(selected_terms), ncol = length(celltypes),
                dimnames = list(selected_terms, celltypes))

  for (i in seq_len(nrow(df))) {
    term <- df$Term[[i]]
    celltype <- df$celltype[[i]]
    if (term %in% selected_terms && celltype %in% celltypes) {
      mat[term, celltype] <- max(mat[term, celltype], df$Odds.Ratio[[i]], na.rm = TRUE)
      sig[term, celltype] <- min(sig[term, celltype], df$Adjusted.P.value[[i]], na.rm = TRUE)
    }
  }
  rownames(mat) <- sub(" \\(.*$", "", rownames(mat))
  rownames(sig) <- rownames(mat)

  finite_values <- mat[is.finite(mat)]
  upper <- if (length(finite_values)) max(finite_values) else 1
  midpoint <- upper / 2
  heatmap <- ComplexHeatmap::Heatmap(
    mat,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = circlize::colorRamp2(c(0, midpoint, upper), c("grey95", "deepskyblue", "dodgerblue4")),
    name = "Odds ratio",
    row_names_side = "left",
    row_names_max_width = ComplexHeatmap::max_text_width(rownames(mat), grid::gpar(fontsize = 10)),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (is.finite(sig[i, j]) && sig[i, j] < 0.05) {
        grid::grid.text("*", x, y, gp = grid::gpar(fontsize = 16, col = "white"))
      }
    }
  )
  grDevices::pdf(file.path(outdir, paste0(prefix, "_GO_enrichment_heatmap.pdf")), width = 13, height = 5)
  ComplexHeatmap::draw(heatmap)
  grDevices::dev.off()
}

save_network <- function(figr, outdir, prefix, score_cutoff, seed) {
  required <- c("Motif", "DORC", "Score", "Corr")
  missing <- setdiff(required, colnames(figr))
  if (length(missing)) stop("FigR table is missing columns: ", paste(missing, collapse = ", "))
  edges <- figr[is.finite(figr$Score) & abs(figr$Score) >= score_cutoff, , drop = FALSE]
  if (!nrow(edges)) stop("No FigR edges pass score cutoff ", score_cutoff)
  edges <- edges[order(edges$Motif, edges$DORC, edges$Score), , drop = FALSE]
  edges$Motif_node <- paste0(edges$Motif, ".")
  edges$weight <- scales::rescale(abs(edges$Score), to = c(0.5, 3))

  net <- network::network(edges[, c("Motif_node", "DORC")], directed = TRUE, matrix.type = "edgelist")
  network::set.edge.attribute(net, "weight", edges$weight)
  network::set.edge.attribute(net, "color", ifelse(edges$Corr > 0, "firebrick", "steelblue"))
  vertex_names <- network::get.vertex.attribute(net, "vertex.names")
  network::set.vertex.attribute(net, "class", ifelse(vertex_names %in% edges$DORC, "DORC", "Motif"))
  network::set.vertex.attribute(net, "vertex.names", sub("\\.$", "", vertex_names))

  set.seed(seed)
  plot <- GGally::ggnet2(
    net = net,
    label = TRUE,
    size = "class",
    color = "class",
    edge.size = "weight",
    edge.color = "color",
    edge.alpha = 1,
    alpha = 0.7,
    size.palette = c(DORC = 4, Motif = 1),
    label.size = 6,
    palette = c(Motif = "grey", DORC = "orange1"),
    legend.position = "none"
  ) + ggplot2::coord_equal() + ggplot2::theme(text = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(
    file.path(outdir, paste0(prefix, "_DORC_network.pdf")),
    plot = plot, device = "pdf", width = 15, height = 15,
    units = "in", limitsize = FALSE
  )
}

save_figr_heatmaps <- function(figr, genes, outdir, prefix, score_cutoff) {
  for (gene in genes) {
    path <- file.path(outdir, paste0(prefix, "_", safe_name(gene), "_FigR_heatmap.pdf"))
    grDevices::pdf(path, width = 7, height = 7)
    plot <- FigR::plotfigRHeatmap(
      figR.d = figr,
      score.cut = score_cutoff,
      DORCs = gene,
      column_names_gp = grid::gpar(fontsize = 6),
      show_row_dend = FALSE
    )
    if (inherits(plot, "Heatmap") || inherits(plot, "HeatmapList")) ComplexHeatmap::draw(plot) else print(plot)
    grDevices::dev.off()
  }
}

main <- function() {
  args <- parse_cli(commandArgs(trailingOnly = TRUE))
  outdir <- require_arg(args, "output_dir")
  prefix <- if (is.null(args$prefix)) "scMINA" else args$prefix
  seed <- if (is.null(args$seed)) 0L else as.integer(args$seed)
  genes <- trimws(strsplit(require_arg(args, "genes"), ",", fixed = TRUE)[[1]])
  plots <- if (is.null(args$plots)) c("coverage", "volcano", "enrichment", "network", "figr_heatmap") else
    trimws(strsplit(args$plots, ",", fixed = TRUE)[[1]])
  allowed <- c("coverage", "volcano", "enrichment", "network", "figr_heatmap")
  unknown <- setdiff(plots, allowed)
  if (length(unknown)) stop("Unknown plot type(s): ", paste(unknown, collapse = ", "))
  if (!length(genes) || any(!nzchar(genes))) stop("--genes must contain at least one gene")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  set.seed(seed)

  if ("coverage" %in% plots) {
    path <- require_arg(args, "seurat_obj")
    require_file(path, "Seurat object")
    save_coverage_plots(readRDS(path), genes, outdir, prefix)
  }
  if ("volcano" %in% plots) {
    path <- require_arg(args, "deg_rds")
    require_file(path, "DEG RDS")
    save_volcano(readRDS(path), outdir, prefix,
                 if (is.null(args$padj_cutoff)) 0.05 else as.numeric(args$padj_cutoff),
                 if (is.null(args$log2fc_cutoff)) 1 else as.numeric(args$log2fc_cutoff))
  }
  if ("enrichment" %in% plots) {
    path <- require_arg(args, "enrichment")
    require_file(path, "Enrichment table")
    save_enrichment_heatmap(read_table_auto(path), outdir, prefix,
                            if (is.null(args$top_terms_per_celltype)) 2L else as.integer(args$top_terms_per_celltype),
                            if (is.null(args$enrichment_direction)) "all" else args$enrichment_direction)
  }
  figr <- NULL
  if (any(c("network", "figr_heatmap") %in% plots)) {
    path <- require_arg(args, "figr_rds")
    require_file(path, "FigR result")
    figr <- readRDS(path)
  }
  if ("network" %in% plots) {
    save_network(figr, outdir, prefix,
                 if (is.null(args$network_score_cutoff)) 1.5 else as.numeric(args$network_score_cutoff), seed)
  }
  if ("figr_heatmap" %in% plots) {
    save_figr_heatmaps(figr, genes, outdir, prefix,
                       if (is.null(args$heatmap_score_cutoff)) 0.8 else as.numeric(args$heatmap_score_cutoff))
  }

  metadata <- data.frame(
    parameter = c("seed", "genes", "plots"),
    value = c(seed, paste(genes, collapse = ","), paste(plots, collapse = ","))
  )
  write.csv(metadata, file.path(outdir, paste0(prefix, "_visualization_parameters.csv")), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(outdir, paste0(prefix, "_sessionInfo.txt")))
  message("Saved visualization outputs to: ", normalizePath(outdir))
}

main()
