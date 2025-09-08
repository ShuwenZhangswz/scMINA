#!/usr/bin/env Rscript
# Example R script for scMINA data visualization

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(plotly)
    library(htmlwidgets)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Default values
input_file <- "processed_data.h5ad"
output_dir <- "."
log_file <- "r_log.txt"

# Parse arguments
for (i in 1:length(args)) {
    if (args[i] == "--input" && i < length(args)) {
        input_file <- args[i + 1]
    } else if (args[i] == "--output" && i < length(args)) {
        output_dir <- args[i + 1]
    } else if (args[i] == "--log" && i < length(args)) {
        log_file <- args[i + 1]
    }
}

# Setup logging
log_path <- file.path(dirname(log_file), basename(log_file))
dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)

cat("Starting scMINA data visualization...\n", file = log_path, append = FALSE)

# Load data
cat(paste("Loading data from", input_file, "...\n"), file = log_path, append = TRUE)

# Note: In a real workflow, you'd use reticulate to load h5ad files
# For this example, we'll create a mock Seurat object
set.seed(42)
n_cells <- 1000
n_genes <- 2000

# Create mock expression matrix
expr_matrix <- matrix(
    rpois(n_cells * n_genes, lambda = 2),
    nrow = n_genes,
    ncol = n_cells
)
rownames(expr_matrix) <- paste0("Gene_", 1:n_genes)
colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = expr_matrix, project = "scMINA")

cat(paste("Loaded", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n"), 
    file = log_path, append = TRUE)

# Basic preprocessing
cat("Performing basic preprocessing...\n", file = log_path, append = TRUE)

# Normalize data
seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)

# Find variable features
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

# Scale data
seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)

# PCA
cat("Performing PCA...\n", file = log_path, append = TRUE)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj), verbose = FALSE)


# Clustering
cat("Performing clustering...\n", file = log_path, append = TRUE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10, verbose = FALSE)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5, verbose = FALSE)

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Generate plots
cat("Generating visualizations...\n", file = log_path, append = TRUE)

# PCA plot
p1 <- DimPlot(seurat_obj, reduction = "pca", group.by = "seurat_clusters") +
    ggtitle("PCA Clustering") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave(file.path(output_dir, "pca_clusters.png"), p1, width = 8, height = 6, dpi = 300)


cat("Visualization completed successfully!\n", file = log_path, append = TRUE)
cat(paste("Plots saved to:", output_dir, "\n"), file = log_path, append = TRUE)

# Print summary
cat("\n=== Analysis Summary ===\n")
cat(paste("Number of cells:", ncol(seurat_obj), "\n"))
cat(paste("Number of genes:", nrow(seurat_obj), "\n"))
cat(paste("Variable features:", length(VariableFeatures(seurat_obj)), "\n"))
cat("=======================\n")
