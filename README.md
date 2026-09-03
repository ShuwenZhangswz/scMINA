[![DOI](https://zenodo.org/badge/1052473576.svg)](https://doi.org/10.5281/zenodo.19582767)


# scMINA: single cell Multimodal Integration and Network Attribution

scMINA is a Nextflow-based framework for single-cell multimodal (RNA + ATAC) analysis. It combines scPair deep generative modeling as an option in Python with Seurat and FigR workflows in R for gene regulatory network (GRN) inference, and explainable downstream analysis.


## Installation

1. **Clone the repository:**
   ```bash
   git clone git@github.com:ShuwenZhangswz/scMINA.git
   cd scMINA
   ```

2. **Set up the conda environment:**
   ```bash
   ./setup.sh
   ```

3. **Activate the environment:**
   ```bash
   conda activate scmina
   ```


## Dependencies

All dependencies are managed through conda and specified in `environment.yml`.

### Python Dependencies
- Python>=3.10,<3.13
- scpair>=0.1.0
- scvi-tools>=1.0.0
- numpy, pandas, scipy, scikit-learn
- matplotlib, seaborn, scanpy, anndata
- See `environment.yml` for the full and up-to-date list.

### R / Bioconductor Dependencies
- R>=4.2.2
- Seurat and SeuratObject
- FigR and supporting Bioconductor packages for ATAC and GRN analysis
- See `environment.yml` for the full list of R/Bioconductor packages.

### Environment files

- `environment.yml` – main conda environment with Python + R
- `nextflow.config` – Nextflow configuration and resource labels

### Extra figure visualizations (R)

The `workflows/extra_visualizations.nf` workflow runs `scripts/run_extra_visualizations.R` to generate publication-oriented supplementary figures from saved scMINA results:

- ATAC coverage plots for one or more genes
- a DEG volcano plot and significant-DEG table
- a GO enrichment heatmap across cell types
- a FigR TF–DORC network
- a FigR heatmap for each selected gene

The workflow also saves its parameters and R `sessionInfo()` alongside the figures. Network layout is initialized with `--seed`; use the same seed and software environment to reproduce a figure.

#### Required inputs

| Parameter | Description |
| --- | --- |
| `--seurat_obj` | Seurat multiome object containing the assay required by `Signac::CoveragePlot` |
| `--deg_rds` | RDS containing differential-expression results |
| `--enrichment` | GO enrichment results as CSV or RDS |
| `--figr_rds` | FigR GRN result table saved as RDS |
| `--genes` | One or more comma-separated gene/DORC names |

The DEG table must contain `p_val_adj` and `avg_log2FC`. The enrichment table must contain `Term`, `celltype`, `Odds.Ratio`, and `Adjusted.P.value`; `Combined.Score` and `direction` are optional. The FigR table must contain `Motif`, `DORC`, `Score`, and `Corr`.

#### Example

```bash
nextflow run workflows/extra_visualizations.nf \
  --seurat_obj /path/to/multiome_seurat.rds \
  --deg_rds /path/to/DEGs.rds \
  --enrichment /path/to/GO_enrichment.csv \
  --figr_rds /path/to/Sample1_FigR_GRN.rds \
  --genes "GENE1,GENE2" \
  --prefix Sample1 \
  --seed 42 \
  -profile local_activated
```

By default, all five plot types are generated. Select a subset with a comma-separated list:

```bash
--plots "volcano,enrichment,network"
```

Available plot names are `coverage`, `volcano`, `enrichment`, `network`, and `figr_heatmap`.

#### Visualization parameters

| Parameter | Default | Description |
| --- | ---: | --- |
| `--seed` | `0` | Seed used for reproducible network layout |
| `--padj_cutoff` | `0.05` | Adjusted P-value cutoff for significant DEGs and the volcano reference line |
| `--log2fc_cutoff` | `1.0` | Absolute log2 fold-change cutoff for significant DEGs |
| `--top_terms_per_celltype` | `2` | Number of top enrichment terms selected per cell type |
| `--enrichment_direction` | `all` | Direction to retain when the enrichment table contains a `direction` column |
| `--network_score_cutoff` | `1.5` | Minimum absolute FigR score for network edges |
| `--heatmap_score_cutoff` | `0.8` | FigR heatmap score cutoff |

Outputs are published under `results/extra_visualizations` by default. They include PDF figures, the significant-DEG CSV, `<prefix>_visualization_parameters.csv`, and `<prefix>_sessionInfo.txt`.

## Reproducibility and random seeds

scMINA exposes a shared `--seed` parameter for stochastic analysis steps. The
default seed is `0`, but users should provide the seed explicitly when running
an analysis intended for comparison, publication, or reproduction.

The seed controls the following steps:

| Workflow | Seeded operations |
| --- | --- |
| `scpair_pipeline.nf` / `scpair_train.nf` | Train/validation/test splitting and scPair model training |
| `scpair_cluster.nf` | Scanpy neighbor graph construction and Leiden clustering |
| `integrate_scpair_multiome.nf` | Seurat clustering and UMAP |
| `figr_pipeline.nf` | FigR background sampling and other stochastic R operations |
| `feature_attribution.nf` | NumPy and PyTorch operations used during feature attribution |

### Using a seed

Pass the seed as a Nextflow parameter:

```bash
nextflow run workflows/scpair_pipeline.nf \
  --input_h5ad /path/to/paired.h5ad \
  --seed 42 \
  -profile local_activated
```

For clustering only:

```bash
nextflow run workflows/scpair_cluster.nf \
  --embeddings_dir /path/to/embeddings \
  --seed 42 \
  -profile local_activated
```

For Seurat integration:

```bash
nextflow run workflows/integrate_scpair_multiome.nf \
  --seurat_obj_path /path/to/multiome_seurat.rds \
  --scpair_csv /path/to/scpair_embeddings.csv \
  --metadata_csv /path/to/metadata.csv \
  --resolution 0.9 \
  --seed 42 \
  --prefix Sample1 \
  -profile local_activated
```

For FigR analysis:

```bash
nextflow run workflows/figr_pipeline.nf \
  --atac_mtx /path/to/ATACmat.mtx \
  --rna_mtx /path/to/RNAmat.mtx \
  --metadata_csv /path/to/metadata.csv \
  --genes_csv /path/to/genes.csv \
  --peaks_csv /path/to/peaks.csv \
  --seurat_scpair_rds /path/to/Sample1_scPair_final_res0.9.rds \
  --genome hg38 \
  --seed 42 \
  -profile local_activated
```

Use the same seed for every workflow stage when reproducing a complete
analysis:

```bash
SEED=42
```

Then pass it to each Nextflow command:

```bash
nextflow run workflows/scpair_pipeline.nf \
  --input_h5ad /path/to/paired.h5ad \
  --seed "${SEED}" \
  -profile local_activated

nextflow run workflows/integrate_scpair_multiome.nf \
  --seurat_obj_path /path/to/multiome_seurat.rds \
  --scpair_csv /path/to/scpair_embeddings.csv \
  --metadata_csv /path/to/metadata.csv \
  --seed "${SEED}" \
  -profile local_activated

nextflow run workflows/figr_pipeline.nf \
  --atac_mtx /path/to/ATACmat.mtx \
  --rna_mtx /path/to/RNAmat.mtx \
  --metadata_csv /path/to/metadata.csv \
  --genes_csv /path/to/genes.csv \
  --peaks_csv /path/to/peaks.csv \
  --seurat_scpair_rds /path/to/seurat_scpair.rds \
  --seed "${SEED}" \
  -profile local_activated
```

## Nextflow integration and overall workflow

scMINA is designed to work seamlessly with Nextflow, orchestrating both Python (scPair) and R (Seurat + FigR) processes in the same environment.

### Quick start with Nextflow

1. **Set up the conda environment:**
   ```bash
   ./setup.sh
   ```

2. **Run the example workflow:**
   ```bash
   nextflow run workflows/example_workflow.nf -profile local_activated
   ```


### Multimodal integration and downstream analysis workflows using scPair and FigR

![Overview of data processing, multimodal integration, and downstream analysis workflows](workflow_fig.png)

The figure above and the schematic below summarize how the main workflows connect across data processing, multimodal integration, and downstream GRN analysis:

```text
scpair_train.nf  -->  scpair_cluster.nf
      |                    |
      v                    v
  scPair embeddings   cluster labels
          |
          v
integrate_scpair_multiome.nf
          |
          v
  Seurat object with scPair embeddings
          |
          v
      figr_pipeline.nf
          |
          v
   FigR GRN results and plots
```


## Workflows

### scPair training, clustering, and feature attribution (Python)

The scPair workflows handle multimodal integration in Python and provide clustering and feature attribution on learned embeddings.

```bash
# Full scPair pipeline: train -> embeddings -> clustering
nextflow run workflows/scpair_pipeline.nf --input_h5ad /path/to/paired.h5ad

# Or with separate inputs (rna + atac + meta + splits):
nextflow run workflows/scpair_pipeline.nf --input_mode separate \
  --input_rna /path/to/rna.h5ad --input_atac /path/to/atac.h5ad \
  --input_meta /path/to/meta.csv --input_index_dir /path/to/splits/

# Clustering only (using existing embeddings from a prior run):
nextflow run workflows/scpair_pipeline.nf --run clustering_only --embeddings_dir /path/to/embeddings/

# Feature attribution (after choosing resolution):
nextflow run workflows/feature_attribution.nf \
  --adata /path/to/adata.h5ad \
  --checkpoint_dir /path/to/checkpoints/ \
  --cluster_labels /path/to/cluster_labels_res0.9.csv \
  --attribution_method both --baseline both --output_ranked --top_n_genes 50
```

See `../scpair_nextflow_pipeline_plan.md` for more details on the scPair pipeline.


### Integration of scPair embeddings into Seurat (R)

The `workflows/integrate_scpair_multiome.nf` workflow runs the R helper `scripts/integrate_scPair_multiome.R` to:

- Load a multiome Seurat object and scPair embeddings exported from the Python pipeline.
- Add scPair embeddings as a Seurat `DimReduc`.
- Perform graph construction, clustering, and UMAP.
- Save updated Seurat objects, markers, and plots.

Example:

```bash
nextflow run workflows/integrate_scpair_multiome.nf \
  --seurat_obj_path /path/to/multiome_seurat.rds \
  --scpair_csv /path/to/scpair_embeddings.csv \
  --metadata_csv /path/to/metadata.csv \
  --resolution 0.9 \
  --prefix Sample1
```


### FigR preprocessing and GRN analysis (R)

The `workflows/figr_pipeline.nf` workflow connects the Seurat object with scPair embeddings to FigR. It runs two R helpers in sequence:

1. `scripts/prep_FigR_inputs.R` – builds FigR-ready inputs from multiome matrices and the Seurat object:
   - ATAC `SummarizedExperiment` (`ATAC.se`)
   - normalized RNA matrix (`RNAmat`)
   - cell kNN index (`cellkNN`) derived from scPair embeddings
2. `scripts/run_FigR_analysis.R` – performs:
   - peak–gene correlation testing
   - DORC identification and visualization
   - DORC and RNA smoothing
   - GRN inference and TF driver ranking

Example:

```bash
nextflow run workflows/figr_pipeline.nf \
  --atac_mtx /path/to/ATACmat.mtx \
  --rna_mtx /path/to/RNAmat.mtx \
  --metadata_csv /path/to/metadata.csv \
  --genes_csv /path/to/genes.csv \
  --peaks_csv /path/to/peaks.csv \
  --seurat_scpair_rds /path/to/Sample1_scPair_final_res0.9.rds \
  --prefix Sample1 \
  --genome hg38
```


## System requirements

- **GPU:** `scpair_train`, `scpair_inference`, and `feature_attribution` require at least 1 GPU (PyTorch/CUDA). `scpair_cluster` is CPU-only. For SLURM, `nextflow.config` requests `--gres=gpu:1` for GPU processes.
- **Local runs:** Use `-profile local_activated`; ensure CUDA is available if running scPair steps locally.
- The same conda environment is used for all Python (scPair) and R (Seurat + FigR) steps; see `environment.yml` and `nextflow.config` for resource labels and profiles.

### Workflow features

- **Unified environment**: Both Python and R in the same conda environment.
- **Cluster support**: Configured for SLURM, PBS, SGE, and local execution.
- **Resource management**: Automatic CPU/memory allocation (GPU for scPair train/inference/attribution).
- **Error handling**: Retry logic and error reporting.
- **Resume capability**: Continue failed workflows from where they stopped.


## Development

### Adding new dependencies

```bash
# Add to environment.yml, then update
conda env update -f environment.yml
```

```bash
conda activate scmina
```

Notes:
- The conda env pins R to 4.2.2 to match the testing cluster.
- If a conda solve fails due to strict R package pins, prefer installing that package via BiocManager after env creation.

When adding new R/Bioconductor dependencies for analysis steps, prefer updating `environment.yml` first. For packages not available on conda, install them inside the environment using BiocManager or remotes.

### Installing Bioconductor/CRAN R packages inside the env

Some Bioconductor packages are not reliably available as conda recipes at specific versions. After activating the env, install them in R directly as needed.

### Verify versions

```bash
R --version
python --version
```


## Authors

- Shuwen Zhang (zhang.shuwen@mayo.edu)
- Hongru Hu (hrhu@ucdavis.edu)
