# scMINA


## Installation

1. **Clone the repository:**
   ```bash
   git clone git@github.com:xyfqwlzoe/scMINA.git
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

All dependencies are managed through conda and specified in `environment.yml`:

### Python Dependencies
- Python>=3.10,<3.13
- scpair>=0.1.0
- scvi-tools>=1.0.0
- numpy, pandas, scipy, scikit-learn
- matplotlib, seaborn, scanpy, anndata
- ...

### R Dependencies
- R>=4.2.2
- Seurat
- ...

### Installation
All dependencies are automatically installed when you run:
```bash
./setup.sh
```

## Nextflow Integration

scMINA is designed to work seamlessly with Nextflow workflows, supporting both Python and R processes in the same environment.

### Quick Start with Nextflow:

1. **Set up the conda environment:**
   ```bash
   ./setup.sh
   ```

2. **Run the example workflow:**
   ```bash
   nextflow run workflows/example_workflow.nf -profile local_activated
   ```

3. **Customize the workflow:**
   - Edit `workflows/example_workflow.nf` for your specific needs
   - Modify `nextflow.config` for cluster-specific settings
   - Add your own Python/R scripts to `scripts/`

### Environment Files:

- `environment.yml` - Main conda environment with Python + R
- `nextflow.config` - Nextflow configuration

### scPair Feature Attribution Pipeline

Modular pipeline for scPair training, clustering, and cluster-specific feature attribution:

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

See `../scpair_nextflow_pipeline_plan.md` for full documentation.

### Workflow Features:

- **Unified Environment**: Both Python and R in the same conda environment
- **Cluster Support**: Configured for SLURM, PBS, SGE, and local execution
- **Resource Management**: Automatic CPU/memory allocation
- **Error Handling**: Retry logic and error reporting
- **Resume Capability**: Continue failed workflows from where they stopped

## Development

### Adding new dependencies:

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

### Installing Bioconductor/CRAN R packages inside the env
Some Bioconductor packages are not reliably available as conda recipes at the exact versions in your session. After activating the env, install them in R directly:

### Verify versions:
```bash
R --version
python --version
```


## Authors

- Shuwen Zhang (zhang.shuwen@mayo.edu)
- Hongru Hu (hrhu@ucdavis.edu)
