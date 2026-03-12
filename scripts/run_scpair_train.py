#!/usr/bin/env python3
"""
scPair training with two input modes + checkpoint save + e_df embedding export.
Mode A: Load paired h5ad with adata.var['modality'].
Mode B: Build from rna_adata, atac_adata + train/val/test split files.
"""

import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd
import scanpy as sc
from scpair import scPair_object, merge_paired_data, training_split


def parse_args():
    p = argparse.ArgumentParser(description="Run scPair training and export embeddings")
    p.add_argument("--input_mode", required=True, choices=["h5ad", "separate"],
                   help="h5ad: single paired h5ad. separate: rna + atac + meta + split files")
    p.add_argument("--input_h5ad", help="Path to paired h5ad (mode=h5ad)")
    p.add_argument("--input_rna", help="Path to RNA h5ad or matrix (mode=separate)")
    p.add_argument("--input_atac", help="Path to ATAC h5ad or matrix (mode=separate)")
    p.add_argument("--input_meta", help="Path to metadata CSV, index=cell_id (mode=separate)")
    p.add_argument("--input_index_dir", help="Dir with train_id.txt, val_id.txt, test_id.txt (mode=separate)")
    p.add_argument("--cov", default=None, help="Comma-separated covariate columns for batch correction")
    p.add_argument("--output_dir", required=True, help="Output directory")
    p.add_argument("--save_adata", action="store_true", help="Save/copy adata_paired.h5ad to output")
    p.add_argument("--export_embeddings", action="store_true", default=True, help="Export e_df embeddings after training")
    p.add_argument("--no_export_embeddings", action="store_false", dest="export_embeddings", help="Skip embedding export (run scpair_inference.nf separately)")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--hidden_layer", default="800,30", help="Comma-separated, e.g. 800,30")
    p.add_argument("--max_epochs", type=int, default=1000)
    return p.parse_args()


def load_rna_atac(rna_path: str, atac_path: str, meta_path: str):
    """Load RNA and ATAC as AnnData with shared obs (meta)."""
    meta = pd.read_csv(meta_path, index_col=0, sep="\t")
    if not meta.index.duplicated().any():
        meta = meta.loc[meta.index.dropna()]
    else:
        meta = meta[~meta.index.duplicated(keep="first")]

    def _load_modality(path: str):
        if path.endswith(".h5ad"):
            adata = sc.read_h5ad(path)
        elif path.endswith(".csv"):
            X = pd.read_csv(path, index_col=0)
            adata = sc.AnnData(X=X)
        else:
            raise ValueError(f"Unsupported format: {path}. Use .h5ad or .csv")
        common = adata.obs.index.intersection(meta.index)
        adata = adata[common].copy()
        adata.obs = meta.loc[common].copy()
        return adata

    rna_adata = _load_modality(rna_path)
    atac_adata = _load_modality(atac_path)
    common = rna_adata.obs.index.intersection(atac_adata.obs.index)
    rna_adata = rna_adata[common]
    atac_adata = atac_adata[common]
    return rna_adata, atac_adata


def load_splits(index_dir: str):
    """Load train, val, test cell IDs from index_dir."""
    def _read_ids(fname: str):
        path = Path(index_dir) / fname
        if not path.exists():
            return []
        df = pd.read_table(path, header=None)
        return df.iloc[:, 0].tolist()

    train_id = _read_ids("train_id.txt")
    val_id = _read_ids("val_id.txt")
    test_id = _read_ids("test_id.txt")
    return train_id, val_id, test_id


def main():
    args = parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    checkpoints_dir = out_dir / "checkpoints"
    checkpoints_dir.mkdir(parents=True, exist_ok=True)

    if args.input_mode == "h5ad":
        if not args.input_h5ad or not Path(args.input_h5ad).exists():
            print("ERROR: --input_h5ad required and must exist", file=sys.stderr)
            sys.exit(1)
        adata_paired = sc.read_h5ad(args.input_h5ad)
        if "scPair_split" not in adata_paired.obs.columns:
            print("WARNING: adata has no scPair_split. Applying training_split with default fracs.", file=sys.stderr)
            adata_paired = training_split(adata_paired, fracs=[0.8, 0.1, 0.1], seed=args.seed)
    else:
        if not all([args.input_rna, args.input_atac, args.input_meta, args.input_index_dir]):
            print("ERROR: separate mode requires --input_rna, --input_atac, --input_meta, --input_index_dir", file=sys.stderr)
            sys.exit(1)
        rna_adata, atac_adata = load_rna_atac(args.input_rna, args.input_atac, args.input_meta)
        train_id, val_id, test_id = load_splits(args.input_index_dir)
        if not train_id:
            print("ERROR: No train IDs found in input_index_dir", file=sys.stderr)
            sys.exit(1)
        adata_paired = merge_paired_data([rna_adata, atac_adata])
        adata_paired = training_split(adata_paired, pre_split=[train_id, val_id, test_id])

    cov_list = [c.strip() for c in args.cov.split(",")] if args.cov else None
    for c in (cov_list or []):
        if c not in adata_paired.obs.columns:
            print(f"WARNING: covariate '{c}' not in adata.obs, skipping", file=sys.stderr)
            cov_list = [x for x in cov_list if x in adata_paired.obs.columns]
            if not cov_list:
                cov_list = None

    hidden = [int(x) for x in args.hidden_layer.split(",")]

    scpair_setup = scPair_object(
        scobj=adata_paired,
        cov=cov_list,
        modalities={"Gene Expression": "zinb", "Peaks": "ber"},
        sample_factor_rna=True,
        sample_factor_atac=False,
        infer_library_size_rna=False,
        infer_library_size_atac=True,
        batchnorm=True,
        layernorm=True,
        SEED=args.seed,
        hidden_layer=hidden,
        dropout_rate=0.1,
        learning_rate_prediction=1e-3,
        max_epochs=args.max_epochs,
        save_model=True,
        save_path=str(checkpoints_dir),
    )
    scpair_setup.run()

    cov_cols = list(scpair_setup.cov_dummy.columns) if scpair_setup.cov_dummy is not None else []
    with open(out_dir / "cov_columns.json", "w") as f:
        json.dump({"cov_columns": cov_cols, "n_batch_cols": len(cov_cols)}, f, indent=2)

    if args.export_embeddings:
        e, e_df = scpair_setup.reference_embeddings()
        for key, df_emb in e_df.items():
            fname = f"embeddings_{key.replace(' ', '_')}.csv"
            df_emb.to_csv(out_dir / fname, sep="\t")

    train_meta = scpair_setup.train_metadata
    train_meta.to_csv(out_dir / "train_metadata.csv", sep="\t")

    meta = {
        "hidden_layer": hidden,
        "seed": args.seed,
        "cov_columns": cov_cols,
        "use_covariates": cov_list is not None and len(cov_list) > 0,
    }
    with open(checkpoints_dir / "checkpoint_meta.json", "w") as f:
        json.dump(meta, f, indent=2)

    run_params = {
        "input_mode": args.input_mode, "cov": args.cov, "hidden_layer": args.hidden_layer,
        "max_epochs": args.max_epochs, "seed": args.seed, "export_embeddings": args.export_embeddings,
    }
    with open(out_dir / "run_params.json", "w") as f:
        json.dump(run_params, f, indent=2)

    if args.save_adata:
        adata_paired.write_h5ad(out_dir / "adata_paired.h5ad")

    print("Done. Checkpoints:", checkpoints_dir)
    if args.export_embeddings:
        print("Embeddings:", [f"embeddings_{k.replace(' ', '_')}.csv" for k in e_df.keys()])
    else:
        print("Embeddings: skipped (use scpair_inference.nf to generate)")


if __name__ == "__main__":
    main()
