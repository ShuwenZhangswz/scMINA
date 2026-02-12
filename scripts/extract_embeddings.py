#!/usr/bin/env python3
"""
Extract e_df embeddings from saved scPair checkpoints + adata.
Use when checkpoints exist but embedding CSVs were not saved (e.g. from a prior run).
"""

import argparse
import json
from pathlib import Path

import pandas as pd
import scanpy as sc
import torch
from scpair import scPair_object, training_split


def parse_args():
    p = argparse.ArgumentParser(description="Extract embeddings from scPair checkpoints + adata")
    p.add_argument("--adata", required=True, help="Path to adata_paired.h5ad (must have scPair_split, var.modality)")
    p.add_argument("--checkpoint_dir", required=True, help="Dir with encoder_Gene Expression_to_Peaks.pt etc")
    p.add_argument("--output_dir", required=True, help="Output directory for embeddings CSVs")
    p.add_argument("--cov", default=None, help="Comma-separated covariate columns (must match training)")
    p.add_argument("--hidden_layer", default="800,30", help="Must match training")
    return p.parse_args()


def main():
    args = parse_args()
    adata = sc.read_h5ad(args.adata)
    ckpt_dir = Path(args.checkpoint_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    cov_list = [c.strip() for c in args.cov.split(",")] if args.cov else None
    hidden = [int(x) for x in args.hidden_layer.split(",")]

    obj = scPair_object(
        scobj=adata,
        cov=cov_list,
        modalities={"Gene Expression": "zinb", "Peaks": "ber"},
        sample_factor_rna=True,
        sample_factor_atac=False,
        infer_library_size_rna=False,
        infer_library_size_atac=True,
        hidden_layer=hidden,
        save_model=False,
    )
    obj.data_loader_builder()

    encoder_re = ckpt_dir / "encoder_Gene Expression_to_Peaks.pt"
    encoder_pe = ckpt_dir / "encoder_Peaks_to_Gene Expression.pt"
    if not encoder_re.exists():
        raise FileNotFoundError(f"Encoder not found: {encoder_re}")
    obj.encoder_dict = {
        "Gene Expression to Peaks": torch.load(encoder_re, map_location=obj.device),
        "Peaks to Gene Expression": torch.load(encoder_pe, map_location=obj.device),
    }
    for enc in obj.encoder_dict.values():
        enc.eval()

    e, e_df = obj.reference_embeddings()
    for key, df_emb in e_df.items():
        fname = f"embeddings_{key.replace(' ', '_')}.csv"
        df_emb.to_csv(out_dir / fname, sep="\t")
    print("Saved embeddings to", out_dir)


if __name__ == "__main__":
    main()
