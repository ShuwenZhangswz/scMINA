#!/usr/bin/env python3
"""
Cluster on embeddings: RNA+Peaks (concat), RNA only, or Peaks only.
"""

import argparse
import json
from pathlib import Path

import pandas as pd
import scanpy as sc
import numpy as np


def parse_args():
    p = argparse.ArgumentParser(description="Cluster on embeddings")
    p.add_argument("--embeddings_dir", required=True, help="Dir with embeddings_Gene_Expression_*.csv, embeddings_Peaks_*.csv")
    p.add_argument("--output_dir", required=True, help="Output directory")
    p.add_argument("--embedding_mode", default="all", choices=["all", "rna_only", "peaks_only"],
                   help="all: concat RNA+Peaks; rna_only: Gene Expression only; peaks_only: Peaks only")
    p.add_argument("--resolutions", default="0.5,0.9,1.2,1.5", help="Comma-separated resolution values for Leiden")
    p.add_argument("--save_concat", action="store_true", help="Save embedding matrix used for clustering")
    return p.parse_args()


def main():
    args = parse_args()
    emb_dir = Path(args.embeddings_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    resolutions = [float(x.strip()) for x in args.resolutions.split(",")]

    rna_keys = ["Gene_Expression_train", "Gene_Expression_val", "Gene_Expression_test"]
    peak_keys = ["Peaks_train", "Peaks_val", "Peaks_test"]
    rna_dfs = {}
    peak_dfs = {}
    for k in rna_keys:
        f = emb_dir / f"embeddings_{k}.csv"
        if f.exists():
            rna_dfs[k] = pd.read_csv(f, index_col=0, sep="\t")
    for k in peak_keys:
        f = emb_dir / f"embeddings_{k}.csv"
        if f.exists():
            peak_dfs[k] = pd.read_csv(f, index_col=0, sep="\t")

    mode = args.embedding_mode
    if mode == "all" and (not rna_dfs or not peak_dfs):
        raise FileNotFoundError(f"Need both RNA and Peaks embeddings for mode=all. Found RNA: {list(rna_dfs)}, Peaks: {list(peak_dfs)}")
    if mode == "rna_only" and not rna_dfs:
        raise FileNotFoundError(f"Need RNA embeddings for mode=rna_only. Found: {list(rna_dfs)}")
    if mode == "peaks_only" and not peak_dfs:
        raise FileNotFoundError(f"Need Peaks embeddings for mode=peaks_only. Found: {list(peak_dfs)}")

    all_cells = []
    rna_concat_list = []
    peak_concat_list = []
    for suf in ["train", "val", "test"]:
        rk = f"Gene_Expression_{suf}"
        pk = f"Peaks_{suf}"
        if mode == "rna_only":
            if rk not in rna_dfs:
                continue
            rdf = rna_dfs[rk]
            common = rdf.index
            for c in common:
                all_cells.append(c)
                rna_concat_list.append(rdf.loc[c].values)
            peak_concat_list = []  # not used
        elif mode == "peaks_only":
            if pk not in peak_dfs:
                continue
            pdf = peak_dfs[pk]
            common = pdf.index
            for c in common:
                all_cells.append(c)
                peak_concat_list.append(pdf.loc[c].values)
            rna_concat_list = []
        else:
            if rk not in rna_dfs or pk not in peak_dfs:
                continue
            rdf, pdf = rna_dfs[rk], peak_dfs[pk]
            common = rdf.index.intersection(pdf.index)
            for c in common:
                all_cells.append(c)
                rna_concat_list.append(rdf.loc[c].values)
                peak_concat_list.append(pdf.loc[c].values)

    if mode == "rna_only":
        concat_emb = np.array(rna_concat_list)
    elif mode == "peaks_only":
        concat_emb = np.array(peak_concat_list)
    else:
        rna_arr = np.array(rna_concat_list)
        peak_arr = np.array(peak_concat_list)
        concat_emb = np.hstack([rna_arr, peak_arr])
    cell_ids = pd.Index(all_cells)

    adata = sc.AnnData(X=concat_emb, obs=pd.DataFrame(index=cell_ids))
    sc.pp.neighbors(adata, use_rep="X", n_neighbors=15)

    for res in resolutions:
        sc.tl.leiden(adata, resolution=res, key_added=f"leiden_res{res}")

    for res in resolutions:
        col = f"leiden_res{res}"
        labels = pd.DataFrame({col: adata.obs[col].astype(str)}, index=adata.obs.index)
        labels.to_csv(out_dir / f"cluster_labels_res{res}.csv", sep="\t")

    if args.save_concat:
        concat_df = pd.DataFrame(concat_emb, index=cell_ids)
        concat_df.to_csv(out_dir / "concat_embeddings.csv", sep="\t")

    cluster_params = {
        "embedding_mode": mode,
        "resolutions": resolutions,
    }
    with open(out_dir / "clustering_params.json", "w") as f:
        json.dump(cluster_params, f, indent=2)

    print("Saved cluster labels to", out_dir)


if __name__ == "__main__":
    main()
