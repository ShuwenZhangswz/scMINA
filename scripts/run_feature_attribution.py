#!/usr/bin/env python3
"""
Cluster-specific feature attribution using Integrated Gradients and/or PathExplainer.
Uses saved scPair RNA encoder (encoder_Gene Expression_to_Peaks.pt).
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import torch
from torch import nn
from tqdm import tqdm


def parse_args():
    p = argparse.ArgumentParser(description="Run cluster-specific feature attribution")
    p.add_argument("--adata", required=True, help="Path to adata (paired h5ad with modality)")
    p.add_argument("--checkpoint_dir", required=True, help="Dir with encoder_Gene Expression_to_Peaks.pt")
    p.add_argument("--cluster_source", default="leiden", choices=["leiden", "obs_column"],
                   help="leiden: use cluster_labels file; obs_column: use adata.obs[--cluster_column]")
    p.add_argument("--cluster_labels", help="Path to cluster_labels_resX.csv (required when cluster_source=leiden)")
    p.add_argument("--cluster_column", help="obs column for cluster/cell type (required when cluster_source=obs_column)")
    p.add_argument("--cov_columns", help="Path to cov_columns.json from scpair_train (required when training used covariates)")
    p.add_argument("--output_dir", required=True, help="Output directory")
    p.add_argument("--attribution_method", default="both", choices=["ig", "pathexplainer", "both"])
    p.add_argument("--baseline", default="both", choices=["zero", "nonzero", "both"])
    p.add_argument("--use_absolute", action="store_true", help="Take abs before aggregating over cells")
    p.add_argument("--rank_aggregation", default="none", choices=["mean", "mean_abs", "none"])
    p.add_argument("--output_ranked", action="store_true", help="Output per-cluster top-N gene list")
    p.add_argument("--top_n_genes", type=int, default=50)
    p.add_argument("--group_column", help="Optional: run attribution separately per group (case/control)")
    return p.parse_args()


def get_model(encoder_path: str, device):
    encoder = torch.load(encoder_path, map_location=device)
    encoder.eval()
    model = nn.Sequential(encoder.hidden_rep, encoder.extra_linear_rep)
    model.eval()
    return model.to(device)


def run_ig(model, rna_input, baseline, latent_nodes, device):
    from captum.attr import IntegratedGradients
    ig = IntegratedGradients(model)
    n_genes_batch = rna_input.size(1)
    df = torch.zeros(latent_nodes, n_genes_batch, device="cpu")
    for node_id in range(latent_nodes):
        for i in range(rna_input.size(0)):
            inp = rna_input[i].reshape(1, -1).to(device).requires_grad_()
            attr = ig.attribute(inp, baseline, target=node_id, return_convergence_delta=False)
            df[node_id] += attr.cpu().detach().numpy().sum(0)
    return df.numpy()


def run_pathexplainer(model, rna_input, baseline, latent_nodes, device):
    import sys
    from pathlib import Path
    _scripts_dir = Path(__file__).resolve().parent
    if str(_scripts_dir) not in sys.path:
        sys.path.insert(0, str(_scripts_dir))
    from pathexplainer import PathExplainerTorch
    n_genes_batch = rna_input.size(1)
    df = torch.zeros(latent_nodes, n_genes_batch, device="cpu")
    baseline_grad = baseline.clone().requires_grad_()
    for node_id in range(latent_nodes):
        def wrapper(x):
            outs = model(x)
            return outs[:, node_id].reshape(-1, 1)
        explainer = PathExplainerTorch(wrapper)
        for i in range(rna_input.size(0)):
            inp = rna_input[i].reshape(1, -1).requires_grad_()
            attr = explainer.attributions(inp, baseline=baseline_grad, num_samples=1, use_expectation=False)
            df[node_id] += attr.cpu().detach().numpy().sum(0)
    return df.numpy()


def attribution_for_cluster(
    model, rna_train, batch_train, train_ids, cluster_df, col, cluster_id,
    genes, batch_cols, latent_nodes, device, method, baseline_opt, use_absolute,
    restrict_ids=None
):
    """Run attribution for one cluster. restrict_ids: if set, only consider these cells as 'other' for nonzero baseline."""
    cells_in_cluster = cluster_df[cluster_df[col].astype(str) == str(cluster_id)].index.tolist()
    mask = np.array([i for i, cid in enumerate(train_ids) if cid in cells_in_cluster])
    if len(mask) == 0:
        return None
    rna_sub = rna_train[mask]
    batch_sub = batch_train[mask] if batch_train is not None else None
    if batch_sub is not None:
        rna_input = torch.cat([torch.FloatTensor(rna_sub), torch.FloatTensor(batch_sub)], dim=1).to(device)
    else:
        rna_input = torch.FloatTensor(rna_sub).to(device)

    universe = set(restrict_ids) if restrict_ids is not None else set(train_ids)
    other_mask = np.array([i for i, cid in enumerate(train_ids) if cid in universe and cid not in cells_in_cluster])
    if baseline_opt == "zero":
        baseline = torch.zeros_like(rna_input[0]).reshape(1, -1).to(device)
    else:
        if len(other_mask) == 0:
            other_rna = rna_train.mean(0)
        else:
            other_rna = rna_train[other_mask].mean(0)
        if batch_sub is not None:
            baseline = torch.cat([
                torch.FloatTensor(other_rna).reshape(1, -1).to(device),
                torch.zeros(1, batch_sub.shape[1], device=device)
            ], dim=1)
        else:
            baseline = torch.FloatTensor(other_rna).reshape(1, -1).to(device)

    if method == "ig":
        arr = run_ig(model, rna_input, baseline, latent_nodes, device)
    else:
        arr = run_pathexplainer(model, rna_input, baseline, latent_nodes, device)

    if use_absolute:
        arr = np.abs(arr)
    cols = list(genes) + list(batch_cols)
    df = pd.DataFrame(arr, columns=cols)
    return df


def save_ranked(df, genes, batch_cols, rank_agg, top_n, out_path):
    """Extract gene cols, aggregate over latent nodes, save top N."""
    gene_cols = [c for c in df.columns if c not in batch_cols]
    df_rna = df[gene_cols]
    if rank_agg == "mean":
        scores = df_rna.mean(0)
    else:
        scores = df_rna.abs().mean(0)
    ranked = scores.sort_values(ascending=False).head(top_n)
    out_df = pd.DataFrame({"gene": ranked.index, "score": ranked.values, "rank": range(1, len(ranked) + 1)})
    out_df.to_csv(out_path, sep="\t", index=False)


def main():
    args = parse_args()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    adata = sc.read_h5ad(args.adata)
    ckpt_dir = Path(args.checkpoint_dir)
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    encoder_path = ckpt_dir / "encoder_Gene Expression_to_Peaks.pt"
    if not encoder_path.exists():
        raise FileNotFoundError(f"Encoder not found: {encoder_path}")
    model = get_model(str(encoder_path), device)

    rna_mask = adata.var["modality"] == "Gene Expression"
    rna_data = adata[:, rna_mask]
    if hasattr(rna_data.X, "toarray"):
        rna_mat = rna_data.X.toarray()
    else:
        rna_mat = np.array(rna_data.X)
    genes = rna_data.var_names.tolist()
    train_mask = adata.obs["scPair_split"] == "train"
    train_ids = adata.obs.index[train_mask].tolist()
    rna_train = rna_mat[train_mask]

    batch_cols = []
    cov_path = args.cov_columns
    if not cov_path and (ckpt_dir.parent / "cov_columns.json").exists():
        cov_path = str(ckpt_dir.parent / "cov_columns.json")
    if cov_path and Path(cov_path).exists():
        with open(cov_path) as f:
            meta = json.load(f)
            batch_cols = meta.get("cov_columns", [])
    batch_train = None
    if batch_cols and all(c in adata.obs.columns for c in batch_cols):
        batch_df = adata.obs[batch_cols].loc[train_ids]
        if np.isin(batch_df.values, [0, 1]).all():
            batch_train = batch_df.values.astype(np.float32)
        else:
            batch_train = pd.get_dummies(batch_df).values
    n_genes = len(genes)
    n_batch = batch_train.shape[1] if batch_train is not None else 0
    expected_in = n_genes + n_batch
    if hasattr(model[0], "in_features") and model[0].in_features != expected_in:
        raise ValueError(f"Encoder expects {model[0].in_features} features but rna+batch has {expected_in} (genes={n_genes}, batch={n_batch}). Pass --cov_columns from training if you used covariates.")

    if args.cluster_source == "leiden":
        if not args.cluster_labels or not Path(args.cluster_labels).exists():
            raise FileNotFoundError("--cluster_labels required when cluster_source=leiden")
        cluster_df = pd.read_csv(args.cluster_labels, index_col=0, sep="\t")
        col = cluster_df.columns[0]
    else:
        if not args.cluster_column or args.cluster_column not in adata.obs.columns:
            raise ValueError("--cluster_column required and must exist in adata.obs when cluster_source=obs_column")
        col = args.cluster_column
        cluster_df = adata.obs[[col]].copy()
    train_clusters = cluster_df.loc[cluster_df.index.isin(train_ids)]
    cluster_ids = train_clusters[col].astype(str).unique().tolist()

    latent_nodes = model[-1].out_features

    groups = [None]
    if args.group_column and args.group_column in adata.obs.columns:
        groups = adata.obs[args.group_column].dropna().unique().tolist()

    method_list = ["ig", "pathexplainer"] if args.attribution_method == "both" else [args.attribution_method]
    baseline_list = ["zero", "nonzero"] if args.baseline == "both" else [args.baseline]

    for grp in groups:
        if grp is not None:
            grp_sel = adata.obs.loc[train_ids, args.group_column] == grp
            sub_train_ids = [train_ids[i] for i in range(len(train_ids)) if grp_sel[i]]
            sub_cluster_df = cluster_df.loc[cluster_df.index.isin(sub_train_ids)]
            cluster_iter = sub_cluster_df[col].astype(str).unique().tolist()
            grp_suffix = f"_{grp}"
            restrict_ids = sub_train_ids
        else:
            sub_cluster_df = cluster_df.loc[cluster_df.index.isin(train_ids)]
            cluster_iter = cluster_ids
            grp_suffix = ""
            restrict_ids = None

        for ct in tqdm(cluster_iter, desc="Clusters"):
            for method in method_list:
                for bl in baseline_list:
                    df = attribution_for_cluster(
                        model, rna_train, batch_train, train_ids, sub_cluster_df, col, ct,
                        genes, batch_cols, latent_nodes, device,
                        method, bl, args.use_absolute, restrict_ids=restrict_ids
                    )
                    if df is None:
                        continue
                    bl_suf = "_nonzero_baseline" if bl == "nonzero" else ""
                    abs_suf = "_abs" if args.use_absolute else ""
                    name = f"rna_{method}_{ct}{bl_suf}{abs_suf}{grp_suffix}.csv"
                    df.to_csv(out_dir / name, sep="\t")
                    if args.output_ranked and args.rank_aggregation != "none":
                        save_ranked(
                            df, genes, batch_cols, args.rank_aggregation,
                            args.top_n_genes, out_dir / f"rna_{method}_{ct}{bl_suf}_top_genes{grp_suffix}.csv"
                        )
    print("Saved attribution results to", out_dir)


if __name__ == "__main__":
    main()
