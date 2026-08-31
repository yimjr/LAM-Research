#!/usr/bin/env python3
"""Run independent PCA/Leiden and State Modeling NMF baselines."""

from __future__ import annotations

import argparse
import gc
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.decomposition import NMF

from state_modeling_utils import PROJECT_ROOT, load_config, recreate_log_normalized_x, write_json


def subset_indices(mask: np.ndarray, max_cells: int, seed: int) -> np.ndarray:
    rows = np.flatnonzero(mask)
    if len(rows) <= max_cells:
        return rows
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(rows, size=max_cells, replace=False))


def choose_hvg(obj: ad.AnnData, n_hvg: int) -> np.ndarray:
    n_hvg = max(1, min(int(n_hvg), obj.n_vars))
    try:
        sc.pp.highly_variable_genes(obj, n_top_genes=n_hvg, flavor="seurat", inplace=True)
        mask = obj.var["highly_variable"].to_numpy().astype(bool)
    except Exception:
        x = obj.X
        if sparse.issparse(x):
            mean = np.asarray(x.mean(axis=0)).ravel()
            mean_sq = np.asarray(x.multiply(x).mean(axis=0)).ravel()
        else:
            dense = np.asarray(x, dtype=np.float32)
            mean = dense.mean(axis=0)
            mean_sq = (dense * dense).mean(axis=0)
        variance = np.maximum(mean_sq - mean * mean, 0.0)
        mask = np.zeros(obj.n_vars, dtype=bool)
        mask[np.argsort(variance)[-n_hvg:]] = True
        obj.var["highly_variable"] = mask
    if not mask.any():
        mask[:] = True
    return mask


def variance_for_features(x) -> np.ndarray:
    if sparse.issparse(x):
        mean = np.asarray(x.mean(axis=0)).ravel()
        mean_sq = np.asarray(x.multiply(x).mean(axis=0)).ravel()
    else:
        dense = np.asarray(x, dtype=np.float32)
        mean = dense.mean(axis=0)
        mean_sq = (dense * dense).mean(axis=0)
    return np.maximum(mean_sq - mean * mean, 0.0)


def run_pca_baseline(prepared: ad.AnnData, config: dict) -> tuple[ad.AnnData, dict]:
    roles = prepared.obs["analysis_role"].astype(str)
    mask = roles.isin({"primary_candidate", "boundary", "normal_reference"}).to_numpy()
    obj = prepared[mask].copy()
    if obj.n_obs < 3 or obj.n_vars < 2:
        obj.uns["state_model_baseline"] = {"status": "insufficient_cells_or_genes"}
        return obj, {"status": "insufficient_cells_or_genes", "n_cells": int(obj.n_obs)}
    hvg_mask = choose_hvg(obj, int(config["preprocess"]["n_hvg"]))
    obj = obj[:, hvg_mask].copy()
    n_pcs = min(int(config["preprocess"]["n_pcs"]), obj.n_obs - 1, obj.n_vars - 1)
    sc.pp.scale(obj, zero_center=True, max_value=10)
    sc.tl.pca(obj, n_comps=max(1, n_pcs), svd_solver="arpack", random_state=int(config["random_seed"]))
    neighbors = min(int(config["preprocess"]["n_neighbors"]), obj.n_obs - 1)
    sc.pp.neighbors(obj, n_neighbors=max(2, neighbors), n_pcs=max(1, n_pcs), random_state=int(config["random_seed"]))
    sc.tl.umap(obj, random_state=int(config["random_seed"]))
    sc.tl.leiden(obj, resolution=float(config["preprocess"]["leiden_resolution"]), key_added="leiden_baseline", random_state=int(config["random_seed"]), flavor="igraph", directed=False)
    obj.uns["state_model_baseline"] = {
        "status": "ready",
        "input_roles": ["primary_candidate", "boundary", "normal_reference"],
        "hvg": int(hvg_mask.sum()),
        "n_pcs": int(n_pcs),
        "n_neighbors": int(neighbors),
        "scaled_for_pca_only": True,
    }
    summary = obj.obs.groupby("leiden_baseline", observed=True).agg(
        cells=("leiden_baseline", "size"),
        patients=("patient_id", "nunique"),
        datasets=("dataset", "nunique"),
        high_confidence=("lam_candidate", "sum"),
    ).reset_index()
    summary.to_csv(PROJECT_ROOT / "results/stage1_6/baseline_cluster_summary.csv", index=False)
    return obj, obj.uns["state_model_baseline"]


def run_nmf(prepared: ad.AnnData, config: dict) -> dict:
    roles = prepared.obs["analysis_role"].astype(str)
    lam_mask = roles.isin({"primary_candidate", "boundary"}).to_numpy()
    rows = subset_indices(lam_mask, int(config["nmf"]["max_cells"]), int(config["random_seed"]))
    if len(rows) < 2:
        return {"status": "insufficient_cells", "n_cells": int(len(rows))}
    obj = prepared[rows].copy()
    # Recreate the NMF branch from counts even if the prepared object's X was
    # later altered by another analysis branch.
    recreate_log_normalized_x(obj, float(config["preprocess"]["target_sum"]))
    hvg_mask = choose_hvg(obj, int(config["preprocess"]["n_hvg"]))
    hvg_indices = np.flatnonzero(hvg_mask)
    variance = variance_for_features(obj.X[:, hvg_indices])
    feature_count = min(int(config["nmf"]["max_features"]), len(hvg_indices))
    selected_local = np.argsort(variance)[-feature_count:]
    selected_indices = hvg_indices[selected_local]
    selected_indices = np.sort(selected_indices)
    matrix = obj.X[:, selected_indices]
    if sparse.issparse(matrix):
        matrix = matrix.tocsr().astype(np.float32)
    else:
        matrix = np.asarray(matrix, dtype=np.float32)
    if float(matrix.min()) < 0:
        raise ValueError("NMF input became negative; NMF must receive non-negative log1p data")
    model = NMF(
        n_components=min(int(config["nmf"]["n_components"]), len(selected_indices), len(rows)),
        init="nndsvda",
        random_state=int(config["random_seed"]),
        max_iter=int(config["nmf"]["max_iter"]),
    )
    scores = model.fit_transform(matrix)
    names = prepared.var_names.astype(str).to_numpy()[selected_indices]
    score_columns = [f"nmf_{i + 1}" for i in range(scores.shape[1])]
    score_frame = pd.DataFrame(scores, columns=score_columns)
    score_frame.insert(0, "analysis_cell_id", obj.obs_names.astype(str).to_numpy())
    score_frame.insert(1, "cell_id", obj.obs["source_cell_id"].astype(str).to_numpy())
    for col in ["dataset", "patient_id", "donor_id", "assay", "lam_candidate", "boundary"]:
        score_frame[col] = obj.obs[col].to_numpy()
    score_frame.to_csv(PROJECT_ROOT / "results/stage1_6/nmf_cell_scores.csv", index=False)

    top_rows = []
    for component, weights in enumerate(model.components_):
        order = np.argsort(weights)[::-1][: int(config["nmf"]["top_genes"])]
        for rank, pos in enumerate(order, start=1):
            top_rows.append({"component": f"nmf_{component + 1}", "gene": names[pos], "weight": float(weights[pos]), "rank": rank})
    pd.DataFrame(top_rows).to_csv(PROJECT_ROOT / "results/stage1_6/nmf_top_genes.csv", index=False)
    manifest = {
        "status": "ready",
        "input_contract": "counts -> normalize_total -> log1p -> HVG -> NMF",
        "n_cells": int(len(rows)),
        "n_hvg": int(len(hvg_indices)),
        "n_features": int(len(selected_indices)),
        "n_components": int(scores.shape[1]),
        "scaled_or_pca_input": False,
    }
    write_json(PROJECT_ROOT / "results/stage1_6/nmf_manifest.json", manifest)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    input_path = PROJECT_ROOT / config["outputs"]["prepared_h5ad"]
    if not input_path.exists():
        print(f"missing input: {input_path}; run scripts 02 and 03 first", file=sys.stderr)
        return 2
    prepared = ad.read_h5ad(input_path)
    # Run the matrix-heavy NMF branch first so its temporary object is
    # released before PCA creates a dense scaled HVG matrix.
    nmf_manifest = run_nmf(prepared, config)
    gc.collect()
    baseline, baseline_manifest = run_pca_baseline(prepared, config)
    output = PROJECT_ROOT / config["outputs"]["baseline_h5ad"]
    output.parent.mkdir(parents=True, exist_ok=True)
    baseline.write_h5ad(output, compression="gzip")
    write_json(PROJECT_ROOT / "results/stage1_6/baseline_manifest.json", {"pca": baseline_manifest, "nmf": nmf_manifest, "output": str(output)})
    print(f"baseline={output} pca={baseline.shape} nmf={nmf_manifest.get('status')}")
    return 0 if nmf_manifest.get("status") == "ready" else 2


if __name__ == "__main__":
    sys.exit(main())
