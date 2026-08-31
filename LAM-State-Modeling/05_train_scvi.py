#!/usr/bin/env python3
"""Train the first dataset-only-batch scVI model from raw counts."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import scanpy as sc

from state_modeling_utils import PROJECT_ROOT, load_config, model_mask, validate_integer_counts, write_json


def choose_scvi_hvg(obj: ad.AnnData, n_hvg: int) -> np.ndarray:
    raw_view = obj.copy()
    raw_view.X = raw_view.layers["counts"].copy()
    n_hvg = min(max(1, int(n_hvg)), raw_view.n_vars)
    try:
        # Seurat v3 is the count-based HVG method.  Do not use the log-space
        # ``seurat`` flavor here: this branch is intentionally independent of
        # the NMF log-normalized matrix.
        sc.pp.highly_variable_genes(
            raw_view,
            n_top_genes=n_hvg,
            flavor="seurat_v3",
            layer="counts",
            batch_key="dataset",
            inplace=True,
        )
        mask = raw_view.var["highly_variable"].to_numpy().astype(bool)
    except Exception:
        mask = np.zeros(raw_view.n_vars, dtype=bool)
        mask[:n_hvg] = True
    if not mask.any():
        mask[:n_hvg] = True
    return mask


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    parser.add_argument("--max-epochs", type=int, default=None)
    parser.add_argument("--no-train", action="store_true", help="validate and prepare the scVI input without training")
    args = parser.parse_args()
    config = load_config(args.config)
    input_path = PROJECT_ROOT / config["outputs"]["prepared_h5ad"]
    if not input_path.exists():
        print(f"missing input: {input_path}; run scripts 02 and 03 first", file=sys.stderr)
        return 2

    prepared = ad.read_h5ad(input_path)
    include_normal = any(prepared.obs["analysis_role"].astype(str).eq("normal_reference"))
    mask = model_mask(prepared, include_normal=include_normal)
    if mask.sum() < 3:
        print("BLOCKED_INPUT: fewer than three candidate/boundary cells for scVI", file=sys.stderr)
        return 2
    obj = prepared[mask].copy()
    counts_audit = validate_integer_counts(obj)
    if not counts_audit["valid"] or "counts" not in obj.layers:
        print(f"BLOCKED_INPUT: invalid scVI raw counts: {counts_audit}", file=sys.stderr)
        return 2

    hvg_mask = choose_scvi_hvg(obj, int(config["preprocess"]["n_hvg"]))
    obj = obj[:, hvg_mask].copy()
    # Make the branch visibly raw-count based.  setup_anndata still points at
    # the layer, which is the source of truth for scVI.
    obj.X = obj.layers["counts"].copy()
    batch_key = config["scvi"]["batch_key"]
    if batch_key != "dataset":
        raise ValueError("scVI v1 must use batch_key=dataset")
    covariates = list(config["scvi"].get("categorical_covariate_keys", []))
    if covariates:
        raise ValueError(f"scVI v1 must not use categorical covariates: {covariates}")

    contract = {
        "layer": "counts",
        "batch_key": "dataset",
        "categorical_covariate_keys": [],
        "assay_retained_in_obs_only": True,
        "input_matrix": "raw integer-valued counts",
        "n_cells": int(obj.n_obs),
        "n_genes": int(obj.n_vars),
        "n_scvi_hvg": int(hvg_mask.sum()),
        "counts_audit": counts_audit,
        "included_roles": sorted(obj.obs["analysis_role"].astype(str).unique()),
    }
    write_json(PROJECT_ROOT / "results/stage1_6/scvi_input_contract.json", contract)
    if args.no_train:
        print("scvi_input_validated", contract)
        return 0

    import scvi
    import torch

    scvi.settings.seed = int(config["random_seed"])
    if hasattr(scvi.settings, "dl_pin_memory"):
        scvi.settings.dl_pin_memory = False
    scvi.model.SCVI.setup_anndata(obj, layer="counts", batch_key="dataset")
    model = scvi.model.SCVI(
        obj,
        n_latent=int(config["scvi"]["n_latent"]),
        n_layers=int(config["scvi"]["n_layers"]),
        n_hidden=int(config["scvi"]["n_hidden"]),
    )
    max_epochs = int(args.max_epochs or config["scvi"]["max_epochs"])
    requested_accelerator = str(config["scvi"].get("accelerator", "auto")).lower()
    cuda_available = bool(torch.cuda.is_available())
    if requested_accelerator == "auto":
        accelerator = "gpu" if cuda_available else "cpu"
    elif requested_accelerator in {"gpu", "cuda"}:
        if not cuda_available:
            print("CUDA requested but unavailable; falling back to CPU", file=sys.stderr)
        accelerator = "gpu" if cuda_available else "cpu"
    else:
        accelerator = "cpu"
    model.train(
        max_epochs=max_epochs,
        early_stopping=bool(config["scvi"]["early_stopping"]),
        accelerator=accelerator,
        devices=1,
    )
    model_dir = PROJECT_ROOT / config["outputs"]["scvi_model_dir"]
    model_dir.parent.mkdir(parents=True, exist_ok=True)
    model.save(model_dir, overwrite=True)

    latent = model.get_latent_representation()
    obj.obsm["X_scVI"] = latent.astype(np.float32)
    neighbors = min(int(config["preprocess"]["n_neighbors"]), obj.n_obs - 1)
    if neighbors >= 2:
        sc.pp.neighbors(obj, n_neighbors=neighbors, use_rep="X_scVI", random_state=int(config["random_seed"]))
        sc.tl.umap(obj, random_state=int(config["random_seed"]))
        sc.tl.leiden(obj, resolution=float(config["preprocess"]["leiden_resolution"]), key_added="leiden_scvi", random_state=int(config["random_seed"]), flavor="igraph", directed=False)
    else:
        obj.obs["leiden_scvi"] = "0"

    contract.update({
        "status": "trained",
        "max_epochs": max_epochs,
        "training_mode": "full" if max_epochs >= int(config["scvi"]["max_epochs"]) else "smoke",
        "model_dir": str(model_dir),
        "n_latent": int(latent.shape[1]),
        "runtime_accelerator": accelerator,
        "cuda_available": cuda_available,
        "n_scvi_clusters": int(obj.obs["leiden_scvi"].nunique()) if "leiden_scvi" in obj.obs else 1,
    })
    obj.uns["state_model_scvi"] = contract
    output = PROJECT_ROOT / config["outputs"]["scvi_h5ad"]
    output.parent.mkdir(parents=True, exist_ok=True)
    obj.write_h5ad(output, compression="gzip")
    write_json(PROJECT_ROOT / "results/stage1_6/scvi_training_manifest.json", contract)
    print(f"scvi={output} shape={obj.shape} latent={latent.shape}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
