#!/usr/bin/env python3
"""Step 12: auxiliary boundary and normal-reference validation."""

from __future__ import annotations

import argparse
import gc
import os
import shutil

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from state_modeling_utils import PROJECT_ROOT, load_config, write_json


def normal_validation(consensus: ad.AnnData, full_path, config: dict, out_path) -> tuple[pd.DataFrame, str]:
    full = ad.read_h5ad(full_path, backed="r")
    normal_mask = full.obs.get("analysis_role", pd.Series("", index=full.obs.index)).astype(str).eq("normal_reference").to_numpy()
    if not normal_mask.any():
        result = pd.DataFrame(columns=["normal_analysis_cell_id", "nearest_state", "mean_distance", "min_distance"])
        result.to_csv(out_path, index=False)
        return result, "not_available"
    full_index = pd.Index(full.obs.index.astype(str))
    positions = pd.Series(np.arange(len(full_index)), index=full_index)
    high_positions = positions.reindex(consensus.obs.index.astype(str)).dropna().astype(int).to_numpy()
    normal_positions = np.flatnonzero(normal_mask)
    high_latent = np.asarray(full.obsm["X_scVI"][high_positions], dtype=np.float32)
    normal_latent = np.asarray(full.obsm["X_scVI"][normal_positions], dtype=np.float32)
    states = consensus.obs["consensus_state"].astype(str).to_numpy()
    k = min(int(config["step12"].get("normal_k", 15)), len(states))
    nn = NearestNeighbors(n_neighbors=max(1, k), metric="euclidean").fit(high_latent)
    rows = []
    batch_size = int(config["step12"].get("query_batch_size", 512))
    for start in range(0, len(normal_positions), batch_size):
        stop = min(start + batch_size, len(normal_positions))
        distances, indices = nn.kneighbors(normal_latent[start:stop])
        for local in range(stop - start):
            neighbor_states = states[indices[local]]
            counts = pd.Series(neighbor_states).value_counts()
            rows.append({
                "normal_analysis_cell_id": str(full.obs.iloc[normal_positions[start + local]].get("analysis_cell_id", full.obs.index[normal_positions[start + local]])),
                "normal_patient_id": str(full.obs.iloc[normal_positions[start + local]].get("patient_id", "")),
                "normal_dataset": str(full.obs.iloc[normal_positions[start + local]].get("dataset", "")),
                "k": int(k),
                "nearest_state": str(counts.index[0]) if len(counts) else "",
                "nearest_state_fraction": float(counts.iloc[0] / k) if len(counts) else np.nan,
                "mean_distance": float(distances[local].mean()),
                "min_distance": float(distances[local].min()),
                "neighbor_states": ";".join(f"{state}:{int(counts[state])}" for state in counts.index),
            })
    del full, high_latent, normal_latent
    gc.collect()
    result = pd.DataFrame(rows)
    result.to_csv(out_path, index=False)
    return result, "available"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    out = PROJECT_ROOT / config["outputs"]["step12_dir"]
    out.mkdir(parents=True, exist_ok=True)
    consensus = ad.read_h5ad(PROJECT_ROOT / config["outputs"]["consensus_h5ad"])
    boundary_source = PROJECT_ROOT / config["outputs"]["step9_dir"] / "boundary_state_transitions.csv"
    if boundary_source.exists():
        boundary = pd.read_csv(boundary_source)
    else:
        boundary = pd.DataFrame()
    boundary.to_csv(out / "boundary_validation.csv", index=False)
    if not boundary.empty:
        boundary["nearest_state"] = boundary["nearest_state"].astype(str)
        boundary_summary = boundary.groupby("nearest_state", observed=True).agg(
            boundary_cells=("nearest_state", "size"),
            mean_state_entropy=("state_entropy", "mean"),
            mean_distance=("mean_distance", "mean"),
            mean_neighbor_states=("n_neighbor_states", "mean"),
        ).reset_index().rename(columns={"nearest_state": "state_id"})
    else:
        boundary_summary = pd.DataFrame(columns=["state_id", "boundary_cells", "mean_state_entropy", "mean_distance", "mean_neighbor_states"])
    normal_path = PROJECT_ROOT / config["outputs"]["scvi_h5ad"]
    normal, normal_status = normal_validation(consensus, normal_path, config, out / "normal_validation.csv")
    if not normal.empty:
        normal["nearest_state"] = normal["nearest_state"].astype(str)
        normal_summary = normal.groupby("nearest_state", observed=True).agg(
            normal_cells=("nearest_state", "size"),
            normal_mean_distance=("mean_distance", "mean"),
            normal_min_distance=("min_distance", "mean"),
            normal_nearest_fraction=("nearest_state_fraction", "mean"),
        ).reset_index().rename(columns={"nearest_state": "state_id"})
    else:
        normal_summary = pd.DataFrame(columns=["state_id", "normal_cells", "normal_mean_distance", "normal_min_distance", "normal_nearest_fraction"])
    state_summary = pd.DataFrame({"state_id": consensus.obs["consensus_state"].astype(str).unique()})
    boundary_summary["state_id"] = boundary_summary["state_id"].astype(str)
    normal_summary["state_id"] = normal_summary["state_id"].astype(str)
    state_summary = state_summary.merge(boundary_summary, on="state_id", how="left").merge(normal_summary, on="state_id", how="left")
    state_summary.to_csv(out / "state_auxiliary_summary.csv", index=False)
    write_json(out / "step12_manifest.json", {
        "boundary_status": "available" if not boundary.empty else "not_available",
        "boundary_cells": int(len(boundary)),
        "normal_reference_status": normal_status,
        "normal_cells": int(len(normal)),
        "states_defined_by": "step7 consensus only",
        "boundary_and_normal_participate_in_state_count": False,
        "new_spatial_or_atac_data": False,
        "scvi_training_called": False,
    })
    (PROJECT_ROOT / "reports/stage12_boundary_normal_validation.md").write_text(
        "# Stage 12 boundary and normal validation\n\n"
        f"- Boundary cells evaluated: {len(boundary)}\n"
        f"- Normal reference status: {normal_status}; cells evaluated: {len(normal)}\n"
        "- These cohorts are auxiliary and do not redefine consensus states or state count.\n",
        encoding="utf-8",
    )
    print(f"Step 12 complete: boundary={len(boundary)}, normal={len(normal)} ({normal_status})")


if __name__ == "__main__":
    main()
