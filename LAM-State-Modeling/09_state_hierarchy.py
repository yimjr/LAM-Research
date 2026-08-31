#!/usr/bin/env python3
"""Step 9: characterize hierarchy, connectivity, and boundary transitions."""

from __future__ import annotations

import argparse
import gc
import os
from collections import Counter

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from scipy.spatial.distance import cdist
from sklearn.neighbors import NearestNeighbors

from state_modeling_utils import PROJECT_ROOT, latent_leiden_labels, load_config, partition_cluster_matches, write_json


def grid_id(item: dict) -> str:
    return f"nn{int(item['n_neighbors'])}_res{float(item['resolution']):g}"


def entropy(probabilities: np.ndarray) -> float:
    values = probabilities[probabilities > 0]
    return float(-(values * np.log2(values)).sum()) if len(values) else 0.0


def boundary_transitions(consensus: ad.AnnData, full_path, config: dict, out_path) -> pd.DataFrame:
    full = ad.read_h5ad(full_path, backed="r")
    boundary_mask = full.obs.get("boundary", pd.Series(False, index=full.obs.index)).astype(bool).to_numpy()
    if not boundary_mask.any():
        return pd.DataFrame(columns=["boundary_analysis_cell_id", "nearest_state", "second_nearest_state", "n_neighbor_states", "state_entropy", "mean_distance"])
    high_index = consensus.obs.index.astype(str)
    full_index = full.obs.index.astype(str)
    positions = pd.Series(np.arange(len(full_index)), index=full_index)
    high_positions = positions.reindex(high_index).dropna().astype(int).to_numpy()
    boundary_positions = np.flatnonzero(boundary_mask)
    high_latent = np.asarray(full.obsm["X_scVI"][high_positions], dtype=np.float32)
    boundary_latent = np.asarray(full.obsm["X_scVI"][boundary_positions], dtype=np.float32)
    high_states = consensus.obs["consensus_state"].astype(str).to_numpy()
    k = min(int(config["step12"].get("boundary_k", 15)), len(high_states))
    nn = NearestNeighbors(n_neighbors=max(1, k), metric="euclidean").fit(high_latent)
    rows = []
    batch_size = int(config["step12"].get("query_batch_size", 512))
    for start in range(0, len(boundary_positions), batch_size):
        stop = min(start + batch_size, len(boundary_positions))
        distances, neighbors = nn.kneighbors(boundary_latent[start:stop])
        for local in range(stop - start):
            counts = Counter(high_states[neighbors[local]])
            total = float(sum(counts.values()))
            ordered = sorted(counts, key=lambda state: (-counts[state], state))
            probs = np.asarray([counts[state] / total for state in ordered], dtype=float)
            rows.append({
                "boundary_analysis_cell_id": str(full.obs.iloc[boundary_positions[start + local]].get("analysis_cell_id", full.obs.index[boundary_positions[start + local]])),
                "boundary_patient_id": str(full.obs.iloc[boundary_positions[start + local]].get("patient_id", "")),
                "boundary_dataset": str(full.obs.iloc[boundary_positions[start + local]].get("dataset", "")),
                "k": int(k),
                "nearest_state": str(ordered[0]) if ordered else "",
                "second_nearest_state": str(ordered[1]) if len(ordered) > 1 else "",
                "n_neighbor_states": int(len(ordered)),
                "state_entropy": entropy(probs),
                "mean_distance": float(distances[local].mean()),
                "neighbor_state_counts": ";".join(f"{state}:{counts[state]}" for state in ordered),
            })
    del full, high_latent, boundary_latent
    gc.collect()
    result = pd.DataFrame(rows)
    result.to_csv(out_path, index=False)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    out = PROJECT_ROOT / config["outputs"]["step9_dir"]
    out.mkdir(parents=True, exist_ok=True)
    consensus = ad.read_h5ad(PROJECT_ROOT / config["outputs"]["consensus_h5ad"])
    states = consensus.obs["consensus_state"].astype(str)
    state_ids = sorted(states.unique(), key=lambda value: (len(value), value))
    latent = np.asarray(consensus.obsm["X_scVI"], dtype=np.float32)

    centroids = np.vstack([latent[states.to_numpy() == state].mean(axis=0) for state in state_ids])
    distance = cdist(centroids, centroids, metric="euclidean")
    distance_df = pd.DataFrame(distance, index=state_ids, columns=state_ids)
    distance_df.index.name = "state_id"
    distance_df.to_csv(out / "state_distance_matrix.csv")

    graph = ad.AnnData(X=np.zeros((len(consensus), 1), dtype=np.float32), obs=consensus.obs[["consensus_state"]].copy())
    graph.obsm["X_scVI"] = latent
    sc.pp.neighbors(graph, n_neighbors=min(int(config["step12"].get("boundary_k", 15)), len(graph) - 1), use_rep="X_scVI", key_added="state_neighbors")
    sc.tl.paga(graph, groups="consensus_state", neighbors_key="state_neighbors")
    connectivities = graph.uns["paga"]["connectivities"]
    if sparse.issparse(connectivities):
        connectivities = connectivities.toarray()
    connectivity_rows = []
    for i, left in enumerate(state_ids):
        for j, right in enumerate(state_ids):
            if i >= j:
                continue
            connectivity_rows.append({
                "state_a": left,
                "state_b": right,
                "paga_connectivity": float(connectivities[i, j]),
                "latent_centroid_distance": float(distance[i, j]),
            })
    pd.DataFrame(connectivity_rows).to_csv(out / "state_connectivity.csv", index=False)

    assignments_path = PROJECT_ROOT / "results/stage1_6/stage6_lam_only_grid_assignments.csv"
    grid_rows = []
    if assignments_path.exists():
        grid = pd.read_csv(assignments_path, dtype=str).set_index("scvi_obs_index")
        aligned = grid.reindex(consensus.obs.index.astype(str))
        for item in config["stage6"]["lam_parameter_grid"]:
            name = grid_id(item)
            if name not in aligned:
                continue
            matches = partition_cluster_matches(aligned[name].to_numpy(), states.to_numpy())
            matches.insert(0, "grid_id", name)
            matches.insert(1, "n_neighbors", int(item["n_neighbors"]))
            matches.insert(2, "resolution", float(item["resolution"]))
            grid_rows.append(matches)
    split_merge = pd.concat(grid_rows, ignore_index=True) if grid_rows else pd.DataFrame()
    split_merge.to_csv(out / "state_split_merge_tree.csv", index=False)

    # A coarse parent is deliberately descriptive: it comes from the coarsest
    # existing grid partition and never changes the consensus state labels.
    parent_name = grid_id({"n_neighbors": 50, "resolution": 0.2})
    if parent_name in aligned:
        parent_by_state = {}
        for state in state_ids:
            values = aligned.loc[states.to_numpy() == state, parent_name].astype(str)
            parent_by_state[state] = values.mode().iloc[0] if not values.empty else ""
    else:
        parent_by_state = {state: "" for state in state_ids}
    consensus.obs["parent_state"] = states.map(parent_by_state).astype(str).to_numpy()
    parent_counts = states.map(parent_by_state).value_counts()
    consensus.obs["substate_role"] = ["substate" if parent_counts.get(parent_by_state[state], 0) > 1 else "parent_or_singleton" for state in states]

    boundary_path = out / "boundary_state_transitions.csv"
    boundary = boundary_transitions(consensus, PROJECT_ROOT / config["outputs"]["scvi_h5ad"], config, boundary_path)
    summary = {
        "n_consensus_states": int(len(state_ids)),
        "n_cells": int(len(consensus)),
        "paga_edges_reported": int(len(connectivity_rows)),
        "boundary_cells": int(len(boundary)),
        "parent_partition": parent_name,
        "scvi_training_called": False,
    }
    consensus.uns["step9_hierarchy"] = summary
    consensus.write(PROJECT_ROOT / config["outputs"]["hierarchy_h5ad"])
    write_json(out / "step9_manifest.json", summary)
    report = PROJECT_ROOT / "reports/stage9_state_hierarchy.md"
    report.write_text(
        "# Stage 9 state hierarchy\n\n"
        f"- Consensus states: {len(state_ids)}\n"
        f"- Cells: {len(consensus)}\n"
        f"- Boundary cells evaluated: {len(boundary)}\n"
        "- State labels remain those from Step 7; parent/substate labels are descriptive only.\n"
        "- Boundary and normal cohorts do not define the number of LAM states.\n",
        encoding="utf-8",
    )
    print(f"Step 9 complete: {len(state_ids)} states, {len(boundary)} boundary transitions")


if __name__ == "__main__":
    main()
