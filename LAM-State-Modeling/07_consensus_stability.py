#!/usr/bin/env python3
"""Step 7: equal-weight grid consensus on high-confidence LAM cells."""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from sklearn.metrics import adjusted_rand_score

from state_modeling_utils import (
    PROJECT_ROOT,
    latent_leiden_labels,
    load_config,
    partition_cluster_matches,
    write_json,
)


def grid_id(item: dict) -> str:
    return f"nn{int(item['n_neighbors'])}_res{float(item['resolution']):g}"


def stage6_grid(config: dict) -> list[dict]:
    return [
        {"n_neighbors": int(item["n_neighbors"]), "resolution": float(item["resolution"])}
        for item in config["stage6"]["lam_parameter_grid"]
    ]


def load_base_assignments(high: ad.AnnData, config: dict) -> dict[str, pd.Series]:
    path = PROJECT_ROOT / "results/stage1_6/stage6_lam_only_grid_assignments.csv"
    result: dict[str, pd.Series] = {}
    if path.exists():
        table = pd.read_csv(path, dtype=str)
        key = pd.Series(high.obs.index.astype(str), index=high.obs.index)
        if "scvi_obs_index" in table.columns:
            table = table.set_index("scvi_obs_index")
            for item in stage6_grid(config):
                name = grid_id(item)
                if name in table.columns:
                    aligned = table[name].reindex(key.to_numpy())
                    if aligned.notna().all():
                        result[name] = pd.Series(aligned.to_numpy(), index=high.obs.index, dtype=str)
    for item in stage6_grid(config):
        name = grid_id(item)
        if name not in result:
            result[name] = latent_leiden_labels(
                np.asarray(high.obsm["X_scVI"]),
                high.obs.index,
                item["n_neighbors"],
                item["resolution"],
                int(config["random_seed"]),
                key_added="leiden",
            )
    return result


def matches_with_context(left: pd.Series, right: pd.Series, left_name: str, right_name: str) -> pd.DataFrame:
    out = partition_cluster_matches(left.to_numpy(), right.to_numpy())
    out.insert(0, "left_partition", left_name)
    out.insert(1, "right_partition", right_name)
    return out


def cell_stability(c_final: np.ndarray, labels: np.ndarray) -> pd.DataFrame:
    rows = []
    for state in np.unique(labels):
        indices = np.flatnonzero(labels == state)
        others = np.flatnonzero(labels != state)
        within = c_final[np.ix_(indices, indices)].mean(axis=1) if len(indices) else np.zeros(0)
        between = c_final[np.ix_(indices, others)].max(axis=1) if len(others) else np.zeros(len(indices))
        for local, cell in enumerate(indices):
            margin = float(within[local] - between[local])
            rows.append({
                "cell_position": int(cell),
                "consensus_state": str(state),
                "within_coassignment": float(within[local]),
                "max_between_coassignment": float(between[local]),
                "margin": margin,
                "edge_score": float(1.0 - margin),
            })
    return pd.DataFrame(rows).sort_values("cell_position").reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    out = PROJECT_ROOT / config["outputs"]["step7_dir"]
    out.mkdir(parents=True, exist_ok=True)
    scvi_path = PROJECT_ROOT / config["outputs"]["scvi_h5ad"]
    obj = ad.read_h5ad(scvi_path)
    high_mask = obj.obs["lam_candidate"].astype(bool).to_numpy()
    if int(high_mask.sum()) == 0:
        raise RuntimeError("BLOCKED_INPUT: no pool_high_confidence cells in scVI artifact")
    high = obj[high_mask].copy()
    if "X_scVI" not in high.obsm:
        raise RuntimeError("BLOCKED_INPUT: X_scVI is missing")
    n_cells = len(high)
    base = load_base_assignments(high, config)
    seeds = [int(x) for x in config["step7"]["seed_values"]]
    repeated = {(30, float(res)) for res in (0.2, 0.4, 0.6)}
    all_partitions: dict[str, pd.Series] = {}
    seed_rows = []
    config_rows = []
    matching_frames = []
    c_final = np.zeros((n_cells, n_cells), dtype=np.float32)
    n_configs = len(stage6_grid(config))

    for item in stage6_grid(config):
        name = grid_id(item)
        config_seeds = seeds if (item["n_neighbors"], item["resolution"]) in repeated else [int(config["random_seed"])]
        seed_labels: dict[str, pd.Series] = {f"{name}__seed{int(config['random_seed'])}": base[name]}
        if len(config_seeds) > 1:
            seed_labels = {}
            for seed in config_seeds:
                seed_labels[f"{name}__seed{seed}"] = latent_leiden_labels(
                    np.asarray(high.obsm["X_scVI"]), high.obs.index,
                    item["n_neighbors"], item["resolution"], seed, key_added="leiden",
                )
        for partition_name, labels in seed_labels.items():
            all_partitions[partition_name] = labels
        pairwise_seed = []
        seed_names = list(seed_labels)
        for i, left_name in enumerate(seed_names):
            for right_name in seed_names[i + 1:]:
                ari = float(adjusted_rand_score(seed_labels[left_name], seed_labels[right_name]))
                pairwise_seed.append(ari)
                seed_rows.append({"grid_id": name, "left_seed": left_name, "right_seed": right_name, "ari": ari})
        # Only this configuration's average is added.  Thus five seeds at
        # nn=30 do not increase its final consensus weight.
        config_coassignment = np.zeros((n_cells, n_cells), dtype=np.float32)
        for labels in seed_labels.values():
            encoded = labels.to_numpy(dtype=str)
            config_coassignment += (encoded[:, None] == encoded[None, :]).astype(np.float32)
        config_coassignment /= float(len(seed_labels))
        c_final += config_coassignment / float(n_configs)
        del config_coassignment
        gc.collect()
        config_rows.append({
            "grid_id": name,
            "n_neighbors": int(item["n_neighbors"]),
            "resolution": float(item["resolution"]),
            "n_seeds": int(len(seed_labels)),
            "seed_ari_mean": float(np.mean(pairwise_seed)) if pairwise_seed else np.nan,
            "seed_ari_median": float(np.median(pairwise_seed)) if pairwise_seed else np.nan,
        })
        # Grid-level matches are retained, including split/merge overlaps.
        for partition_name, labels in seed_labels.items():
            reference = base[name]
            matching_frames.append(matches_with_context(reference, labels, f"{name}__base", partition_name))

    # Also retain direct matching among the nine configuration-level
    # partitions.  These rows are descriptive split/merge correspondence;
    # they do not change the equal-weight co-assignment calculation above.
    base_names = list(base)
    for i, left_name in enumerate(base_names):
        for right_name in base_names[i + 1:]:
            matching_frames.append(matches_with_context(base[left_name], base[right_name], left_name, right_name))

    np.savez_compressed(out / "coassignment_matrix.npz", coassignment=c_final, cell_ids=high.obs.index.astype(str).to_numpy())
    pd.DataFrame(config_rows).to_csv(out / "grid_configuration_summary.csv", index=False)
    pd.DataFrame(seed_rows).to_csv(out / "seed_stability.csv", index=False)
    pd.concat(matching_frames, ignore_index=True).to_csv(out / "cluster_matching_across_grid.csv", index=False)

    distance = squareform(1.0 - c_final, checks=False)
    linkage_matrix = linkage(distance, method="average")
    np.savez_compressed(out / "consensus_dendrogram.npz", linkage=linkage_matrix)
    del distance

    rng = np.random.default_rng(int(config["random_seed"]))
    pair_count = min(int(config["step7"].get("stability_pair_sample_size", 1_000_000)), n_cells * (n_cells - 1) // 2)
    pair_i = rng.integers(0, n_cells, size=pair_count, dtype=np.int32)
    pair_j = rng.integers(0, n_cells, size=pair_count, dtype=np.int32)
    keep = pair_i != pair_j
    pair_i, pair_j = pair_i[keep], pair_j[keep]
    del keep
    candidate_rows = []
    max_k = min(int(config["step7"].get("consensus_max_clusters", 100)), n_cells - 1)
    for k in range(int(config["step7"].get("consensus_min_clusters", 2)), max_k + 1):
        labels = fcluster(linkage_matrix, t=k, criterion="maxclust").astype(str)
        same = labels[pair_i] == labels[pair_j]
        if same.any():
            within = float(c_final[pair_i[same], pair_j[same]].mean())
        else:
            within = 0.0
        if (~same).any():
            between = float(c_final[pair_i[~same], pair_j[~same]].mean())
        else:
            between = 0.0
        # The cut-selection diagnostic follows the same nine-configuration
        # weighting as C_final.  Raw seed partitions are retained for audit
        # and seed stability, but never receive extra consensus weight.
        grid_ari = [float(adjusted_rand_score(labels, part.to_numpy())) for part in base.values()]
        candidate_rows.append({
            "n_clusters_requested": int(k),
            "n_clusters_observed": int(pd.Series(labels).nunique()),
            "within_coassignment": within,
            "between_coassignment": between,
            "separation": within - between,
            "mean_ari_to_configurations": float(np.mean(grid_ari)),
        })
    candidates = pd.DataFrame(candidate_rows)
    candidates["separation_rank"] = candidates["separation"].rank(method="average", ascending=False)
    candidates["configuration_ari_rank"] = candidates["mean_ari_to_configurations"].rank(method="average", ascending=False)
    candidates["consensus_score"] = (
        candidates["separation"].rank(pct=True) + candidates["mean_ari_to_configurations"].rank(pct=True)
    ) / 2.0
    best_score = float(candidates["consensus_score"].max())
    tolerance = float(config["step7"].get("consensus_score_tolerance", 0.05))
    eligible = candidates[candidates["consensus_score"] >= best_score - tolerance]
    selected_k = int(eligible.sort_values("n_clusters_observed").iloc[0]["n_clusters_requested"])
    consensus_labels = fcluster(linkage_matrix, t=selected_k, criterion="maxclust").astype(str)
    candidates["selected"] = candidates["n_clusters_requested"].eq(selected_k)
    candidates.to_csv(out / "consensus_cut_candidates.csv", index=False)

    stability = cell_stability(c_final, consensus_labels)
    assignments = pd.DataFrame({
        "scvi_obs_index": high.obs.index.astype(str),
        "analysis_cell_id": high.obs.get("analysis_cell_id", pd.Series(high.obs.index, index=high.obs.index)).astype(str).to_numpy(),
        "patient_id": high.obs["patient_id"].astype(str).to_numpy(),
        "dataset": high.obs["dataset"].astype(str).to_numpy(),
        "assay": high.obs["assay"].astype(str).to_numpy(),
        "consensus_state": consensus_labels,
    })
    assignments = assignments.join(stability.drop(columns=["consensus_state"]))
    assignments.to_csv(out / "state_consensus_assignments.csv", index=False)
    state_summary = assignments.groupby("consensus_state", observed=True).agg(
        cells=("consensus_state", "size"),
        patients=("patient_id", "nunique"),
        datasets=("dataset", "nunique"),
        assays=("assay", "nunique"),
        mean_within_coassignment=("within_coassignment", "mean"),
        median_within_coassignment=("within_coassignment", "median"),
        mean_margin=("margin", "mean"),
        median_margin=("margin", "median"),
        mean_edge_score=("edge_score", "mean"),
    ).reset_index()
    state_summary.to_csv(out / "state_stability_summary.csv", index=False)

    consensus = high.copy()
    consensus.obs["consensus_state"] = pd.Categorical(consensus_labels)
    for column in ["within_coassignment", "max_between_coassignment", "margin", "edge_score"]:
        consensus.obs[column] = stability[column].to_numpy()
    consensus.uns["step7_consensus"] = {
        "n_cells": n_cells,
        "n_grid_configurations": n_configs,
        "configuration_equal_weight": True,
        "seed_values": seeds,
        "selected_n_clusters": selected_k,
        "selected_cut_rule": "max combined separation/partition-ARI rank within tolerance; no fixed K",
        "scvi_training_called": False,
    }
    consensus.write(PROJECT_ROOT / config["outputs"]["consensus_h5ad"])
    write_json(out / "step7_manifest.json", consensus.uns["step7_consensus"] | {
        "grid_configurations": config_rows,
        "n_partitions": len(all_partitions),
        "repeated_configuration_seed_averaging": True,
        "coassignment_dtype": "float32",
        "full_distance_used_for_average_linkage": True,
    })
    print(f"Step 7 complete: {n_cells} high-confidence cells, {selected_k} consensus states")


if __name__ == "__main__":
    main()
