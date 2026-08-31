#!/usr/bin/env python3
"""Evaluate Stage 6 without making upstream state correspondence a gate."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors

from state_modeling_utils import PROJECT_ROOT, load_config, write_json


def ari_if_variable(labels: pd.Series, clusters: pd.Series) -> float | None:
    values = labels.astype(str)
    if values.nunique() < 2 or clusters.astype(str).nunique() < 2:
        return None
    return float(adjusted_rand_score(values, clusters.astype(str)))


def posthoc_summary(obj: ad.AnnData, cluster_col: str) -> pd.DataFrame:
    # Historical upstream tables may contain FIGF.  Preserve those raw
    # annotations in AnnData, but do not expose the legacy alias in new
    # marker/program reports; the canonical symbol is VEGFD.
    upstream = [
        col for col in obj.obs.columns
        if str(col).startswith("upstream_") and "FIGF" not in str(col).upper()
    ]
    if not upstream or cluster_col not in obj.obs:
        return pd.DataFrame()
    rows = []
    for cluster, frame in obj.obs.groupby(cluster_col, observed=True):
        row = {cluster_col: str(cluster), "cells": int(len(frame))}
        for col in upstream:
            values = frame[col]
            if pd.api.types.is_bool_dtype(values) or pd.api.types.is_numeric_dtype(values):
                numeric = pd.to_numeric(values, errors="coerce")
                if numeric.notna().any():
                    row[f"mean__{col}"] = float(numeric.mean())
        rows.append(row)
    return pd.DataFrame(rows)


def normal_neighbors(obj: ad.AnnData, high_obj: ad.AnnData, config: dict) -> tuple[str, dict]:
    output = PROJECT_ROOT / "results/stage1_6/normal_to_lam_neighbors.csv"
    normal_mask = obj.obs["analysis_role"].astype(str).eq("normal_reference").to_numpy()
    high_mask = obj.obs["lam_candidate"].astype(bool).to_numpy()
    if not normal_mask.any():
        pd.DataFrame().to_csv(output, index=False)
        return "not_available", {"status": "not_available", "reason": "normal reference is absent"}
    if len(high_obj) == 0 or "X_scVI" not in high_obj.obsm:
        pd.DataFrame().to_csv(output, index=False)
        return "not_evaluable", {"status": "not_evaluable", "reason": "missing high-confidence cells or scVI latent"}
    k = min(int(config["stage6"]["normal_k"]), int(high_mask.sum()))
    high_latent = np.asarray(high_obj.obsm["X_scVI"])
    normal_latent = np.asarray(obj.obsm["X_scVI"])[normal_mask]
    nn = NearestNeighbors(n_neighbors=max(1, k), metric="euclidean").fit(high_latent)
    high_obs = high_obj.obs.reset_index(drop=True)
    normal_obs = obj.obs.loc[normal_mask].reset_index(drop=True)
    rows = []
    # Query normal cells in bounded batches.  A single brute-force query for
    # ~33k normal cells against all high-confidence cells can otherwise
    # materialize a multi-GB distance matrix and obscure the core checkpoint.
    batch_size = int(config["stage6"].get("normal_query_batch_size", 512))
    for start in range(0, len(normal_obs), batch_size):
        stop = min(start + batch_size, len(normal_obs))
        distances, indices = nn.kneighbors(normal_latent[start:stop])
        for local_i in range(stop - start):
            i = start + local_i
            nearest = high_obs.iloc[indices[local_i]]
            rows.append({
                "normal_analysis_cell_id": str(normal_obs.iloc[i].get("analysis_cell_id", normal_obs.index[i])),
                "normal_patient_id": str(normal_obs.iloc[i].get("patient_id", "")),
                "normal_dataset": str(normal_obs.iloc[i].get("dataset", "")),
                "k": int(k),
                "mean_distance": float(distances[local_i].mean()),
                "nearest_lam_patient_ids": ";".join(sorted(set(nearest["patient_id"].astype(str)))),
                "nearest_lam_clusters": ";".join(sorted(set(nearest["leiden_scvi_lam_only"].astype(str)))),
            })
    pd.DataFrame(rows).to_csv(output, index=False)
    return "available", {"status": "available", "k": int(k), "n_normal_cells": int(len(rows)), "output": str(output)}


def boundary_connectivity(obj: ad.AnnData, high_obj: ad.AnnData, config: dict) -> dict:
    """Describe how boundary cells connect to LAM-only latent clusters.

    This is deliberately auxiliary: no boundary statistic is used in the
    Stage 6 Go/No-Go decision.
    """
    output = PROJECT_ROOT / "results/stage1_6/stage6_boundary_state_connectivity.csv"
    boundary_mask = obj.obs["boundary"].astype(bool).to_numpy() if "boundary" in obj.obs else np.zeros(obj.n_obs, dtype=bool)
    columns = [
        "boundary_analysis_cell_id", "boundary_dataset", "boundary_patient_id",
        "k", "mean_distance", "n_nearest_lam_clusters", "nearest_lam_clusters",
    ]
    if not boundary_mask.any() or len(high_obj) == 0 or "X_scVI" not in high_obj.obsm or "leiden_scvi_lam_only" not in high_obj.obs:
        pd.DataFrame(columns=columns).to_csv(output, index=False)
        return {"status": "not_evaluable", "reason": "missing boundary cells, high-confidence latent, or LAM-only clusters", "output": str(output)}

    k = min(int(config["stage6"]["normal_k"]), len(high_obj))
    nn = NearestNeighbors(n_neighbors=max(1, k), metric="euclidean").fit(np.asarray(high_obj.obsm["X_scVI"]))
    boundary_obs = obj.obs.loc[boundary_mask].reset_index(drop=True)
    boundary_latent = np.asarray(obj.obsm["X_scVI"])[boundary_mask]
    high_obs = high_obj.obs.reset_index(drop=True)
    rows = []
    batch_size = int(config["stage6"].get("normal_query_batch_size", 512))
    for start in range(0, len(boundary_obs), batch_size):
        stop = min(start + batch_size, len(boundary_obs))
        distances, indices = nn.kneighbors(boundary_latent[start:stop])
        for local_i in range(stop - start):
            i = start + local_i
            nearest = high_obs.iloc[indices[local_i]]
            clusters = sorted(set(nearest["leiden_scvi_lam_only"].astype(str)))
            rows.append({
                "boundary_analysis_cell_id": str(boundary_obs.iloc[i].get("analysis_cell_id", boundary_obs.index[i])),
                "boundary_dataset": str(boundary_obs.iloc[i].get("dataset", "")),
                "boundary_patient_id": str(boundary_obs.iloc[i].get("patient_id", "")),
                "k": int(k),
                "mean_distance": float(distances[local_i].mean()),
                "n_nearest_lam_clusters": int(len(clusters)),
                "nearest_lam_clusters": ";".join(clusters),
            })
    pd.DataFrame(rows, columns=columns).to_csv(output, index=False)
    return {"status": "available", "k": int(k), "n_boundary_cells": int(len(rows)), "output": str(output)}


def parameter_grid(config: dict) -> list[dict]:
    stage6 = config["stage6"]
    default = {
        "n_neighbors": int(stage6["lam_n_neighbors"]),
        "resolution": float(stage6["lam_leiden_resolution"]),
    }
    raw_grid = stage6.get("lam_parameter_grid", [default])
    grid: list[dict] = []
    seen: set[tuple[int, float]] = set()
    for item in raw_grid:
        current = {
            "n_neighbors": int(item["n_neighbors"]),
            "resolution": float(item["resolution"]),
        }
        if current["n_neighbors"] < 2 or current["resolution"] <= 0:
            raise ValueError(f"invalid Stage 6 clustering parameters: {current}")
        key = (current["n_neighbors"], current["resolution"])
        if key not in seen:
            grid.append(current)
            seen.add(key)
    default_key = (default["n_neighbors"], default["resolution"])
    if default_key not in seen:
        grid.append(default)
    return grid


def grid_id(n_neighbors: int, resolution: float) -> str:
    return f"nn{n_neighbors}_res{resolution:g}"


def lam_only_clusters(high_obj: ad.AnnData, n_neighbors: int, resolution: float, seed: int) -> pd.Series:
    """Cluster only high-confidence cells using the existing scVI latent.

    The temporary AnnData intentionally contains no expression matrix.  This
    prevents accidental use of counts/log-normalized X and keeps each grid
    run bounded in memory.
    """
    if len(high_obj) == 0 or "X_scVI" not in high_obj.obsm:
        return pd.Series("", index=high_obj.obs.index, dtype=str)
    graph = ad.AnnData(
        X=np.zeros((len(high_obj), 1), dtype=np.float32),
        obs=pd.DataFrame(index=high_obj.obs.index.copy()),
    )
    graph.obsm["X_scVI"] = np.asarray(high_obj.obsm["X_scVI"], dtype=np.float32)
    effective_neighbors = min(int(n_neighbors), len(high_obj) - 1)
    if effective_neighbors < 2:
        return pd.Series("0", index=high_obj.obs.index, dtype=str)
    sc.pp.neighbors(
        graph,
        n_neighbors=effective_neighbors,
        use_rep="X_scVI",
        random_state=seed,
        key_added="neighbors_scvi_lam_only",
    )
    sc.tl.leiden(
        graph,
        resolution=float(resolution),
        key_added="leiden_scvi_lam_only",
        neighbors_key="neighbors_scvi_lam_only",
        random_state=seed,
        flavor="igraph",
        directed=False,
    )
    return pd.Series(graph.obs["leiden_scvi_lam_only"].astype(str).to_numpy(), index=high_obj.obs.index)


def summarize_clusters(obs: pd.DataFrame, labels: pd.Series, config: dict) -> tuple[dict, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cluster_col = "leiden_scvi_lam_only"
    frame = obs.copy()
    frame[cluster_col] = labels.astype(str).to_numpy()
    coverage = (
        frame.groupby(cluster_col, observed=True)
        .agg(cells=(cluster_col, "size"), patients=("patient_id", "nunique"), datasets=("dataset", "nunique"), assays=("assay", "nunique"))
        .reset_index()
    )
    patient_counts = (
        frame.groupby([cluster_col, "patient_id"], observed=True)
        .size()
        .reset_index(name="cells")
    )
    qualified_patients = (
        patient_counts[patient_counts["cells"] >= 5]
        .groupby(cluster_col, observed=True)["patient_id"]
        .nunique()
    )
    driver_rows = []
    for field in ["patient_id", "dataset", "assay"]:
        variable = bool(field in frame and frame[field].astype(str).nunique() >= 2)
        value = ari_if_variable(frame[field], frame[cluster_col]) if variable else None
        driver_rows.append({"driver": field, "ari": value, "variable": variable})
    driver_table = pd.DataFrame(driver_rows)
    summary = {
        "n_latent_clusters": int(len(coverage)),
        "singleton_cluster_count": int((coverage["cells"] == 1).sum()) if len(coverage) else 0,
        "cluster_count_with_fewer_than_5_cells": int((coverage["cells"] < 5).sum()) if len(coverage) else 0,
        "largest_cluster_cells": int(coverage["cells"].max()) if len(coverage) else 0,
        "median_cluster_size": float(coverage["cells"].median()) if len(coverage) else 0.0,
        "shared_cluster_count": int((coverage["patients"] >= int(config["stage6"]["min_shared_patients"])).sum()) if len(coverage) else 0,
        "shared_qualified_cluster_count": int((qualified_patients >= int(config["stage6"]["min_shared_patients"])).sum()) if len(qualified_patients) else 0,
        "n_patients": int(frame["patient_id"].astype(str).nunique()) if len(frame) else 0,
    }
    return summary, coverage, patient_counts, driver_table


def obs_strings(obs: pd.DataFrame, column: str) -> np.ndarray:
    if column not in obs:
        return np.full(len(obs), "", dtype=object)
    return obs[column].astype(str).to_numpy()


def write_grid_assignments(high_obj: ad.AnnData, assignments: dict[str, pd.Series]) -> None:
    table = pd.DataFrame({
        "scvi_obs_index": high_obj.obs.index.astype(str),
        "analysis_cell_id": obs_strings(high_obj.obs, "analysis_cell_id"),
        "patient_id": obs_strings(high_obj.obs, "patient_id"),
        "dataset": obs_strings(high_obj.obs, "dataset"),
    })
    for name, labels in assignments.items():
        table[name] = labels.astype(str).to_numpy()
    table.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_lam_only_grid_assignments.csv", index=False)


def pairwise_grid_ari(assignments: dict[str, pd.Series]) -> pd.DataFrame:
    names = list(assignments)
    rows = []
    for i, left in enumerate(names):
        for right in names[i + 1:]:
            rows.append({
                "grid_id_a": left,
                "grid_id_b": right,
                "partition_ari": float(adjusted_rand_score(assignments[left].astype(str), assignments[right].astype(str))),
            })
    return pd.DataFrame(rows, columns=["grid_id_a", "grid_id_b", "partition_ari"])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    input_path = PROJECT_ROOT / config["outputs"]["scvi_h5ad"]
    if not input_path.exists():
        print(f"BLOCKED_INPUT: missing {input_path}; run script 05 first", file=sys.stderr)
        return 2
    obj = ad.read_h5ad(input_path)
    high = obj.obs["lam_candidate"].astype(bool).to_numpy()
    high_obj = obj[high].copy()
    blockers: list[str] = []
    if len(high_obj) == 0:
        blockers.append("no pool_high_confidence cells in scVI output")
    for field in ["patient_id", "donor_id"]:
        if field not in high_obj.obs or high_obj.obs[field].astype(str).isin({"", "nan", "None", "unknown"}).any():
            blockers.append(f"unresolved {field} in pool_high_confidence")
    if "X_scVI" not in high_obj.obsm:
        blockers.append("scVI latent representation is missing")

    grid = parameter_grid(config)
    assignments: dict[str, pd.Series] = {}
    grid_rows: list[dict] = []
    coverage_frames: list[pd.DataFrame] = []
    patient_frames: list[pd.DataFrame] = []
    driver_frames: list[pd.DataFrame] = []
    for params in grid:
        name = grid_id(params["n_neighbors"], params["resolution"])
        labels = lam_only_clusters(high_obj, params["n_neighbors"], params["resolution"], int(config["random_seed"]))
        assignments[name] = labels
        summary, coverage, patient_counts, driver_table = summarize_clusters(high_obj.obs, labels, config)
        grid_rows.append({
            "grid_id": name,
            "n_neighbors": params["n_neighbors"],
            "resolution": params["resolution"],
            **summary,
            **{f"{row['driver']}_ari": row["ari"] for row in driver_table.to_dict("records")},
        })
        coverage_frames.append(coverage.assign(grid_id=name, n_neighbors=params["n_neighbors"], resolution=params["resolution"]))
        patient_frames.append(patient_counts.assign(grid_id=name, n_neighbors=params["n_neighbors"], resolution=params["resolution"]))
        driver_frames.append(driver_table.assign(grid_id=name, n_neighbors=params["n_neighbors"], resolution=params["resolution"]))

    grid_summary = pd.DataFrame(grid_rows)
    grid_coverage = pd.concat(coverage_frames, ignore_index=True) if coverage_frames else pd.DataFrame()
    grid_patient_counts = pd.concat(patient_frames, ignore_index=True) if patient_frames else pd.DataFrame()
    grid_drivers = pd.concat(driver_frames, ignore_index=True) if driver_frames else pd.DataFrame()
    pairwise = pairwise_grid_ari(assignments)
    grid_summary.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_parameter_grid.csv", index=False)
    grid_coverage.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_grid_cluster_patient_coverage.csv", index=False)
    grid_patient_counts.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_grid_cluster_patient_counts.csv", index=False)
    grid_drivers.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_grid_driver_ari.csv", index=False)
    pairwise.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_grid_pairwise_ari.csv", index=False)
    write_grid_assignments(high_obj, assignments)

    reference_n = int(config["stage6"]["lam_n_neighbors"])
    reference_resolution = float(config["stage6"]["lam_leiden_resolution"])
    reference_name = grid_id(reference_n, reference_resolution)
    if reference_name not in assignments:
        raise RuntimeError("reference Stage 6 parameters were not present in the parameter grid")
    reference_labels = assignments[reference_name]
    lam_cluster_col = "leiden_scvi_lam_only"
    high_obj.obs[lam_cluster_col] = reference_labels.astype(str).to_numpy()
    reference_summary, coverage, cluster_patient, driver_table = summarize_clusters(high_obj.obs, reference_labels, config)
    coverage.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_cluster_patient_coverage.csv", index=False)
    cluster_patient.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_cluster_patient_counts.csv", index=False)
    driver_table.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_driver_ari.csv", index=False)

    full_cohort_cluster_count = int(obj.obs["leiden_scvi"].nunique()) if "leiden_scvi" in obj.obs else None
    max_ari = float(config["stage6"]["max_driver_ari"])
    variable_driver_rows = driver_table[driver_table["variable"] & driver_table["ari"].notna()]
    dominant_drivers = variable_driver_rows[variable_driver_rows["ari"] > max_ari]["driver"].astype(str).tolist()
    shared_cluster = reference_summary["shared_cluster_count"] > 0
    only_one_patient = reference_summary["n_patients"] <= 1

    if blockers:
        status = "BLOCKED_INPUT"
        reasons = blockers
    elif reference_summary["n_latent_clusters"] < int(config["stage6"]["min_latent_clusters"]):
        status = "NO_GO"
        reasons = ["no LAM internal latent structure"]
    elif only_one_patient or not shared_cluster:
        status = "NO_GO"
        reasons = ["latent structure is confined to a single independent patient"]
    elif dominant_drivers:
        status = "NO_GO"
        reasons = [f"latent structure is mainly driven by: {', '.join(dominant_drivers)}"]
    else:
        status = "GO"
        reasons = ["high-confidence LAM latent structure is multi-cluster and cross-patient without a dominant measured driver"]

    posthoc = posthoc_summary(high_obj, lam_cluster_col)
    posthoc.to_csv(PROJECT_ROOT / "results/stage1_6/stage6_posthoc_upstream_summary.csv", index=False)
    boundary_audit = boundary_connectivity(obj, high_obj, config)
    _, normal_audit = normal_neighbors(obj, high_obj, config)
    scvi_training = obj.uns.get("state_model_scvi", {})
    training_mode = str(scvi_training.get("training_mode", "unknown"))
    pairwise_mean = float(pairwise["partition_ari"].mean()) if len(pairwise) else None
    pairwise_min = float(pairwise["partition_ari"].min()) if len(pairwise) else None
    pairwise_max = float(pairwise["partition_ari"].max()) if len(pairwise) else None
    payload = {
        "status": status,
        "reasons": reasons,
        "n_high_confidence_cells": int(len(high_obj)),
        "n_high_confidence_patients": reference_summary["n_patients"],
        "n_latent_clusters": reference_summary["n_latent_clusters"],
        "cluster_label": lam_cluster_col,
        "reference_parameters": {"n_neighbors": reference_n, "resolution": reference_resolution, "grid_id": reference_name},
        "full_cohort_cluster_count_retained_for_audit": full_cohort_cluster_count,
        **{key: reference_summary[key] for key in ["singleton_cluster_count", "cluster_count_with_fewer_than_5_cells", "largest_cluster_cells", "median_cluster_size", "shared_cluster_count", "shared_qualified_cluster_count"]},
        "shared_cluster_across_patients": shared_cluster,
        "dominant_drivers": dominant_drivers,
        "driver_ari_threshold": max_ari,
        "parameter_grid_output": "results/stage1_6/stage6_parameter_grid.csv",
        "parameter_grid_size": len(grid),
        "parameter_grid_pairwise_partition_ari": {"mean": pairwise_mean, "min": pairwise_min, "max": pairwise_max, "output": "results/stage1_6/stage6_grid_pairwise_ari.csv"},
        "boundary_auxiliary": boundary_audit,
        "normal_reference": normal_audit,
        "scvi_training_mode": training_mode,
        "scvi_max_epochs": scvi_training.get("max_epochs"),
        "posthoc_interpretation": "upstream state/program correspondence is descriptive only; unmatched clusters remain novel_or_unexplained candidates",
        "tuning_interpretation": "the grid is reported for stability assessment; no configuration is selected by cluster count alone",
    }
    write_json(PROJECT_ROOT / config["outputs"]["stage6_json"], payload)
    write_json(PROJECT_ROOT / "results/stage1_6/stage6_parameter_grid_summary.json", {"reference_parameters": payload["reference_parameters"], "grid": grid_summary.to_dict("records"), "pairwise_partition_ari": pairwise.to_dict("records"), "selection_rule": payload["tuning_interpretation"]})
    report = [
        "# Stage 6 checkpoint",
        "",
        f"- Status: **{status}**",
        f"- High-confidence candidate cells: {payload['n_high_confidence_cells']}",
        f"- Independent patients: {payload['n_high_confidence_patients']}",
        f"- Reference LAM-only parameters: n_neighbors={reference_n}, resolution={reference_resolution}",
        f"- Reference LAM-only scVI latent clusters: {payload['n_latent_clusters']}",
        f"- Full-cohort clusters retained for audit only: {payload['full_cohort_cluster_count_retained_for_audit']}",
        f"- Reference singleton clusters: {payload['singleton_cluster_count']}",
        f"- Reference clusters with fewer than 5 cells: {payload['cluster_count_with_fewer_than_5_cells']}",
        f"- Reference median/largest cluster size: {payload['median_cluster_size']}/{payload['largest_cluster_cells']}",
        f"- Reference clusters shared by at least {config['stage6']['min_shared_patients']} patients: {payload['shared_cluster_count']}",
        f"- Reference clusters shared by at least {config['stage6']['min_shared_patients']} patients with ≥5 cells per patient: {payload['shared_qualified_cluster_count']}",
        f"- Shared cluster across patients: {payload['shared_cluster_across_patients']}",
        f"- Dominant drivers: {', '.join(dominant_drivers) if dominant_drivers else 'none detected'}",
        f"- scVI training mode: {training_mode}",
        "",
        "## Parameter grid",
        "The nine configurations are reported as a stability analysis. Cluster count alone is not used to choose a configuration.",
        "",
        "| grid_id | n_neighbors | resolution | clusters | singletons | <5 cells | largest | median | shared ≥2 patients | shared ≥2 patients and ≥5 cells/patient | patient ARI | dataset ARI | assay ARI |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in grid_summary.to_dict("records"):
        report.append(
            f"| {row['grid_id']} | {row['n_neighbors']} | {row['resolution']:g} | {row['n_latent_clusters']} | {row['singleton_cluster_count']} | {row['cluster_count_with_fewer_than_5_cells']} | {row['largest_cluster_cells']} | {row['median_cluster_size']:g} | {row['shared_cluster_count']} | {row['shared_qualified_cluster_count']} | {row.get('patient_id_ari')} | {row.get('dataset_ari')} | {row.get('assay_ari')} |"
        )
    report.extend([
        "",
        f"Pairwise partition ARI across grid: mean={pairwise_mean}, min={pairwise_min}, max={pairwise_max}.",
        "",
        "## Decision reasons",
        *[f"- {reason}" for reason in reasons],
        "",
        "## Interpretation boundary",
        "The full-cohort Leiden labels are retained for audit only. The Go/No-Go decision uses the reference LAM-only clustering of high-confidence cells; boundary and normal are auxiliary analyses.",
        "Upstream candidate/state/program correspondence is post-hoc interpretation and is not a Go/No-Go criterion. Any unmatched scVI cluster is retained as a `novel_or_unexplained` candidate for later stages.",
        "The parameter grid is intended to identify a stable parameter interval, not to optimize for a preferred number of clusters. Singleton/small-cluster counts and pairwise partition ARI should be reviewed before treating clusters as biological states.",
        "",
        f"Boundary auxiliary status: {boundary_audit.get('status', 'unknown')}; normal reference status: {normal_audit.get('status', 'unknown')}.",
    ])
    report_path = PROJECT_ROOT / config["outputs"]["stage6_report"]
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"stage6={status} reference={reference_name} grid={len(grid)} report={report_path}")
    return 0 if status != "BLOCKED_INPUT" else 2


if __name__ == "__main__":
    sys.exit(main())
