#!/usr/bin/env python3
"""Decompose local branches around the frozen State 15 anchor.

This stage is geometry-first and diagnostic only.  It reuses the Stage 20
30-nearest-neighbor scope, the frozen State 15 labels, and the Stage 21 score
table.  It does not train scVI, recluster, change the candidate gate, or
redefine any State 1--20 label.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.stats import spearmanr, t
from sklearn.neighbors import NearestNeighbors


PROJECT_ROOT = Path(__file__).resolve().parent
TARGET_STATE = "State_15"
EXPECTED_ANCHOR_CELLS = 200
EXPECTED_MAIN_CELLS = 22261
EXPECTED_CANDIDATE_CELLS = 5378
EXPECTED_BOUNDARY_CELLS = 16883
GRAPH_K = 30
LOCAL_HOPS = [1, 2, 3]
DISTANCE_MATCH_BINS = 5
DISTANCE_SEGMENTS = ["near", "mid", "far"]
LAMCORE_FEATURES = [
    "LAMCORE_full", "LAMCORE_no_gate", "LAMCORE_outside_scVI", "LAMCORE_independent",
]
LAM_FEATURES = [
    *LAMCORE_FEATURES, "CORE1", "CORE2", "CORE3", "melanocytic", "myogenic",
    "VEGFD", "CTSK", "lam_support", "HOX_PBX", "hormone", "ECM", "protease",
    "mTOR_translation",
]
LINEAGE_FEATURES = [
    "T_NK", "macrophage", "AT2", "endothelial", "lymphatic_endothelial",
    "fibroblast", "pericyte_VSMC", "mesothelial", "ciliated",
]
BRANCH_FEATURES = [*LAM_FEATURES, *LINEAGE_FEATURES]


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def load_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_helpers() -> tuple[Any, Any, Any]:
    stage18 = load_module(PROJECT_ROOT / "18_validate_state15_anchor.py", "stage18_for_stage22")
    stage21 = load_module(PROJECT_ROOT / "21_validate_state15_manifold.py", "stage21_for_stage22")
    return stage18, stage21, load_module(PROJECT_ROOT / "20_state15_centered_manifold.py", "stage20_for_stage22")


def read_frozen_inputs(
    distance_path: Path,
    stage21_score_path: Path,
    scvi_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    distances = pd.read_csv(distance_path)
    distances["analysis_cell_id"] = distances["analysis_cell_id"].astype(str)
    distances["current_state"] = distances["current_state"].fillna("").astype(str)
    if distances["analysis_cell_id"].duplicated().any():
        raise ValueError("Stage 20 distance file contains duplicate cell IDs")
    anchor = distances[distances["current_state"].eq(TARGET_STATE)].copy()
    main = distances.copy()
    main["state"] = main["current_state"].replace("", "boundary")
    validation = main[~main["current_state"].eq(TARGET_STATE)].copy()
    scores = pd.read_csv(stage21_score_path)
    scores["analysis_cell_id"] = scores["analysis_cell_id"].astype(str)
    if scores["analysis_cell_id"].duplicated().any():
        raise ValueError("Stage 21 score table contains duplicate cell IDs")
    if len(main) != EXPECTED_MAIN_CELLS or len(anchor) != EXPECTED_ANCHOR_CELLS:
        raise ValueError(f"Unexpected Stage 20 scope: main={len(main)}, anchor={len(anchor)}")
    if len(validation) != EXPECTED_MAIN_CELLS - EXPECTED_ANCHOR_CELLS:
        raise ValueError("Stage 20 non-State15 scope is not frozen at 22,061 cells")
    if len(scores) != len(validation) or set(scores["analysis_cell_id"]) != set(validation["analysis_cell_id"]):
        raise ValueError("Stage 21 score table is not exactly the Stage 20 non-State15 validation object")
    anchor_ids = sorted(anchor["analysis_cell_id"].tolist())
    manifest = {
        "stage20_distance_file": str(distance_path),
        "stage21_score_file": str(stage21_score_path),
        "scvi_artifact": str(scvi_path),
        "latent_key": "X_scVI",
        "graph_k": GRAPH_K,
        "main_cell_count": len(main),
        "candidate_cell_count": int(main["analysis_role"].eq("primary_candidate").sum()),
        "boundary_cell_count": int(main["analysis_role"].eq("boundary").sum()),
        "anchor_state": TARGET_STATE,
        "anchor_cell_count": len(anchor_ids),
        "anchor_cell_id_sha256": hashlib.sha256("\n".join(anchor_ids).encode("utf-8")).hexdigest(),
        "anchor_cell_ids": anchor_ids,
        "no_scvi_training": True,
        "no_reclustering": True,
        "no_candidate_gate_change": True,
        "no_state15_redefinition": True,
    }
    return main, anchor, validation, {"manifest": manifest, "scores": scores}


def load_latent(stage21: Any, scvi_path: Path, ids: pd.Series) -> np.ndarray:
    latent, _ = stage21.load_latent_for_ids(scvi_path, ids)
    if latent.shape[1] != 20:
        raise ValueError(f"Expected 20-dimensional X_scVI, found {latent.shape[1]}")
    return latent


def prepare_scores(
    main: pd.DataFrame,
    anchor: pd.DataFrame,
    validation: pd.DataFrame,
    stage21_scores: pd.DataFrame,
    prepared_path: Path,
    stage18: Any,
    stage21: Any,
    stage20: Any,
    config: dict[str, Any],
    scvi_hvg: set[str],
    block_size: int,
) -> pd.DataFrame:
    validation_scores = validation[["analysis_cell_id", "cell_id", "current_state", "patient", "dataset", "analysis_role", "distance_to_state15", "nearest_state15_distance"]].copy()
    validation_scores = validation_scores.merge(
        stage21_scores.drop(columns=["cell_id", "current_state", "patient", "dataset", "analysis_role", "distance_to_state15", "nearest_state15_distance"], errors="ignore"),
        on="analysis_cell_id",
        how="left",
        validate="one_to_one",
    )
    formal_genes, _ = stage18.resolve_formal_signature(config)
    score_modules, _ = stage21.build_score_modules(stage18, stage20, config, formal_genes, scvi_hvg)
    anchor_selected = anchor[["analysis_cell_id"]].copy()
    anchor_scores, _ = stage21.score_selected_cells(prepared_path, anchor_selected, score_modules, stage18, block_size)
    anchor_scores = anchor[["analysis_cell_id", "cell_id", "current_state", "patient", "dataset", "analysis_role", "distance_to_state15", "nearest_state15_distance"]].merge(
        anchor_scores.drop(columns=["cell_id"], errors="ignore"),
        on="analysis_cell_id",
        how="left",
        validate="one_to_one",
    )
    # Preserve Stage 20 marker-level VEGFD/CTSK values for branch characterization.
    marker_cols = ["analysis_cell_id", "marker_VEGFD", "marker_CTSK"]
    stage20_markers = main[[column for column in marker_cols if column in main.columns]].copy()
    stage20_markers = stage20_markers.rename(columns={"marker_VEGFD": "VEGFD", "marker_CTSK": "CTSK"})
    combined = pd.concat([anchor_scores, validation_scores], ignore_index=True, sort=False)
    combined = combined.drop(columns=["VEGFD", "CTSK"], errors="ignore").merge(stage20_markers, on="analysis_cell_id", how="left")
    combined["current_state"] = combined["current_state"].fillna("").astype(str)
    combined["state"] = combined["current_state"].replace("", "boundary")
    combined["analysis_role"] = combined["analysis_role"].astype(str)
    combined["patient"] = combined["patient"].astype(str)
    combined["dataset"] = combined["dataset"].astype(str)
    combined = combined.drop_duplicates("analysis_cell_id").reset_index(drop=True)
    if len(combined) != len(main):
        raise ValueError("Combined Stage 20/21 score table does not cover the frozen main cohort")
    combined = combined.set_index("analysis_cell_id").loc[main["analysis_cell_id"].tolist()].reset_index()
    return combined


def build_local_graph(latent: np.ndarray, table: pd.DataFrame) -> tuple[sparse.csr_matrix, np.ndarray]:
    if len(latent) != len(table):
        raise ValueError("Latent and table row counts differ")
    knn = NearestNeighbors(n_neighbors=min(GRAPH_K + 1, len(latent)), metric="euclidean", n_jobs=1).fit(latent)
    _, directed_indices = knn.kneighbors(latent)
    rows = np.repeat(np.arange(len(latent)), directed_indices.shape[1] - 1)
    cols = directed_indices[:, 1:].ravel()
    directed = sparse.csr_matrix((np.ones(len(rows), dtype=np.uint8), (rows, cols)), shape=(len(latent), len(latent)))
    undirected = directed.maximum(directed.T).tocsr()
    return undirected, directed_indices


def graph_hops(adjacency: sparse.csr_matrix, anchor_mask: np.ndarray, max_hop: int = 3) -> np.ndarray:
    hops = np.full(adjacency.shape[0], -1, dtype=np.int16)
    frontier = np.flatnonzero(anchor_mask)
    hops[frontier] = 0
    for hop in range(1, max_hop + 1):
        if len(frontier) == 0:
            break
        neighbors = np.unique(adjacency[frontier].indices)
        frontier = neighbors[hops[neighbors] < 0]
        hops[frontier] = hop
    return hops


def local_graph_table(table: pd.DataFrame, adjacency: sparse.csr_matrix, directed_indices: np.ndarray, hops: np.ndarray) -> pd.DataFrame:
    anchor_indices = np.flatnonzero(table["state"].astype(str).eq(TARGET_STATE).to_numpy())
    state15_counts = np.asarray(adjacency[:, anchor_indices].getnnz(axis=1)).ravel()
    graph_degree = np.diff(adjacency.indptr)
    output = table[["cell_id", "analysis_cell_id", "current_state", "state", "analysis_role", "patient", "dataset", "distance_to_state15", "state15_neighbor_count_k30", "state15_neighbor_fraction"]].copy()
    output = output.rename(
        columns={
            "current_state": "current_state",
            "distance_to_state15": "latent_distance_to_State15",
            "state15_neighbor_count_k30": "stage20_directed_State15_neighbor_count",
            "state15_neighbor_fraction": "stage20_directed_State15_neighbor_fraction",
        }
    )
    output["number_of_State15_neighbors"] = state15_counts
    output["fraction_of_State15_neighbors"] = state15_counts / np.maximum(graph_degree, 1)
    output["min_graph_hop_to_State15"] = hops
    output["within_1_3_hops"] = np.isin(hops, [1, 2, 3])
    output["is_State15_anchor"] = output["state"].eq(TARGET_STATE)
    output["directed_knn_neighbor_count"] = directed_indices.shape[1] - 1
    return output


def state_connectivity(table: pd.DataFrame, adjacency: sparse.csr_matrix, hops: np.ndarray) -> pd.DataFrame:
    state = table["state"].astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    state15_indices = np.flatnonzero(state == TARGET_STATE)
    anchor_neighbors = np.concatenate([adjacency[index].indices for index in state15_indices])
    external_targets = state[anchor_neighbors]
    state_counts = pd.Series(state).value_counts().to_dict()
    for target_state, total in sorted(state_counts.items()):
        if target_state == TARGET_STATE:
            continue
        edge_count = int(np.sum(external_targets == target_state))
        direct_cells = int(np.sum((state == target_state) & (hops == 1)))
        hop_counts = [int(np.sum((state == target_state) & (hops == hop))) for hop in [1, 2, 3]]
        target_indices = np.flatnonzero(state == target_state)
        direct_target_indices = np.flatnonzero((state == target_state) & (hops == 1))
        # Branch eligibility must use the direct State15 neighborhood.  Keep
        # whole-state coverage separately for descriptive context.
        patient_count = int(table.iloc[direct_target_indices]["patient"].astype(str).nunique()) if len(direct_target_indices) else 0
        dataset_count = int(table.iloc[direct_target_indices]["dataset"].astype(str).nunique()) if len(direct_target_indices) else 0
        state_patient_count = int(table.iloc[target_indices]["patient"].astype(str).nunique()) if len(target_indices) else 0
        state_dataset_count = int(table.iloc[target_indices]["dataset"].astype(str).nunique()) if len(target_indices) else 0
        rows.append(
            {
                "state": target_state,
                "state_total_cells": int(total),
                "hop_1_cells": hop_counts[0],
                "hop_2_cells": hop_counts[1],
                "hop_3_cells": hop_counts[2],
                "state15_edge_count": edge_count,
                "state15_edge_count_over_state_cells": float(edge_count / max(int(total), 1)),
                "state15_direct_neighbor_cell_count": direct_cells,
                "patient_count": patient_count,
                "dataset_count": dataset_count,
                "state_patient_count": state_patient_count,
                "state_dataset_count": state_dataset_count,
            }
        )
    return pd.DataFrame(rows).sort_values("state").reset_index(drop=True)


def select_branches(connectivity: pd.DataFrame) -> pd.DataFrame:
    if len(connectivity) == 0:
        return pd.DataFrame(columns=["branch_id", "source_state", "1hop_cells", "2hop_cells", "3hop_cells", "patient_count", "dataset_count", "state_patient_count", "state_dataset_count", "connectivity"])
    eligible = connectivity[(connectivity["state"] != "boundary") & (connectivity["hop_1_cells"] >= 10) & (connectivity["patient_count"] >= 2)].copy()
    eligible = eligible.sort_values(["hop_1_cells", "state"], ascending=[False, True]).reset_index(drop=True)
    eligible["branch_id"] = [f"branch_{i + 1:02d}" for i in range(len(eligible))]
    return eligible.rename(
        columns={
            "state": "source_state",
            "hop_1_cells": "1hop_cells",
            "hop_2_cells": "2hop_cells",
            "hop_3_cells": "3hop_cells",
            "state15_edge_count_over_state_cells": "connectivity",
        }
    )[["branch_id", "source_state", "1hop_cells", "2hop_cells", "3hop_cells", "patient_count", "dataset_count", "state_patient_count", "state_dataset_count", "connectivity"]]


def assign_segments(values: pd.Series) -> pd.Series:
    if len(values) == 0:
        return pd.Series([], index=values.index, dtype="string")
    ranks = values.rank(method="first", pct=True).to_numpy(dtype=float)
    indices = np.clip(np.ceil(ranks * 3).astype(int) - 1, 0, 2)
    return pd.Series(np.asarray(DISTANCE_SEGMENTS, dtype=object)[indices], index=values.index, dtype="string")


def numeric_slope(distance: pd.Series, score: pd.Series, patient: pd.Series | None = None) -> dict[str, float]:
    frame = pd.DataFrame({"distance": pd.to_numeric(distance, errors="coerce"), "score": pd.to_numeric(score, errors="coerce")})
    if patient is not None:
        frame["patient"] = patient.astype(str)
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < 4 or frame["distance"].nunique() < 2:
        return {"n": len(frame), "slope": np.nan, "ci95_low": np.nan, "ci95_high": np.nan, "pvalue": np.nan}
    if patient is not None and frame["patient"].nunique() > 1:
        dummy = pd.get_dummies(frame["patient"], drop_first=True, dtype=float).to_numpy()
        design = np.column_stack([np.ones(len(frame)), frame["distance"].to_numpy(dtype=float), dummy])
    else:
        design = np.column_stack([np.ones(len(frame)), frame["distance"].to_numpy(dtype=float)])
    response = frame["score"].to_numpy(dtype=float)
    beta, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
    residual = response - design @ beta
    df_resid = len(response) - int(rank)
    slope = float(beta[1])
    if df_resid <= 0:
        return {"n": len(frame), "slope": slope, "ci95_low": np.nan, "ci95_high": np.nan, "pvalue": np.nan}
    sigma2 = float(np.dot(residual, residual) / df_resid)
    covariance = sigma2 * np.linalg.pinv(design.T @ design)
    se = float(np.sqrt(max(covariance[1, 1], 0.0)))
    if se == 0:
        return {"n": len(frame), "slope": slope, "ci95_low": slope, "ci95_high": slope, "pvalue": 0.0 if slope else np.nan}
    critical = float(t.ppf(0.975, df_resid))
    return {
        "n": len(frame),
        "slope": slope,
        "ci95_low": slope - critical * se,
        "ci95_high": slope + critical * se,
        "pvalue": float(2 * t.sf(abs(slope / se), df_resid)),
    }


def score_summary(table: pd.DataFrame, group_field: str, group_value: str, features: list[str]) -> dict[str, Any]:
    sub = table[table[group_field].astype(str).eq(str(group_value))]
    row: dict[str, Any] = {group_field: str(group_value), "n_cells": len(sub)}
    for feature in features:
        values = pd.to_numeric(sub[feature], errors="coerce").dropna() if feature in sub else pd.Series(dtype=float)
        row[f"{feature}_median"] = float(values.median()) if len(values) else np.nan
        row[f"{feature}_q25"] = float(values.quantile(0.25)) if len(values) else np.nan
        row[f"{feature}_q75"] = float(values.quantile(0.75)) if len(values) else np.nan
    return row


def local_branch_scope(table: pd.DataFrame, source_state: str) -> pd.DataFrame:
    """Return a branch's fixed local 1–3-hop scope."""
    return table[
        table["state"].eq(source_state)
        & table["min_graph_hop_to_State15"].isin(LOCAL_HOPS)
    ].copy()


def branch_gradient(branch_id: str, source_state: str, table: pd.DataFrame, features: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sub = local_branch_scope(table, source_state)
    sub["distance_segment"] = assign_segments(sub["latent_distance_to_State15"])
    segment_rows: list[dict[str, Any]] = []
    for segment in DISTANCE_SEGMENTS:
        part = sub[sub["distance_segment"].eq(segment)]
        row: dict[str, Any] = {"row_type": "segment", "branch_id": branch_id, "source_state": source_state, "distance_segment": segment, "n_cells": len(part)}
        for feature in features:
            values = pd.to_numeric(part[feature], errors="coerce").dropna() if feature in part else pd.Series(dtype=float)
            row[f"{feature}_median"] = float(values.median()) if len(values) else np.nan
            row[f"{feature}_q25"] = float(values.quantile(0.25)) if len(values) else np.nan
            row[f"{feature}_q75"] = float(values.quantile(0.75)) if len(values) else np.nan
        segment_rows.append(row)
    segment = pd.DataFrame(segment_rows)
    model_rows: list[dict[str, Any]] = []
    for metric, distance in [("latent_distance", sub["latent_distance_to_State15"]), ("graph_hop", sub["min_graph_hop_to_State15"])]:
        for feature in features:
            valid = sub[["patient", feature]].copy()
            valid["distance"] = pd.to_numeric(distance, errors="coerce")
            valid[feature] = pd.to_numeric(valid[feature], errors="coerce")
            valid = valid.replace([np.inf, -np.inf], np.nan).dropna()
            if len(valid) >= 3 and valid["distance"].nunique() > 1:
                rho, rho_p = spearmanr(valid["distance"], valid[feature])
                fit = numeric_slope(valid["distance"], valid[feature], valid["patient"])
            else:
                rho, rho_p = np.nan, np.nan
                fit = {"n": len(valid), "slope": np.nan, "ci95_low": np.nan, "ci95_high": np.nan, "pvalue": np.nan}
            model_rows.append(
                {
                    "row_type": "model",
                    "branch_id": branch_id,
                    "source_state": source_state,
                    "distance_metric": metric,
                    "score_name": feature,
                    "n_cells": len(valid),
                    "n_patients": int(valid["patient"].nunique()),
                    "patient_adjusted_slope": fit["slope"],
                    "ci95_low": fit["ci95_low"],
                    "ci95_high": fit["ci95_high"],
                    "pvalue": fit["pvalue"],
                    "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                    "spearman_pvalue": float(rho_p) if np.isfinite(rho_p) else np.nan,
                }
            )
    models = pd.DataFrame(model_rows)
    positions = sub[["analysis_cell_id", "cell_id", "state", "patient", "dataset", "latent_distance_to_State15", "min_graph_hop_to_State15", "fraction_of_State15_neighbors", "distance_segment", *[feature for feature in features if feature in sub.columns]]].sort_values("latent_distance_to_State15")
    return segment, models, positions


def patient_branch_consistency(branch_id: str, source_state: str, table: pd.DataFrame, features: list[str], min_cells: int = 10) -> pd.DataFrame:
    branch = local_branch_scope(table, source_state)
    rows: list[dict[str, Any]] = []
    for patient, sub in branch.groupby("patient", observed=True):
        if len(sub) < min_cells:
            continue
        sub = sub.copy()
        sub["distance_segment"] = assign_segments(sub["latent_distance_to_State15"])
        near = sub[sub["distance_segment"].eq("near")]
        far = sub[sub["distance_segment"].eq("far")]
        patient_fit = numeric_slope(sub["latent_distance_to_State15"], sub["LAMCORE_independent"])
        row: dict[str, Any] = {
            "branch_id": branch_id,
            "source_state": source_state,
            "patient": str(patient),
            "n_branch_cells": len(sub),
            "near_cells": len(near),
            "far_cells": len(far),
            "patient_LAMCORE_slope": patient_fit["slope"],
            "patient_LAMCORE_slope_ci95_low": patient_fit["ci95_low"],
            "patient_LAMCORE_slope_ci95_high": patient_fit["ci95_high"],
            "patient_LAMCORE_slope_pvalue": patient_fit["pvalue"],
            "patient_LAMCORE_slope_direction": (
                "decrease" if np.isfinite(patient_fit["slope"]) and patient_fit["slope"] < 0
                else "increase_or_flat" if np.isfinite(patient_fit["slope"])
                else "not_estimable"
            ),
            "patient_distance_range": float(sub["latent_distance_to_State15"].max() - sub["latent_distance_to_State15"].min()),
        }
        for feature in features:
            near_value = pd.to_numeric(near[feature], errors="coerce").median() if feature in near else np.nan
            far_value = pd.to_numeric(far[feature], errors="coerce").median() if feature in far else np.nan
            row[f"{feature}_near_median"] = float(near_value) if pd.notna(near_value) else np.nan
            row[f"{feature}_far_median"] = float(far_value) if pd.notna(far_value) else np.nan
            row[f"{feature}_far_minus_near"] = float(far_value - near_value) if pd.notna(near_value) and pd.notna(far_value) else np.nan
            if feature == "LAMCORE_independent":
                row["LAMCORE_direction"] = "decrease" if pd.notna(near_value) and pd.notna(far_value) and far_value < near_value else "increase_or_flat" if pd.notna(near_value) and pd.notna(far_value) else "not_estimable"
            if feature in LINEAGE_FEATURES:
                row[f"{feature}_direction"] = "increase" if pd.notna(near_value) and pd.notna(far_value) and far_value > near_value else "decrease_or_flat" if pd.notna(near_value) and pd.notna(far_value) else "not_estimable"
        rows.append(row)
    return pd.DataFrame(rows)


def patient_lopo_robustness(branches: pd.DataFrame, table: pd.DataFrame) -> pd.DataFrame:
    """Fit the branch slope after leaving each represented patient out."""
    rows: list[dict[str, Any]] = []
    for _, branch in branches.iterrows():
        branch_id = str(branch["branch_id"])
        source_state = str(branch["source_state"])
        local = local_branch_scope(table, source_state)
        patients = sorted(local["patient"].astype(str).unique())
        for omitted_patient in patients:
            remaining = local[~local["patient"].astype(str).eq(omitted_patient)].copy()
            fit = numeric_slope(
                remaining["latent_distance_to_State15"],
                remaining["LAMCORE_independent"],
                remaining["patient"],
            )
            rows.append(
                {
                    "branch_id": branch_id,
                    "source_state": source_state,
                    "omitted_patient": omitted_patient,
                    "n_remaining_cells": len(remaining),
                    "n_remaining_patients": int(remaining["patient"].astype(str).nunique()),
                    "LOPO_LAMCORE_slope": fit["slope"],
                    "LOPO_ci95_low": fit["ci95_low"],
                    "LOPO_ci95_high": fit["ci95_high"],
                    "LOPO_pvalue": fit["pvalue"],
                    "LOPO_slope_direction": (
                        "decrease" if np.isfinite(fit["slope"]) and fit["slope"] < 0
                        else "increase_or_flat" if np.isfinite(fit["slope"])
                        else "not_estimable"
                    ),
                }
            )
    return pd.DataFrame(rows)


def boundary_assignment(table: pd.DataFrame, adjacency: sparse.csr_matrix, branches: pd.DataFrame) -> pd.DataFrame:
    # Only local boundary cells are eligible for branch projection. Farther
    # boundary cells remain in the frozen cohort but have no local direction.
    boundary = table[
        table["state"].eq("boundary")
        & table["min_graph_hop_to_State15"].isin(LOCAL_HOPS)
    ].copy().reset_index(drop=True)
    branch_states = branches["source_state"].astype(str).tolist()
    state = table["state"].astype(str).to_numpy()
    rows: list[dict[str, Any]] = []
    for row_index, row in boundary.iterrows():
        original_index = table.index[table["analysis_cell_id"].eq(row["analysis_cell_id"])][0]
        neighbors = state[adjacency[original_index].indices]
        counts = {branch_state: int(np.sum(neighbors == branch_state)) for branch_state in branch_states}
        ordered = sorted(counts.items(), key=lambda item: (-item[1], item[0]))
        best_state, best_count = ordered[0] if ordered else ("", 0)
        second_count = ordered[1][1] if len(ordered) > 1 else 0
        assignment = best_state if best_count >= 2 and best_count > second_count else "unresolved"
        record = row.to_dict()
        record.update(counts)
        record.update(
            {
                "nearest_branch": best_state,
                "nearest_branch_neighbor_count": best_count,
                "second_branch_neighbor_count": second_count,
                "branch_assignment": assignment,
                "branch_assignment_rule": "best direct branch-state neighbor count >=2 and strictly greater than second",
            }
        )
        rows.append(record)
    return pd.DataFrame(rows)


def boundary_extension(boundary: pd.DataFrame, branches: pd.DataFrame, table: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    anchor = table[table["state"].eq(TARGET_STATE)]
    rows.append({"branch_id": "all", "source_state": TARGET_STATE, "extension_stage": "State15_anchor", **score_summary(anchor, "state", TARGET_STATE, features)})
    for _, branch in branches.iterrows():
        branch_id = str(branch["branch_id"])
        state = str(branch["source_state"])
        branch_table = local_branch_scope(table, state)
        branch_table["distance_segment"] = assign_segments(branch_table["latent_distance_to_State15"])
        for segment in DISTANCE_SEGMENTS:
            part = branch_table[branch_table["distance_segment"].eq(segment)]
            rows.append({"branch_id": branch_id, "source_state": state, "extension_stage": f"{state}_{segment}", **score_summary(part, "state", state, features)})
        extension = boundary[boundary["branch_assignment"].eq(state)].copy()
        rows.append({"branch_id": branch_id, "source_state": state, "extension_stage": "boundary_assigned_extension", **score_summary(extension, "branch_assignment", state, features)})
    return pd.DataFrame(rows)


def empirical_tail_probabilities(null_values: np.ndarray, real_slope: float) -> dict[str, float]:
    """Compute direct empirical tails without assuming a zero-centered null."""
    if len(null_values) == 0 or not np.isfinite(real_slope):
        return {"empirical_left_p": np.nan, "empirical_right_p": np.nan, "empirical_two_sided_p": np.nan}
    n = len(null_values)
    left = float((1 + np.sum(null_values <= real_slope)) / (n + 1))
    right = float((1 + np.sum(null_values >= real_slope)) / (n + 1))
    return {
        "empirical_left_p": left,
        "empirical_right_p": right,
        "empirical_two_sided_p": float(min(1.0, 2.0 * min(left, right))),
    }


def bh_adjust(values: pd.Series) -> np.ndarray:
    """Benjamini–Hochberg adjustment, preserving NaN for unavailable tests."""
    p = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    q = np.full(len(p), np.nan, dtype=float)
    valid = np.isfinite(p)
    if not valid.any():
        return q
    indices = np.flatnonzero(valid)
    order = indices[np.argsort(p[valid], kind="mergesort")]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1, dtype=float)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(ranked, 1.0)
    return q


def _distance_bins(table: pd.DataFrame, n_bins: int = DISTANCE_MATCH_BINS) -> pd.Series:
    distances = pd.to_numeric(table["distance_to_state15"], errors="coerce")
    ranks = distances.rank(method="first", pct=False)
    bins = np.floor((ranks - 1) * n_bins / max(len(table), 1)).astype("Int64")
    return bins.clip(lower=0, upper=n_bins - 1)


def _sample_distance_matched(
    real: pd.DataFrame,
    pool: pd.DataFrame,
    rng: np.random.Generator,
) -> tuple[list[int], int, float, set[str]]:
    """Sample null cells matching patient×dataset, distance-bin and count.

    Distance bins preserve the coarse local geometry. If a joint stratum/bin
    is absent, the method falls back to the same joint stratum, then the same
    patient, choosing cells closest to the target distances and recording the
    relaxation level.
    """
    real = real.copy()
    pool = pool.copy()
    real["_key"] = real["patient"].astype(str) + "||" + real["dataset"].astype(str)
    pool["_key"] = pool["patient"].astype(str) + "||" + pool["dataset"].astype(str)
    real["_distance_bin"] = real["_distance_bin"].astype(int)
    pool["_distance_bin"] = pool["_distance_bin"].astype(int)
    by_key_bin = {
        (str(key), int(distance_bin)): values.index.to_numpy(dtype=int)
        for (key, distance_bin), values in pool.groupby(["_key", "_distance_bin"], observed=True)
    }
    by_key = {str(key): values.index.to_numpy(dtype=int) for key, values in pool.groupby("_key", observed=True)}
    by_patient = {str(key): values.index.to_numpy(dtype=int) for key, values in pool.groupby("patient", observed=True)}
    sampled: list[int] = []
    replacements = 0
    abs_distance_differences: list[float] = []
    levels: set[str] = set()
    grouped = real.groupby(["_key", "_distance_bin"], observed=True, sort=False)
    for (key, distance_bin), target in grouped:
        key = str(key)
        distance_bin = int(distance_bin)
        target_distances = pd.to_numeric(target["distance_to_state15"], errors="coerce").to_numpy(dtype=float)
        eligible = by_key_bin.get((key, distance_bin))
        level = "patient_dataset_distance_bin"
        if eligible is None or len(eligible) == 0:
            eligible = by_key.get(key)
            level = "patient_dataset_nearest_distance"
        if eligible is None or len(eligible) == 0:
            eligible = by_patient.get(key.split("||", 1)[0])
            level = "patient_nearest_distance"
        if eligible is None or len(eligible) == 0:
            raise ValueError(f"No matched null pool for patient/dataset stratum {key}")
        replace = len(eligible) < len(target)
        if level == "patient_dataset_distance_bin":
            choices = rng.choice(eligible, size=len(target), replace=replace)
            chosen_distances = pd.to_numeric(pool.loc[choices, "distance_to_state15"], errors="coerce").to_numpy(dtype=float)
        else:
            choices_list: list[int] = []
            chosen_distances_list: list[float] = []
            for target_distance in target_distances:
                candidate_distances = pd.to_numeric(pool.loc[eligible, "distance_to_state15"], errors="coerce").to_numpy(dtype=float)
                delta = np.abs(candidate_distances - target_distance)
                finite = np.isfinite(delta)
                if not finite.any():
                    choice_pool = eligible
                else:
                    finite_indices = np.flatnonzero(finite)
                    nearest_order = finite_indices[np.argsort(delta[finite], kind="mergesort")]
                    choice_pool = eligible[nearest_order[: min(10, len(nearest_order))]]
                choice = int(rng.choice(choice_pool))
                choices_list.append(choice)
                chosen_distances_list.append(float(pool.loc[choice, "distance_to_state15"]))
            choices = np.asarray(choices_list, dtype=int)
            chosen_distances = np.asarray(chosen_distances_list, dtype=float)
        sampled.extend(int(value) for value in choices)
        replacements += int(len(target)) if replace else 0
        abs_distance_differences.extend(np.abs(chosen_distances - target_distances).tolist())
        levels.add(level)
    mean_abs_delta = float(np.nanmean(abs_distance_differences)) if abs_distance_differences else np.nan
    return sampled, replacements, mean_abs_delta, levels


def matched_branch_null(
    branches: pd.DataFrame,
    table: pd.DataFrame,
    reps: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    rows: list[dict[str, Any]] = []
    summaries: dict[str, Any] = {}
    # Real and null branches share exactly this local scope.
    local = table[table["min_graph_hop_to_State15"].isin(LOCAL_HOPS)].copy()
    local["_distance_bin"] = _distance_bins(local)
    for _, branch in branches.iterrows():
        branch_id = str(branch["branch_id"])
        source_state = str(branch["source_state"])
        real = local[local["state"].eq(source_state)].copy()
        pool = local[(local["state"].ne(source_state)) & (local["state"].ne(TARGET_STATE))].copy()
        real_keys = real["patient"].astype(str) + "||" + real["dataset"].astype(str)
        pool_keys = pool["patient"].astype(str) + "||" + pool["dataset"].astype(str)
        requested = real_keys.value_counts().to_dict()
        pool_key_set = set(pool_keys)
        missing = sorted(set(requested).difference(pool_key_set))
        by_patient = set(pool["patient"].astype(str))
        missing_patient = sorted({key.split("||", 1)[0] for key in missing if key.split("||", 1)[0] not in by_patient})
        if missing_patient:
            summaries[branch_id] = {"status": "insufficient_matched_null_pool", "missing_strata": missing, "missing_patients": missing_patient}
            continue
        real_valid = real[["patient", "distance_to_state15", "LAMCORE_independent"]].copy()
        real_fit = numeric_slope(real_valid["distance_to_state15"], real_valid["LAMCORE_independent"], real_valid["patient"])
        null_slopes: list[float] = []
        branch_summary: list[dict[str, Any]] = []
        for replicate in range(reps):
            sampled, replacements, mean_abs_delta, levels = _sample_distance_matched(real, pool, rng)
            fake = pool.loc[sampled]
            fit = numeric_slope(fake["distance_to_state15"], fake["LAMCORE_independent"], fake["patient"])
            null_slopes.append(float(fit["slope"]) if np.isfinite(fit["slope"]) else np.nan)
            branch_summary.append(
                {
                    "branch_id": branch_id,
                    "source_state": source_state,
                    "replicate": replicate,
                    "n_real_branch_cells": len(real),
                    "n_fake_branch_cells": len(fake),
                    "matched_strata": len(requested),
                    "sampling_with_replacement_cells": replacements,
                    "distance_match_bins": DISTANCE_MATCH_BINS,
                    "distance_match_mean_abs_delta": mean_abs_delta,
                    "distance_match_levels": ";".join(sorted(levels)),
                    "real_slope": real_fit["slope"],
                    "null_slope": fit["slope"],
                    "distance_metric": "Stage20 distance_to_state15",
                    "scope": "real and null: non-State15 local 1–3-hop cells",
                }
            )
        null_values = np.asarray([value for value in null_slopes if np.isfinite(value)], dtype=float)
        real_slope = float(real_fit["slope"]) if np.isfinite(real_fit["slope"]) else np.nan
        tails = empirical_tail_probabilities(null_values, real_slope)
        null_median = float(np.median(null_values)) if len(null_values) else np.nan
        null_q05 = float(np.quantile(null_values, 0.05)) if len(null_values) else np.nan
        null_q95 = float(np.quantile(null_values, 0.95)) if len(null_values) else np.nan
        for row in branch_summary:
            row["null_median_slope"] = null_median
            row["null_q05_slope"] = null_q05
            row["null_q95_slope"] = null_q95
            row.update(tails)
        rows.extend(branch_summary)
        summaries[branch_id] = {
            "status": "available",
            "real_slope": real_slope,
            "null_median_slope": null_median,
            "null_q05_slope": null_q05,
            "null_q95_slope": null_q95,
            **tails,
            "repetitions": len(null_values),
            "matched_strata": requested,
            "scope": "real and null: non-State15 local 1–3-hop cells",
            "distance_match_bins": DISTANCE_MATCH_BINS,
            "distance_match_method": "patient×dataset plus global distance-bin matching with nearest-distance fallback",
        }
    return pd.DataFrame(rows), summaries


def evidence_summary(
    branches: pd.DataFrame,
    all_models: pd.DataFrame,
    patient_consistency: pd.DataFrame,
    patient_lopo: pd.DataFrame,
    null_manifest: dict[str, Any],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for _, branch in branches.iterrows():
        branch_id = str(branch["branch_id"])
        source_state = str(branch["source_state"])
        models = all_models[(all_models["branch_id"].eq(branch_id)) & all_models["row_type"].eq("model")]
        independent = models[(models["score_name"].eq("LAMCORE_independent")) & models["distance_metric"].eq("latent_distance")]
        core3 = models[(models["score_name"].eq("CORE3")) & models["distance_metric"].eq("latent_distance")]
        lineages = models[(models["distance_metric"].eq("latent_distance")) & models["score_name"].isin(LINEAGE_FEATURES)].copy()
        best_lineage = ""
        best_lineage_slope = np.nan
        if len(lineages):
            lineages = lineages.sort_values("patient_adjusted_slope", ascending=False)
            best_lineage = str(lineages.iloc[0]["score_name"])
            best_lineage_slope = float(lineages.iloc[0]["patient_adjusted_slope"])
        patient = patient_consistency[patient_consistency["branch_id"].eq(branch_id)]
        lam_direction_fraction = float((patient["LAMCORE_direction"].eq("decrease")).mean()) if len(patient) else np.nan
        patient_slopes = pd.to_numeric(patient.get("patient_LAMCORE_slope", pd.Series(dtype=float)), errors="coerce").dropna()
        patient_slope_negative_fraction = float((patient_slopes < 0).mean()) if len(patient_slopes) else np.nan
        lopo = patient_lopo[patient_lopo["branch_id"].eq(branch_id)] if len(patient_lopo) else pd.DataFrame()
        lopo_slopes = pd.to_numeric(lopo.get("LOPO_LAMCORE_slope", pd.Series(dtype=float)), errors="coerce").dropna()
        lopo_negative_fraction = float((lopo_slopes < 0).mean()) if len(lopo_slopes) else np.nan
        null = null_manifest.get(branch_id, {})
        lam_slope = float(independent["patient_adjusted_slope"].iloc[0]) if len(independent) and pd.notna(independent["patient_adjusted_slope"].iloc[0]) else np.nan
        core3_slope = float(core3["patient_adjusted_slope"].iloc[0]) if len(core3) and pd.notna(core3["patient_adjusted_slope"].iloc[0]) else np.nan
        null_p = float(null.get("empirical_two_sided_p", np.nan)) if null else np.nan
        null_q = np.nan
        lam_like = bool(
            np.isfinite(lam_slope)
            and lam_slope < 0
            and np.isfinite(null_p)
            and null_p <= 0.05
            and (not np.isfinite(lam_direction_fraction) or lam_direction_fraction > 0.5)
            and (not np.isfinite(patient_slope_negative_fraction) or patient_slope_negative_fraction > 0.5)
        )
        transition = bool(lam_like and np.isfinite(best_lineage_slope) and best_lineage_slope > 0)
        if transition:
            label = "LAM_to_lineage_transition_candidate"
        elif lam_like:
            label = "LAM_like_branch_candidate"
        elif np.isfinite(best_lineage_slope) and best_lineage_slope > 0:
            label = "ordinary_lineage_adjacency"
        else:
            label = "ambiguous_branch"
        rows.append(
            {
                "branch_id": branch_id,
                "source_state": source_state,
                "1hop_cells": branch["1hop_cells"],
                "patient_count": branch["patient_count"],
                "dataset_count": branch["dataset_count"],
                "LAMCORE_independent_latent_slope": lam_slope,
                "CORE3_latent_slope": core3_slope,
                "patient_LAMCORE_decrease_fraction": lam_direction_fraction,
                "patient_LAMCORE_slope_median": float(patient_slopes.median()) if len(patient_slopes) else np.nan,
                "patient_LAMCORE_slope_negative_fraction": patient_slope_negative_fraction,
                "patient_LAMCORE_slope_n": int(len(patient_slopes)),
                "LOPO_LAMCORE_slope_negative_fraction": lopo_negative_fraction,
                "LOPO_LAMCORE_slope_n": int(len(lopo_slopes)),
                "dominant_competing_lineage": best_lineage,
                "dominant_competing_lineage_slope": best_lineage_slope,
                "matched_null_empirical_left_p": float(null.get("empirical_left_p", np.nan)) if null else np.nan,
                "matched_null_empirical_right_p": float(null.get("empirical_right_p", np.nan)) if null else np.nan,
                "matched_null_empirical_p": null_p,
                "matched_null_q_value": null_q,
                "evidence_label": label,
            }
        )
    output = pd.DataFrame(rows)
    if len(output):
        output["matched_null_q_value"] = bh_adjust(output["matched_null_empirical_p"])
        # The corrected multiple-testing result, not an individual raw p,
        # controls whether a branch can receive a LAM-preserving label.
        for index, row in output.iterrows():
            q_value = row["matched_null_q_value"]
            raw_label = str(row["evidence_label"])
            lam_like = bool(
                np.isfinite(row["LAMCORE_independent_latent_slope"])
                and row["LAMCORE_independent_latent_slope"] < 0
                and np.isfinite(q_value)
                and q_value <= 0.05
                and (not np.isfinite(row["patient_LAMCORE_slope_negative_fraction"]) or row["patient_LAMCORE_slope_negative_fraction"] > 0.5)
            )
            lineage_positive = np.isfinite(row["dominant_competing_lineage_slope"]) and row["dominant_competing_lineage_slope"] > 0
            corrected_label = (
                "LAM_to_lineage_transition_candidate" if lam_like and lineage_positive
                else "LAM_like_branch_candidate" if lam_like
                else "ordinary_lineage_adjacency" if lineage_positive
                else "ambiguous_branch"
            )
            output.at[index, "evidence_label"] = corrected_label
            output.at[index, "raw_evidence_label_before_fdr"] = raw_label
    return output


def checkpoint(summary: pd.DataFrame) -> tuple[str, str]:
    if len(summary) == 0:
        return "no_reproducible_local_branch", "No existing State 1–20 state met the fixed direct-connectivity rule of at least 10 one-hop cells from at least 2 patients."
    labels = set(summary["evidence_label"].astype(str))
    if "LAM_to_lineage_transition_candidate" in labels:
        return "local_branched_lam_manifold_candidate", "At least one local branch retains independent LAM identity while moving toward a distinct lineage direction; interpret as a local branched candidate, not a new state."
    if "LAM_like_branch_candidate" in labels:
        return "local_lam_like_branch_candidates", "At least one local branch retains independent LAM identity, but no clear competing-lineage transition was established."
    if "ordinary_lineage_adjacency" in labels:
        return "ordinary_lineage_adjacency_dominates", "Directly connected branches are better described as ordinary lineage adjacency than LAM-preserving branches."
    return "ambiguous_local_branch_structure", "Local branches exist, but the current evidence does not distinguish LAM-preserving structure from ordinary adjacency."


def write_report(output_dir: Path, manifest: dict[str, Any], connectivity: pd.DataFrame, branches: pd.DataFrame, evidence: pd.DataFrame, null_manifest: dict[str, Any], cp: str, interpretation: str) -> None:
    lines = [
        "# Stage 22：State 15 局部分支分解",
        "",
        "本阶段固定使用 Stage 20/21 的 22,261 个细胞、State 15 的 200 个 anchor、既有 `X_scVI` 和 Stage 21 scores；不重训 scVI、不重新聚类、不修改 candidate gate、不修改 State 1–20 标签。",
        "",
        "## Frozen scope",
        "",
        f"- Main cells: {manifest['main_cell_count']} ({manifest['candidate_cell_count']} candidates + {manifest['boundary_cell_count']} boundary).",
        f"- State 15 anchor: {manifest['anchor_cell_count']} cells; ID SHA-256 `{manifest['anchor_cell_id_sha256']}`.",
        f"- Local graph: undirected form of the same `X_scVI` k={manifest['graph_k']} neighbor scope; no Leiden.",
        "",
        "## State 15 connectivity and branch selection",
        "",
        f"- Fixed branch candidates selected: {len(branches)}.",
        "- Selection rule: external existing state, at least 10 one-hop cells, and at least 2 patients among those one-hop cells; boundary is not promoted to a branch.",
        "详见 `state15_state_connectivity.csv` 和 `branch_candidates.csv`。",
        "",
        "## Branch evidence",
        "",
        evidence.to_string(index=False) if len(evidence) else "No branch candidate was selected.",
        "",
        "## Boundary",
        "",
        "Boundary cells within 1–3 hops are assigned only to a local direction when the direct branch-neighbor count is at least 2 and strictly exceeds the second branch; unresolved cells remain unresolved and no new LAM label is produced.",
        "",
        "## Matched null",
        "",
        f"- Null repetitions per available branch: {manifest['branch_null_repetitions']}.",
        "- Real and null scopes are identical: non-State15 local 1–3-hop cells; null cells match patient×dataset, cell count and five-bin local distance structure.",
        "- Empirical p-values use direct left/right tails of the observed null distribution; the two-sided p is `2*min(left,right)`, with no zero-centered or symmetry assumption.",
        "- Benjamini–Hochberg q-values are computed across all selected branches with available null tests; labels use q rather than raw p alone.",
        "- Per-patient slopes and leave-one-patient-out fits are recorded in `branch_patient_consistency.csv` and `branch_patient_lopo.csv`.",
        "",
        "## Stage 22 checkpoint",
        "",
        f"- `{cp}`",
        f"- {interpretation}",
        "",
        "## Outputs",
        "",
        *[f"- {path.name}" for path in sorted(output_dir.iterdir()) if path.is_file()],
    ]
    (output_dir / "stage22_local_branch_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/state_modeling.yaml"))
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "results/stage22"))
    parser.add_argument("--block-size", type=int, default=4096)
    parser.add_argument("--null-reps", type=int, default=None)
    parser.add_argument("--null-seed", type=int, default=None)
    args = parser.parse_args()
    config = load_config(Path(args.config).resolve())
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stage18, stage21, stage20 = load_helpers()
    outputs = config["outputs"]
    distance_path = PROJECT_ROOT / "results/stage20/state15_cell_distances.csv"
    stage21_score_path = PROJECT_ROOT / "results/stage21/independent_lamcore_scores.csv"
    prepared_path = PROJECT_ROOT / str(outputs["prepared_h5ad"])
    scvi_path = PROJECT_ROOT / str(outputs["scvi_h5ad"])
    for path in [distance_path, stage21_score_path, prepared_path, scvi_path]:
        if not path.exists():
            raise FileNotFoundError(path)

    main_table, anchor, validation, frozen = read_frozen_inputs(distance_path, stage21_score_path, scvi_path)
    manifest = frozen["manifest"]
    latent = load_latent(stage21, scvi_path, main_table["analysis_cell_id"])
    stage21_audit_path = PROJECT_ROOT / "results/stage21/lamcore_independence_gene_audit.csv"
    if not stage21_audit_path.exists():
        raise FileNotFoundError(stage21_audit_path)
    stage21_audit = pd.read_csv(stage21_audit_path)
    scvi_hvg = set(stage21_audit.loc[stage21_audit["in_scvi_4000_hvg"].astype(bool), "gene"].astype(str))
    score_table = prepare_scores(main_table, anchor, validation, frozen["scores"], prepared_path, stage18, stage21, stage20, config, scvi_hvg, int(args.block_size))
    adjacency, directed_indices = build_local_graph(latent, main_table)
    anchor_mask = main_table["current_state"].astype(str).eq(TARGET_STATE).to_numpy()
    hops = graph_hops(adjacency, anchor_mask, max_hop=3)
    graph_table = local_graph_table(main_table, adjacency, directed_indices, hops)
    graph_table.to_csv(output_dir / "state15_local_graph_cells.csv", index=False)
    connectivity = state_connectivity(main_table.assign(state=score_table["state"].to_numpy()), adjacency, hops)
    connectivity.to_csv(output_dir / "state15_state_connectivity.csv", index=False)
    branches = select_branches(connectivity)
    branches.to_csv(output_dir / "branch_candidates.csv", index=False)

    full = score_table.copy()
    full["latent_distance_to_State15"] = pd.to_numeric(main_table["distance_to_state15"], errors="coerce").to_numpy()
    full["min_graph_hop_to_State15"] = hops
    full["fraction_of_State15_neighbors"] = graph_table["fraction_of_State15_neighbors"].to_numpy()
    full["number_of_State15_neighbors"] = graph_table["number_of_State15_neighbors"].to_numpy()
    if len(full) != len(main_table):
        raise ValueError("Graph/score row alignment failed")

    all_segment_rows: list[pd.DataFrame] = []
    all_model_rows: list[pd.DataFrame] = []
    position_tables: dict[str, pd.DataFrame] = {}
    patient_tables: list[pd.DataFrame] = []
    for _, branch in branches.iterrows():
        branch_id = str(branch["branch_id"])
        source_state = str(branch["source_state"])
        segment, models, positions = branch_gradient(branch_id, source_state, full, BRANCH_FEATURES)
        all_segment_rows.append(segment)
        all_model_rows.append(models)
        position_tables[source_state] = positions
        patient_tables.append(patient_branch_consistency(branch_id, source_state, full, BRANCH_FEATURES, min_cells=10))
    segments = pd.concat(all_segment_rows, ignore_index=True) if all_segment_rows else pd.DataFrame()
    models = pd.concat(all_model_rows, ignore_index=True) if all_model_rows else pd.DataFrame()
    positions = pd.concat(position_tables.values(), ignore_index=True) if position_tables else pd.DataFrame()
    patients = pd.concat(patient_tables, ignore_index=True) if patient_tables else pd.DataFrame()
    patient_lopo = patient_lopo_robustness(branches, full)
    segments.to_csv(output_dir / "all_branch_gradients.csv", index=False)
    models.to_csv(output_dir / "branch_gradient_models.csv", index=False)
    patients.to_csv(output_dir / "branch_patient_consistency.csv", index=False)
    patient_lopo.to_csv(output_dir / "branch_patient_lopo.csv", index=False)
    state16 = local_branch_scope(full, "State_16")
    if len(state16):
        state16_segment, state16_model, state16_positions = branch_gradient("state16", "State_16", full, BRANCH_FEATURES)
        state16_output = pd.concat([state16_segment, state16_model], ignore_index=True, sort=False)
        state16_output.to_csv(output_dir / "state16_branch_gradient.csv", index=False)
        state16_positions.to_csv(output_dir / "state16_branch_position.csv", index=False)
        state16_patients = patient_branch_consistency("state16", "State_16", full, BRANCH_FEATURES, min_cells=10)
    else:
        pd.DataFrame().to_csv(output_dir / "state16_branch_gradient.csv", index=False)
        pd.DataFrame().to_csv(output_dir / "state16_branch_position.csv", index=False)
        state16_patients = pd.DataFrame()
    state16_patients.to_csv(output_dir / "state16_patient_branch_consistency.csv", index=False)

    boundary = boundary_assignment(full, adjacency, branches)
    boundary.to_csv(output_dir / "boundary_local_branch_assignment.csv", index=False)
    extension = boundary_extension(boundary, branches, full, BRANCH_FEATURES)
    extension.to_csv(output_dir / "boundary_branch_extension.csv", index=False)

    stage22_config = config.get("stage22", {})
    null_reps = int(args.null_reps if args.null_reps is not None else stage22_config.get("matched_branch_null_reps", 500))
    null_seed = int(args.null_seed if args.null_seed is not None else stage22_config.get("matched_branch_null_seed", 20260831))
    null, null_manifest = matched_branch_null(branches, full, null_reps, null_seed)
    null.to_csv(output_dir / "branch_matched_null.csv", index=False)
    evidence = evidence_summary(branches, pd.concat([segments, models], ignore_index=True, sort=False), patients, patient_lopo, null_manifest)
    evidence.to_csv(output_dir / "branch_evidence_summary.csv", index=False)
    cp, interpretation = checkpoint(evidence)
    manifest.update(
        {
            "stage": 22,
            "branch_count": len(branches),
            "branch_null_repetitions": null_reps,
            "branch_null_seed": null_seed,
            "branch_candidates": branches.to_dict(orient="records"),
            "branch_null_manifest": null_manifest,
            "graph_hop_definition": "shortest hop on symmetrized directed kNN graph; no clustering",
            "branch_selection_patient_scope": "patient_count is the number of unique patients among state-specific hop==1 cells only",
            "branch_analysis_scope": "real branch and matched-null both restricted to non-State15 local 1–3-hop cells",
            "matched_null_distance_matching": "patient×dataset strata plus five global distance bins; nearest-distance fallback recorded per replicate",
            "matched_null_empirical_p": "direct left/right empirical tails with two-sided 2*min tail probability; no zero-centered/symmetry assumption",
            "matched_null_multiple_testing": "Benjamini-Hochberg q over all selected branches with available empirical p",
            "patient_level_robustness": "per-patient slope for branches with at least 10 local cells and LOPO patient-adjusted slopes",
            "boundary_projection_scope": "boundary cells within local 1–3-hop scope only; farther boundary cells remain unresolved/not projected",
            "score_source": "Stage21 scores for non-State15 cells; same Stage21 score function backfilled only for 200 anchor cells; Stage20 marker_VEGFD/marker_CTSK preserved as VEGFD/CTSK",
            "checkpoint": cp,
            "checkpoint_interpretation": interpretation,
        }
    )
    (output_dir / "stage22_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    write_report(output_dir, manifest, connectivity, branches, evidence, null_manifest, cp, interpretation)
    print(f"Stage 22 main cells: {len(full)}")
    print(f"Selected branch candidates: {len(branches)}")
    print(f"Checkpoint: {cp}")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
