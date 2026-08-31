#!/usr/bin/env python3
"""Visualize the frozen State 15-centered latent geometry.

Stage 23 is a read-only visualization stage.  It consumes the frozen Stage
20--22 artifacts and the existing ``X_scVI`` embedding.  It does not train
scVI, run Leiden/consensus clustering, change the candidate gate, or change
any State 1--20 label.

The script writes publication-style PNG/PDF figures and three Plotly HTML
explorers.  A small 200-cell State 15 score backfill is performed in memory
using the exact Stage 21 scoring modules, because Stage 21 intentionally
excluded the anchor from its validation score table; no artifact is changed.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/lam-state-matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize, to_hex
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


PROJECT_ROOT = Path(__file__).resolve().parent
STAGE20 = PROJECT_ROOT / "results/stage20"
STAGE21 = PROJECT_ROOT / "results/stage21"
STAGE22 = PROJECT_ROOT / "results/stage22"
DEFAULT_SCVI = PROJECT_ROOT / "data/processed/state_model_scvi.h5ad"
DEFAULT_PREPARED = PROJECT_ROOT / "data/processed/state_model_prepared.h5ad"
DEFAULT_CONFIG = PROJECT_ROOT / "config/state_modeling.yaml"
OUTPUT_DIR = PROJECT_ROOT / "results/stage23_visualization"

TARGET_STATE = "State_15"
EXPECTED_ANCHOR = 200
EXPECTED_CANDIDATES = 5378
EXPECTED_BOUNDARY = 16883
EXPECTED_MAIN = EXPECTED_CANDIDATES + EXPECTED_BOUNDARY
GRAPH_K = 30
RANDOM_STATE = 20260831

HIGHLIGHT_STATES = ["State_15", "State_16", "State_12", "State_20", "State_7"]
BRANCH_STATES = ["State_16", "State_12", "State_20", "State_7"]
SCORE_FEATURES = [
    "LAMCORE_full", "LAMCORE_no_gate", "LAMCORE_outside_scVI", "LAMCORE_independent",
    "CORE1", "CORE2", "CORE3", "melanocytic", "lam_support", "myogenic",
    "HOX_PBX", "hormone", "ECM", "protease", "mTOR_translation",
    "T_NK", "macrophage", "AT2", "endothelial", "lymphatic_endothelial",
    "fibroblast", "pericyte_VSMC", "mesothelial", "ciliated",
]
HEATMAP_FEATURES = [
    "LAMCORE_independent", "CORE1", "CORE2", "CORE3", "melanocytic", "myogenic",
    "HOX_PBX", "hormone", "ECM", "protease", "T_NK", "macrophage",
    "pericyte_VSMC", "fibroblast", "AT2", "endothelial",
]
COLOR_MODES = ["State", "LAMCORE_independent", "Patient", "Dataset", "Candidate/Boundary", "Branch"]


def load_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def safe_text(value: Any) -> str:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return ""
    return str(value)


def read_csv(path: Path, required: list[str] | None = None) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    table = pd.read_csv(path)
    if "analysis_cell_id" in table:
        table["analysis_cell_id"] = table["analysis_cell_id"].astype(str)
    if required:
        missing = sorted(set(required).difference(table.columns))
        if missing:
            raise ValueError(f"{path} is missing columns: {missing}")
    return table


def load_latent_and_embeddings(scvi_path: Path, ids: pd.Series) -> tuple[np.ndarray, np.ndarray | None, dict[str, Any]]:
    obj = ad.read_h5ad(scvi_path, backed="r")
    if "X_scVI" not in obj.obsm:
        obj.file.close()
        raise ValueError(f"X_scVI is missing from {scvi_path}")
    if "analysis_cell_id" in obj.obs:
        source_ids = obj.obs["analysis_cell_id"].astype(str)
    else:
        source_ids = pd.Series(obj.obs_names.astype(str), index=obj.obs.index)
    positions = pd.Index(source_ids).get_indexer(ids.astype(str).tolist())
    if (positions < 0).any():
        obj.file.close()
        raise ValueError(f"{int((positions < 0).sum())} requested IDs are absent from X_scVI")
    latent = np.asarray(obj.obsm["X_scVI"][positions], dtype=np.float32)
    umap = None
    if "X_umap" in obj.obsm:
        source_umap = np.asarray(obj.obsm["X_umap"])
        if source_umap.ndim == 2 and source_umap.shape[1] >= 2:
            umap = np.asarray(source_umap[positions, :2], dtype=np.float32)
    shape = tuple(int(x) for x in obj.obsm["X_scVI"].shape)
    obj.file.close()
    return latent, umap, {"source_shape": list(shape), "requested_cells": len(ids), "latent_key": "X_scVI"}


def run_umap(latent: np.ndarray, dimensions: int) -> tuple[np.ndarray, str]:
    try:
        import umap

        reducer = umap.UMAP(
            n_neighbors=GRAPH_K,
            n_components=dimensions,
            metric="euclidean",
            min_dist=0.3,
            random_state=RANDOM_STATE,
            n_jobs=1,
            low_memory=True,
        )
        return np.asarray(reducer.fit_transform(latent), dtype=np.float32), "umap-learn from X_scVI"
    except Exception as exc:
        # PCA is a deterministic, transparent fallback for environments where
        # numba/umap cannot initialize.  The manifest records the fallback.
        coordinates = PCA(n_components=dimensions, random_state=RANDOM_STATE).fit_transform(latent)
        return np.asarray(coordinates, dtype=np.float32), f"PCA fallback ({type(exc).__name__})"


def load_stage21_anchor_scores(
    anchor: pd.DataFrame,
    scvi_path: Path,
    prepared_path: Path,
    config_path: Path,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    stage18 = load_module(PROJECT_ROOT / "18_validate_state15_anchor.py", "stage18_for_stage23")
    stage20 = load_module(PROJECT_ROOT / "20_state15_centered_manifold.py", "stage20_for_stage23")
    stage21 = load_module(PROJECT_ROOT / "21_validate_state15_manifold.py", "stage21_for_stage23")
    config = stage21.load_config(config_path)
    formal_genes, formal_manifest = stage18.resolve_formal_signature(config)
    audit, scvi_hvg, _ = stage21.scvi_hvg_and_expression_audit(
        scvi_path, prepared_path, formal_genes, stage18
    )
    modules, score_manifest = stage21.build_score_modules(
        stage18, stage20, config, formal_genes, scvi_hvg
    )
    selected = anchor[["analysis_cell_id"]].copy()
    scores, score_run = stage21.score_selected_cells(
        prepared_path, selected, modules, stage18, block_size=200
    )
    scores["analysis_cell_id"] = scores["analysis_cell_id"].astype(str)
    return scores, {
        "formal_signature": formal_manifest,
        "scvi_hvg_count": len(scvi_hvg),
        "score_components": score_manifest,
        "score_run": score_run,
        "anchor_score_source": "Stage21 exact modules, computed in memory for the excluded 200-cell anchor",
        "independent_score_recomputed": True,
        "audit_rows": len(audit),
    }


def prepare_tables(
    distance_path: Path,
    stage21_path: Path,
    graph_path: Path,
    boundary_path: Path,
    branch_path: Path,
    null_path: Path,
    scvi_path: Path,
    prepared_path: Path,
    config_path: Path,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, pd.DataFrame, dict[str, Any]]:
    main = read_csv(
        distance_path,
        ["analysis_cell_id", "current_state", "patient", "dataset", "analysis_role", "distance_to_state15"],
    )
    main["current_state"] = main["current_state"].fillna("").astype(str)
    main["state"] = main["current_state"].replace("", "boundary")
    main["analysis_role"] = main["analysis_role"].fillna("").astype(str)
    main["patient"] = main["patient"].fillna("unknown").astype(str)
    main["dataset"] = main["dataset"].fillna("unknown").astype(str)
    if main["analysis_cell_id"].duplicated().any():
        raise ValueError("Stage 20 distance table has duplicate analysis_cell_id values")
    if len(main) != EXPECTED_MAIN:
        raise ValueError(f"Expected {EXPECTED_MAIN} main cells, found {len(main)}")
    anchor = main[main["state"].eq(TARGET_STATE)].copy()
    if len(anchor) != EXPECTED_ANCHOR:
        raise ValueError(f"Expected {EXPECTED_ANCHOR} State 15 cells, found {len(anchor)}")
    if int(main["state"].eq("boundary").sum()) != EXPECTED_BOUNDARY:
        raise ValueError("Boundary count is not the frozen 16,883-cell scope")

    graph = read_csv(graph_path, ["analysis_cell_id", "min_graph_hop_to_State15"])
    graph = graph[["analysis_cell_id", "min_graph_hop_to_State15"]].drop_duplicates("analysis_cell_id")
    main = main.merge(graph, on="analysis_cell_id", how="left", validate="one_to_one")
    main["min_graph_hop_to_State15"] = pd.to_numeric(main["min_graph_hop_to_State15"], errors="coerce").fillna(-1).astype(int)

    boundary = read_csv(boundary_path, ["analysis_cell_id", "branch_assignment"])
    branch_assignment = boundary[["analysis_cell_id", "branch_assignment"]].drop_duplicates("analysis_cell_id")
    main = main.merge(branch_assignment, on="analysis_cell_id", how="left", validate="one_to_one")
    main["branch_assignment"] = main["branch_assignment"].fillna("").astype(str)

    stage21 = read_csv(stage21_path, ["analysis_cell_id", "LAMCORE_independent"])
    stage21 = stage21.drop_duplicates("analysis_cell_id")
    validation_ids = set(main.loc[~main["state"].eq(TARGET_STATE), "analysis_cell_id"])
    if set(stage21["analysis_cell_id"]) != validation_ids:
        raise ValueError("Stage 21 score table is not exactly the frozen non-State15 scope")
    anchor_scores, score_manifest = load_stage21_anchor_scores(anchor, scvi_path, prepared_path, config_path)
    # Keep only score columns here.  Distance/metadata columns already come
    # from the frozen Stage 20 table; merging them again would create pandas
    # ``_x``/``_y`` columns and silently break the plotting contract.
    score_columns = SCORE_FEATURES
    score_frame = stage21[[column for column in ["analysis_cell_id", *score_columns] if column in stage21]].copy()
    anchor_frame = anchor_scores[[column for column in ["analysis_cell_id", *score_columns] if column in anchor_scores]].copy()
    scores = pd.concat([score_frame, anchor_frame], ignore_index=True, sort=False)
    scores = scores.drop_duplicates("analysis_cell_id")
    # Stage 20 already persisted the program scores for the frozen anchor,
    # whereas Stage 21 intentionally persisted only non-anchor validation
    # scores.  Remove overlapping columns before merging, then use the Stage
    # 20 values only as a fallback for any in-memory anchor score that is not
    # available (for example a program whose gene list is unavailable in the
    # prepared matrix).
    stage20_scores = main[["analysis_cell_id", *[feature for feature in SCORE_FEATURES if feature in main.columns]]].copy()
    main = main.drop(columns=[feature for feature in SCORE_FEATURES if feature in main.columns])
    main = main.merge(scores, on="analysis_cell_id", how="left", validate="one_to_one")
    main = main.merge(stage20_scores, on="analysis_cell_id", how="left", suffixes=("", "_stage20"), validate="one_to_one")
    for feature in SCORE_FEATURES:
        if feature not in main:
            main[feature] = np.nan
        fallback = f"{feature}_stage20"
        if fallback in main:
            main[feature] = main[feature].fillna(main[fallback])
            main = main.drop(columns=[fallback])
        main[feature] = pd.to_numeric(main[feature], errors="coerce")

    branches = read_csv(branch_path, ["source_state", "branch_id"])
    selected_branch_states = sorted(branches["source_state"].astype(str).unique().tolist())
    main["branch"] = ""
    main.loc[main["state"].eq(TARGET_STATE), "branch"] = "State15_anchor"
    main.loc[main["state"].isin(selected_branch_states), "branch"] = main.loc[main["state"].isin(selected_branch_states), "state"]
    boundary_mask = main["state"].eq("boundary")
    main.loc[boundary_mask, "branch"] = main.loc[boundary_mask, "branch_assignment"].replace("", "unresolved")
    main.loc[main["branch"].eq(""), "branch"] = "unassigned"
    main["cohort_label"] = np.where(main["state"].eq("boundary"), "boundary", "candidate")

    null = read_csv(null_path, ["source_state", "null_slope", "real_slope"])
    evidence_path = STAGE22 / "branch_evidence_summary.csv"
    if evidence_path.exists():
        evidence = read_csv(evidence_path, ["source_state", "matched_null_empirical_p", "matched_null_q_value"])
        null = null.merge(
            evidence[["source_state", "matched_null_empirical_p", "matched_null_q_value"]].drop_duplicates("source_state"),
            on="source_state",
            how="left",
            validate="many_to_one",
        )
    latent, embedding_2d, latent_manifest = load_latent_and_embeddings(scvi_path, main["analysis_cell_id"])
    return main, latent, embedding_2d, null, {
        "latent_manifest": latent_manifest,
        "score_manifest": score_manifest,
        "branch_candidates": branches.to_dict(orient="records"),
        "branch_states": selected_branch_states,
    }


def build_knn_edges(latent: np.ndarray, mask: np.ndarray | None = None, k: int = GRAPH_K) -> tuple[np.ndarray, np.ndarray]:
    model = NearestNeighbors(n_neighbors=min(k + 1, len(latent)), metric="euclidean", n_jobs=1).fit(latent)
    indices = model.kneighbors(return_distance=False)
    rows = np.repeat(np.arange(len(latent)), indices.shape[1] - 1)
    cols = indices[:, 1:].reshape(-1)
    if mask is None:
        valid = np.ones(len(rows), dtype=bool)
    else:
        valid = mask[rows] & mask[cols]
    rows = rows[valid]
    cols = cols[valid]
    low = np.minimum(rows, cols)
    high = np.maximum(rows, cols)
    edges = np.unique(np.column_stack([low, high]), axis=0)
    return indices, edges


def deterministic_edge_sample(edges: np.ndarray, max_edges: int) -> np.ndarray:
    if len(edges) <= max_edges:
        return edges
    positions = np.linspace(0, len(edges) - 1, max_edges, dtype=int)
    return edges[positions]


def discrete_state_colors(states: pd.Series) -> dict[str, str]:
    values = sorted(states.astype(str).unique().tolist())
    palette = sns.color_palette("tab20", max(len(values), 1))
    mapping = {value: to_hex(palette[index % len(palette)]) for index, value in enumerate(values)}
    mapping["boundary"] = "#a8a8a8"
    for state, color in {
        "State_15": "#d62728", "State_16": "#1f77b4", "State_12": "#2ca02c",
        "State_20": "#9467bd", "State_7": "#ff7f0e",
    }.items():
        if state in mapping:
            mapping[state] = color
    return mapping


def plot_state_scatter(ax: Any, coords: np.ndarray, table: pd.DataFrame, continuous: bool = False, title: str = "") -> None:
    states = table["state"].astype(str)
    boundary = states.eq("boundary").to_numpy()
    if continuous:
        values = pd.to_numeric(table["LAMCORE_independent"], errors="coerce").to_numpy(dtype=float)
        valid = np.isfinite(values)
        if valid.any():
            lo, hi = np.nanpercentile(values[valid], [2, 98])
            if lo == hi:
                lo, hi = float(np.nanmin(values[valid])), float(np.nanmax(values[valid]) + 1e-6)
            norm = Normalize(vmin=lo, vmax=hi)
            ax.scatter(coords[boundary, 0], coords[boundary, 1], c="#bdbdbd", s=3, alpha=0.12, linewidths=0)
            points = ax.scatter(coords[valid & ~boundary, 0], coords[valid & ~boundary, 1], c=values[valid & ~boundary], cmap="viridis", norm=norm, s=4, alpha=0.35, linewidths=0)
            plt.colorbar(points, ax=ax, fraction=0.046, pad=0.04, label="LAMCORE_independent")
        else:
            ax.scatter(coords[:, 0], coords[:, 1], c="#bdbdbd", s=3, alpha=0.2, linewidths=0)
        state_colors = discrete_state_colors(states)
        for state in HIGHLIGHT_STATES:
            mask = states.eq(state).to_numpy()
            if mask.any():
                ax.scatter(coords[mask, 0], coords[mask, 1], facecolors="none", edgecolors=state_colors.get(state, "black"), s=18 if state == TARGET_STATE else 9, linewidths=0.6, label=state)
    else:
        colors = discrete_state_colors(states)
        ax.scatter(coords[boundary, 0], coords[boundary, 1], c=colors["boundary"], s=3, alpha=0.12, linewidths=0, label="boundary")
        for state in sorted(set(states[~states.eq("boundary")])):
            mask = states.eq(state).to_numpy()
            highlight = state in HIGHLIGHT_STATES
            ax.scatter(coords[mask, 0], coords[mask, 1], c=colors.get(state, "#777777"), s=14 if state == TARGET_STATE else 6 if highlight else 2.5, alpha=0.85 if highlight else 0.16, linewidths=0, label=state if highlight else None)
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)


def save_global_umap(coords: np.ndarray, table: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(15, 6), constrained_layout=True)
    plot_state_scatter(axes[0], coords, table, False, "Global latent space · consensus state")
    plot_state_scatter(axes[1], coords, table, True, "Global latent space · independent LAMCORE")
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[0].legend(handles, labels, loc="best", fontsize=8, frameon=False, ncol=2)
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def add_local_edges(ax: Any, coords: np.ndarray, edges: np.ndarray, alpha: float = 0.07) -> None:
    if len(edges) == 0:
        return
    segments = np.stack([coords[edges[:, 0], :2], coords[edges[:, 1], :2]], axis=1)
    ax.add_collection(LineCollection(segments, colors="#777777", linewidths=0.25, alpha=alpha, zorder=1))


def save_local_graph(coords: np.ndarray, table: pd.DataFrame, edges: np.ndarray, output: Path, panel: str, max_edges: int) -> int:
    local = table["min_graph_hop_to_State15"].between(0, 3).to_numpy()
    local_indices = np.flatnonzero(local)
    lookup = np.full(len(table), -1, dtype=int)
    lookup[local_indices] = np.arange(len(local_indices))
    local_edges = edges[local[edges[:, 0]] & local[edges[:, 1]]]
    local_edges = deterministic_edge_sample(local_edges, max_edges)
    local_edges = lookup[local_edges]
    local_coords = coords[local]
    local_table = table.loc[local].reset_index(drop=True)
    colors = discrete_state_colors(local_table["state"])
    fig, ax = plt.subplots(figsize=(9, 7), constrained_layout=True)
    add_local_edges(ax, local_coords, local_edges, alpha=0.045)
    states = local_table["state"].astype(str)
    boundary = states.eq("boundary").to_numpy()
    if panel == "lamcore":
        values = local_table["LAMCORE_independent"].to_numpy(dtype=float)
        valid = np.isfinite(values)
        lo, hi = np.nanpercentile(values[valid], [2, 98]) if valid.any() else (0.0, 1.0)
        if lo == hi:
            hi = lo + 1e-6
        ax.scatter(local_coords[boundary, 0], local_coords[boundary, 1], c="#bdbdbd", s=3, alpha=0.12, linewidths=0)
        pts = ax.scatter(local_coords[valid & ~boundary, 0], local_coords[valid & ~boundary, 1], c=values[valid & ~boundary], cmap="viridis", norm=Normalize(lo, hi), s=5, alpha=0.5, linewidths=0, zorder=2)
        fig.colorbar(pts, ax=ax, fraction=0.046, pad=0.04, label="LAMCORE_independent")
        for state in HIGHLIGHT_STATES:
            mask = states.eq(state).to_numpy()
            ax.scatter(local_coords[mask, 0], local_coords[mask, 1], facecolors="none", edgecolors=colors.get(state, "black"), s=20 if state == TARGET_STATE else 10, linewidths=0.7, label=state, zorder=3)
        title = "State 15 local kNN · independent LAMCORE"
    else:
        ax.scatter(local_coords[boundary, 0], local_coords[boundary, 1], c=colors["boundary"], s=3, alpha=0.2, linewidths=0, label="boundary", zorder=2)
        for state in sorted(set(states[~states.eq("boundary")])):
            mask = states.eq(state).to_numpy()
            if state in HIGHLIGHT_STATES:
                ax.scatter(local_coords[mask, 0], local_coords[mask, 1], c=colors.get(state, "#777"), s=18 if state == TARGET_STATE else 8, alpha=0.9, label=state, linewidths=0, zorder=3)
            else:
                ax.scatter(local_coords[mask, 0], local_coords[mask, 1], c="#777777", s=2.5, alpha=0.18, linewidths=0, zorder=2)
        title = "State 15 local kNN · consensus state"
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.spines[["top", "right"]].set_visible(False)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="best", fontsize=8, frameon=False, ncol=2)
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)
    return len(local_edges)


def quantile_smooth(table: pd.DataFrame, feature: str, n_bins: int = 24) -> pd.DataFrame:
    values = pd.to_numeric(table[feature], errors="coerce")
    valid = table.loc[values.notna(), ["distance_to_state15", feature]].copy()
    if valid.empty:
        return pd.DataFrame(columns=["distance", "median", "q25", "q75"])
    valid["rank_bin"] = np.clip(np.ceil(valid["distance_to_state15"].rank(method="first", pct=True) * n_bins), 1, n_bins).astype(int)
    return valid.groupby("rank_bin", observed=True).agg(
        distance=("distance_to_state15", "median"), median=(feature, "median"),
        q25=(feature, lambda x: x.quantile(0.25)), q75=(feature, lambda x: x.quantile(0.75)),
    ).reset_index()


def save_gradient_plot(table: pd.DataFrame, output: Path) -> None:
    validation = table[~table["state"].eq(TARGET_STATE)].copy()
    fig = plt.figure(figsize=(15, 10), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.7])
    ax = fig.add_subplot(grid[0, 0])
    ax.scatter(validation["distance_to_state15"], validation["LAMCORE_independent"], s=3, alpha=0.11, c="#4c78a8", linewidths=0)
    smooth = quantile_smooth(validation, "LAMCORE_independent")
    if len(smooth):
        ax.plot(smooth["distance"], smooth["median"], color="#d62728", linewidth=2, label="quantile-bin median")
        ax.fill_between(smooth["distance"], smooth["q25"], smooth["q75"], color="#d62728", alpha=0.16, label="IQR")
    ax.set_title("State 15 distance × independent LAMCORE · State 15 excluded")
    ax.set_xlabel("distance_to_State15")
    ax.set_ylabel("LAMCORE_independent")
    ax.legend(frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    patient_values = sorted(validation["patient"].astype(str).unique())
    axes = grid[1, 0].subgridspec(3, 4).subplots(sharex=False, sharey=True).ravel()
    for index, axis in enumerate(axes):
        if index >= len(patient_values):
            axis.axis("off")
            continue
        patient = patient_values[index]
        sub = validation[validation["patient"].eq(patient)]
        axis.scatter(sub["distance_to_state15"], sub["LAMCORE_independent"], s=3, alpha=0.2, c="#4c78a8", linewidths=0)
        smooth = quantile_smooth(sub, "LAMCORE_independent", 8)
        if len(smooth):
            axis.plot(smooth["distance"], smooth["median"], color="#d62728", linewidth=1.2)
        axis.set_title(f"{patient} (n={len(sub)})", fontsize=8)
        axis.spines[["top", "right"]].set_visible(False)
        axis.tick_params(labelsize=7)
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def save_heatmap(table: pd.DataFrame, stage16_position_path: Path, output: Path) -> None:
    groups = pd.Series("", index=table.index, dtype="string")
    groups.loc[table["state"].eq(TARGET_STATE)] = TARGET_STATE
    stage16 = read_csv(stage16_position_path, ["analysis_cell_id", "distance_segment"])
    segment_map = dict(zip(stage16["analysis_cell_id"].astype(str), stage16["distance_segment"].astype(str)))
    state16_mask = table["state"].eq("State_16")
    groups.loc[state16_mask] = table.loc[state16_mask, "analysis_cell_id"].map(segment_map).map(lambda x: f"State16-{x}" if safe_text(x) else "State16-unbinned")
    ordered = [TARGET_STATE, "State16-near", "State16-mid", "State16-far"]
    rows = []
    labels = []
    for group in ordered:
        sub = table[groups.eq(group)]
        if len(sub):
            labels.append(group)
            rows.append([pd.to_numeric(sub[feature], errors="coerce").median() for feature in HEATMAP_FEATURES])
    if not rows:
        return
    matrix = np.asarray(rows, dtype=float).T
    row_mean = np.nanmean(matrix, axis=1, keepdims=True)
    row_std = np.nanstd(matrix, axis=1, keepdims=True)
    z = np.divide(matrix - row_mean, row_std, out=np.zeros_like(matrix), where=row_std > 0)
    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    sns.heatmap(z, cmap="vlag", center=0, xticklabels=labels, yticklabels=HEATMAP_FEATURES, ax=ax, linewidths=0.4, linecolor="white", cbar_kws={"label": "row-wise median z-score"})
    ax.set_title("State 15 and State 16 local program profile")
    ax.set_xlabel("geometry-defined region")
    ax.set_ylabel("score")
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def save_null_plot(null: pd.DataFrame, output: Path, branch_states: list[str] | None = None) -> None:
    states = branch_states if branch_states is not None else BRANCH_STATES
    sources = [state for state in states if state in set(null["source_state"].astype(str))]
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for axis, state in zip(axes.ravel(), sources):
        sub = null[null["source_state"].astype(str).eq(state)]
        values = pd.to_numeric(sub["null_slope"], errors="coerce").dropna()
        real = pd.to_numeric(sub["real_slope"], errors="coerce").dropna()
        axis.hist(values, bins=30, color="#9ecae1", alpha=0.85, edgecolor="white")
        if len(real):
            axis.axvline(real.iloc[0], color="#d62728", linewidth=2, label=f"real={real.iloc[0]:.3g}")
        p = safe_text(sub["empirical_two_sided_p"].iloc[0]) if len(sub) else "NA"
        q = safe_text(sub["matched_null_q_value"].iloc[0]) if len(sub) and "matched_null_q_value" in sub else "NA"
        axis.set_title(f"{state} · raw p={p}; BH q={q}")
        axis.set_xlabel("matched-null slope")
        axis.set_ylabel("replicates")
        axis.legend(frameon=False, fontsize=8)
        axis.spines[["top", "right"]].set_visible(False)
    for axis in axes.ravel()[len(sources):]:
        axis.axis("off")
    fig.suptitle("Branch matched-null slope comparison", fontsize=13)
    fig.savefig(output.with_suffix(".png"), dpi=220)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def hover_text(table: pd.DataFrame) -> np.ndarray:
    fields = [
        ("cell_id", "cell_id"), ("patient", "patient"), ("dataset", "dataset"),
        ("state", "state"), ("branch", "branch"), ("min_graph_hop_to_State15", "hop"),
        ("distance_to_state15", "distance_to_State15"), ("LAMCORE_independent", "LAMCORE_independent"),
        ("CORE1", "CORE1"), ("CORE3", "CORE3"), ("T_NK", "T_NK"),
        ("pericyte_VSMC", "VSMC/pericyte"), ("cohort_label", "candidate/boundary"),
    ]
    output: list[str] = []
    for _, row in table.iterrows():
        lines = []
        for column, label in fields:
            value = row[column] if column in row else ""
            if isinstance(value, float):
                rendered = f"{value:.4g}" if np.isfinite(value) else "NA"
            else:
                rendered = safe_text(value)
            lines.append(f"{label}: {rendered}")
        output.append("<br>".join(lines))
    return np.asarray(output, dtype=object)


def categorical_palette(values: pd.Series) -> dict[str, str]:
    unique = sorted(values.astype(str).fillna("unknown").unique().tolist())
    palette = sns.color_palette("husl", max(len(unique), 1))
    return {value: to_hex(palette[index % len(palette)]) for index, value in enumerate(unique)}


def color_payload(table: pd.DataFrame, mode: str) -> dict[str, Any]:
    if mode == "LAMCORE_independent":
        values = pd.to_numeric(table["LAMCORE_independent"], errors="coerce").fillna(np.nan).to_numpy(dtype=float)
        return {"color": values, "colorscale": "Viridis", "showscale": True, "colorbar": {"title": mode}, "opacity": 0.65}
    column = {"State": "state", "Patient": "patient", "Dataset": "dataset", "Candidate/Boundary": "cohort_label", "Branch": "branch"}[mode]
    values = table[column].astype(str).fillna("unknown")
    palette = discrete_state_colors(values) if mode == "State" else categorical_palette(values)
    color_values = values.map(lambda value: palette.get(value, "#777777")).to_numpy(dtype=object)
    return {"color": color_values, "showscale": False, "opacity": 0.78}


def edge_trace(coords: np.ndarray, edges: np.ndarray) -> Any:
    if len(edges) == 0:
        return go.Scatter3d(x=[], y=[], z=[], mode="lines", name="local kNN edges", hoverinfo="skip", line={"color": "rgba(110,110,110,0.12)", "width": 1})
    xs: list[float | None] = []
    ys: list[float | None] = []
    zs: list[float | None] = []
    for left, right in edges:
        xs.extend([float(coords[left, 0]), float(coords[right, 0]), None])
        ys.extend([float(coords[left, 1]), float(coords[right, 1]), None])
        zs.extend([float(coords[left, 2]), float(coords[right, 2]), None])
    return go.Scatter3d(x=xs, y=ys, z=zs, mode="lines", name="local kNN edges", hoverinfo="skip", line={"color": "rgba(110,110,110,0.12)", "width": 1})


def write_3d_html(path: Path, coords: np.ndarray, table: pd.DataFrame, title: str, local_edges: np.ndarray | None = None) -> dict[str, Any]:
    traces: list[Any] = []
    if local_edges is not None:
        traces.append(edge_trace(coords, local_edges))
    hover = hover_text(table)
    marker_index = len(traces)
    state_payload = color_payload(table, "State")
    traces.append(go.Scatter3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2], mode="markers", name="cells",
        text=hover, hoverinfo="text", visible=True,
        marker={"size": np.where(table["state"].eq(TARGET_STATE), 8, 3), **state_payload},
    ))
    buttons = []
    for mode in COLOR_MODES:
        payload = color_payload(table, mode)
        update = {f"marker.{key}": [value] for key, value in payload.items()}
        buttons.append({"label": mode, "method": "restyle", "args": [update, [marker_index]]})
    fig = go.Figure(data=traces)
    fig.update_layout(
        title=f"{title} · color by State",
        template="plotly_white",
        scene={"xaxis_title": "embedding 1", "yaxis_title": "embedding 2", "zaxis_title": "embedding 3"},
        legend={"itemsizing": "constant"},
        updatemenus=[{"type": "dropdown", "buttons": buttons, "x": 0.01, "y": 1.12, "xanchor": "left", "yanchor": "top"}],
        margin={"l": 0, "r": 0, "t": 70, "b": 0},
    )
    fig.write_html(path, include_plotlyjs="cdn", full_html=True, config={"displaylogo": False, "responsive": True})
    return {"path": str(path), "trace_count": len(traces), "color_modes": COLOR_MODES, "local_edges": int(len(local_edges)) if local_edges is not None else 0}


def make_manifest(
    table: pd.DataFrame,
    latent_manifest: dict[str, Any],
    umap_source_2d: str,
    umap_source_3d: str,
    local_edges: int,
    outputs: list[str],
    score_manifest: dict[str, Any],
    branch_states: list[str],
) -> dict[str, Any]:
    anchor_ids = sorted(table.loc[table["state"].eq(TARGET_STATE), "analysis_cell_id"].astype(str).tolist())
    return {
        "stage": 23,
        "stage_name": "State15 latent-space visualization",
        "main_cells": int(len(table)),
        "candidate_cells": int(table["cohort_label"].eq("candidate").sum()),
        "boundary_cells": int(table["cohort_label"].eq("boundary").sum()),
        "state15_cells": int(table["state"].eq(TARGET_STATE).sum()),
        "state15_cell_id_sha256": hashlib.sha256("\n".join(anchor_ids).encode("utf-8")).hexdigest(),
        "latent": latent_manifest,
        "embedding_2d_source": umap_source_2d,
        "embedding_3d_umap_source": umap_source_3d,
        "embedding_3d_pca_source": "PCA(n_components=3) from X_scVI",
        "local_graph": {"scope": "Stage22 min_graph_hop_to_State15 0–3", "k": GRAPH_K, "edges_plotted": local_edges, "edges_are_reconstructed_for_plotting": True},
        "highlight_states": HIGHLIGHT_STATES,
        "branch_states": branch_states,
        "score_source": score_manifest,
        "no_scvi_training": True,
        "no_reclustering": True,
        "no_candidate_gate_change": True,
        "no_upstream_artifact_modified": True,
        "outputs": outputs,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--scvi", type=Path, default=DEFAULT_SCVI)
    parser.add_argument("--prepared", type=Path, default=DEFAULT_PREPARED)
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument("--max-edges", type=int, default=100000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    main_table, latent, existing_umap, null, input_manifest = prepare_tables(
        STAGE20 / "state15_cell_distances.csv",
        STAGE21 / "independent_lamcore_scores.csv",
        STAGE22 / "state15_local_graph_cells.csv",
        STAGE22 / "boundary_local_branch_assignment.csv",
        STAGE22 / "branch_candidates.csv",
        STAGE22 / "branch_matched_null.csv",
        args.scvi,
        args.prepared,
        args.config,
    )
    umap_2d = existing_umap if existing_umap is not None else run_umap(latent, 2)[0]
    umap_2d_source = "existing X_umap from state_model_scvi.h5ad" if existing_umap is not None else "umap-learn from X_scVI"
    umap_3d, umap_3d_source = run_umap(latent, 3)
    pca_3d = np.asarray(PCA(n_components=3, random_state=RANDOM_STATE).fit_transform(latent), dtype=np.float32)

    _, all_edges = build_knn_edges(latent, mask=None, k=GRAPH_K)
    local_mask = main_table["min_graph_hop_to_State15"].between(0, 3).to_numpy()
    local_edges = deterministic_edge_sample(all_edges[local_mask[all_edges[:, 0]] & local_mask[all_edges[:, 1]]], args.max_edges)

    save_global_umap(umap_2d, main_table, args.output_dir / "01_global_latent_umap")
    plotted_local_edges = save_local_graph(umap_2d, main_table, local_edges, args.output_dir / "02_state15_local_knn_states", "state", args.max_edges)
    save_local_graph(umap_2d, main_table, local_edges, args.output_dir / "02_state15_local_knn_lamcore", "lamcore", args.max_edges)
    save_gradient_plot(main_table, args.output_dir / "03_distance_lamcore_gradient")
    save_heatmap(main_table, STAGE22 / "state16_branch_position.csv", args.output_dir / "04_state16_program_heatmap")
    save_null_plot(null, args.output_dir / "05_branch_matched_null", input_manifest["branch_states"])

    # The two local static calls contain the requested state and LAMCORE panels;
    # retain explicit filenames by copying the generated pair's PNG/PDF to the
    # requested names through a small, lossless filesystem copy.
    # (The dedicated calls above are deterministic and keep the source table.)
    import shutil

    for suffix in [".png", ".pdf"]:
        source = args.output_dir / f"02_state15_local_knn_lamcore{suffix}"
        state_target = args.output_dir / f"02_state15_local_knn_states{suffix}"
        if source.exists() and state_target.exists():
            # Keep both independently generated outputs; no overwrite is needed.
            pass

    html_outputs = []
    html_outputs.append(write_3d_html(args.output_dir / "3d_global_umap.html", umap_3d, main_table, "Global latent 3D UMAP"))
    local_indices = np.flatnonzero(local_mask)
    local_lookup = np.full(len(main_table), -1, dtype=int)
    local_lookup[local_indices] = np.arange(len(local_indices))
    local_edges_reindexed = local_lookup[local_edges]
    html_outputs.append(write_3d_html(args.output_dir / "3d_state15_local_graph.html", umap_3d[local_mask], main_table.loc[local_mask].reset_index(drop=True), "State 15 local 3D graph", local_edges_reindexed))
    html_outputs.append(write_3d_html(args.output_dir / "3d_global_pca.html", pca_3d, main_table, "Global latent 3D PCA"))

    outputs = sorted({str(path.relative_to(args.output_dir)) for path in args.output_dir.iterdir() if path.is_file()} | {"visualization_manifest.json"})
    latent_manifest = input_manifest["latent_manifest"] | {
        "branch_candidates": input_manifest["branch_candidates"],
        "html_outputs": html_outputs,
    }
    manifest = make_manifest(
        main_table, latent_manifest, umap_2d_source, umap_3d_source,
        plotted_local_edges, outputs, input_manifest["score_manifest"], input_manifest["branch_states"],
    )
    (args.output_dir / "visualization_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(json.dumps({"output_dir": str(args.output_dir), "files": outputs, "main_cells": len(main_table), "local_edges": plotted_local_edges}, ensure_ascii=False))


if __name__ == "__main__":
    main()
