#!/usr/bin/env python3
"""Map the existing transcriptomic neighborhood around frozen State 15.

Stage 20 is geometry-first and read-only with respect to the existing model:
it uses the current ``X_scVI`` embedding, the current State 1--20 labels, and
the inherited boundary labels.  It does not train scVI, run Leiden, rebuild a
consensus, or change a candidate gate.
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
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.neighbors import NearestNeighbors


PROJECT_ROOT = Path(__file__).resolve().parent
TARGET_STATE = "15"
TARGET_PATIENT = "LAM1163"
EXPECTED_STATE15_CELLS = 200
EXPECTED_CANDIDATE_CELLS = 5378
EXPECTED_BOUNDARY_CELLS = 16883
DISTANCE_BIN_LABELS = ["0-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-100%"]
DISTANCE_BIN_EDGES = np.asarray([0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])
COMPARISON_STATES = ["18", "20", "12", "7", "5"]


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def load_stage18_module() -> Any:
    path = PROJECT_ROOT / "18_validate_state15_anchor.py"
    spec = importlib.util.spec_from_file_location("stage18_anchor", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Stage 18 helpers: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_modules(stage18: Any, config: dict[str, Any]) -> tuple[dict[str, list[str]], dict[str, Any], dict[str, Any]]:
    programs, program_manifest = stage18.load_programs(config)
    formal_genes, formal_manifest = stage18.resolve_formal_signature(config)
    stage16 = stage18.load_stage16_module()
    modules: dict[str, list[str]] = {
        "melanocytic": ["PMEL", "MLANA", "MITF"],
        "lam_support": ["VEGFD", "CTSK", "ESR1"],
        "myogenic": ["ACTA2", "ACTG2", "MYH11"],
        "HOX_PBX_markers": ["EMX2", "HOXA11"],
        "T_NK": [
            "CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2", "IL7R", "LTB",
            "NKG7", "GNLY", "GZMB", "PRF1", "KLRD1", "CD8A", "CD8B",
        ],
        "endothelial": ["PECAM1", "EMCN", "VWF", "KDR", "ESAM", "ENG"],
        "lymphatic_endothelial": ["CCL21", "FLT4", "PDPN", "LYVE1", "PROX1"],
    }
    for name in [
        "CORE1", "CORE2", "CORE3_identity", "LAM_myogenic_contractile",
        "ECM_remodeling", "mTOR_translation", "hormone_related",
        "protease_ECM_niche", "HOX_PBX",
    ]:
        if name in programs:
            output_name = {
                "CORE3_identity": "CORE3",
                "LAM_myogenic_contractile": "LAM_myogenic",
                "hormone_related": "hormone",
                "ECM_remodeling": "ECM",
                "protease_ECM_niche": "protease",
            }.get(name, name)
            modules[output_name] = programs[name]
    if formal_genes:
        modules["LAMCORE_777"] = formal_genes
    modules.update({name: list(genes) for name, genes in stage16.COMPETING_GENES.items()})
    for gene in [
        "PMEL", "MLANA", "MITF", "VEGFD", "CTSK", "ESR1", "ACTA2", "ACTG2",
        "MYH11", "NKG7", "GNLY", "GZMB", "CD3D", "CD3E",
    ]:
        modules[f"marker_{gene}"] = [gene]
    modules = {name: stage18.unique_genes(genes) for name, genes in modules.items()}
    return modules, formal_manifest, program_manifest


def bool_series(stage18: Any, values: pd.Series) -> pd.Series:
    return values.map(stage18.as_bool).astype(bool)


def prepared_metadata(stage18: Any, prepared: ad.AnnData) -> pd.DataFrame:
    return stage18.prepared_obs(prepared)


def fixed_state15_manifest(
    consensus: pd.DataFrame,
    scvi_path: Path,
    latent_shape: tuple[int, int],
) -> dict[str, Any]:
    state15_ids = sorted(
        consensus.loc[consensus["consensus_state"].astype(str).eq(TARGET_STATE), "analysis_cell_id"].astype(str).tolist()
    )
    return {
        "state": TARGET_STATE,
        "state15_cell_count": len(state15_ids),
        "state15_cell_id_sha256": hashlib.sha256("\n".join(state15_ids).encode("utf-8")).hexdigest(),
        "state15_cell_ids": state15_ids,
        "latent_artifact": str(scvi_path),
        "latent_key": "X_scVI",
        "latent_shape": list(latent_shape),
        "neighbor_graph_scope": "5,378 high-confidence candidates + 16,883 inherited boundary cells",
        "normal_reference_scope": "remote summary only; excluded from State 15 neighborhood graph",
        "no_scvi_training": True,
        "no_reclustering": True,
        "no_candidate_gate_change": True,
    }


def load_scvi_cohort(
    scvi_path: Path,
    consensus: pd.DataFrame,
    stage18: Any,
) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame, np.ndarray, dict[str, Any]]:
    obj = ad.read_h5ad(scvi_path, backed="r")
    if "X_scVI" not in obj.obsm:
        obj.file.close()
        raise ValueError(f"X_scVI not found in {scvi_path}")
    obs = obj.obs.copy().reset_index(names="obs_name")
    obs["analysis_cell_id"] = obs["analysis_cell_id"].astype(str) if "analysis_cell_id" in obs else obs["obs_name"].astype(str)
    consensus_ids = set(consensus["analysis_cell_id"].astype(str))
    boundary = bool_series(stage18, obs["boundary"]) if "boundary" in obs else pd.Series(False, index=obs.index)
    lam_condition = obs["condition"].astype(str).eq("LAM") if "condition" in obs else pd.Series(False, index=obs.index)
    main_mask = obs["analysis_cell_id"].isin(consensus_ids) | (boundary & lam_condition & ~obs["analysis_cell_id"].isin(consensus_ids))
    main_obs = obs.loc[main_mask].copy().reset_index(drop=True)
    main_positions = np.flatnonzero(main_mask.to_numpy())
    latent_main = np.asarray(obj.obsm["X_scVI"][main_positions], dtype=np.float32)
    normal_positions = np.flatnonzero((~lam_condition).to_numpy())
    latent_normal = np.asarray(obj.obsm["X_scVI"][normal_positions], dtype=np.float32)
    obj.file.close()
    state_map = dict(
        zip(
            consensus["analysis_cell_id"].astype(str),
            consensus["consensus_state"].astype(str).map(lambda value: f"State_{value}"),
        )
    )
    main_obs["current_state"] = main_obs["analysis_cell_id"].map(state_map).fillna("")
    main_obs["analysis_role"] = np.where(main_obs["current_state"].ne(""), "primary_candidate", "boundary")
    main_obs["patient"] = main_obs.get("patient_id", "").astype(str)
    main_obs["cell_id"] = main_obs["analysis_cell_id"].astype(str)
    if len(main_obs) != EXPECTED_CANDIDATE_CELLS + EXPECTED_BOUNDARY_CELLS:
        raise ValueError(f"Expected {EXPECTED_CANDIDATE_CELLS + EXPECTED_BOUNDARY_CELLS} main cells, found {len(main_obs)}")
    if int(main_obs["current_state"].eq(f"State_{TARGET_STATE}").sum()) != EXPECTED_STATE15_CELLS:
        raise ValueError("State 15 is not exactly 200 cells in the frozen input")
    manifest = {
        "main_cell_count": len(main_obs),
        "candidate_cell_count": int(main_obs["current_state"].ne("").sum()),
        "boundary_cell_count": int(main_obs["current_state"].eq("").sum()),
        "normal_remote_cell_count": len(latent_normal),
        "latent_shape": list(latent_main.shape),
    }
    return main_obs, latent_main, obs, latent_normal, manifest


def compute_geometry(main_obs: pd.DataFrame, latent: np.ndarray, reference_mask: np.ndarray) -> pd.DataFrame:
    reference = latent[reference_mask]
    if len(reference) != EXPECTED_STATE15_CELLS:
        raise ValueError(f"State 15 reference must contain {EXPECTED_STATE15_CELLS} cells")
    reference_neighbors = NearestNeighbors(n_neighbors=min(15, len(reference)), metric="euclidean").fit(reference)
    ref_distances, _ = reference_neighbors.kneighbors(latent)
    graph_k = min(31, len(latent))
    graph_neighbors = NearestNeighbors(n_neighbors=graph_k, metric="euclidean").fit(latent)
    graph_distances, graph_indices = graph_neighbors.kneighbors(latent)
    state_values = main_obs["current_state"].astype(str).to_numpy()
    state15_neighbor_count = (state_values[graph_indices[:, 1:]] == f"State_{TARGET_STATE}").sum(axis=1)
    reference_centroid = reference.mean(axis=0)
    result = main_obs[["cell_id", "analysis_cell_id", "current_state", "patient", "dataset", "analysis_role"]].copy()
    result["nearest_state15_distance"] = ref_distances[:, 0]
    result["mean_5_nearest_state15_distance"] = ref_distances[:, :5].mean(axis=1)
    result["mean_15_nearest_state15_distance"] = ref_distances[:, :15].mean(axis=1)
    result["state15_centroid_distance"] = np.linalg.norm(latent - reference_centroid, axis=1)
    result["state15_neighbor_count_k30"] = state15_neighbor_count
    result["state15_neighbor_fraction"] = state15_neighbor_count / max(graph_k - 1, 1)
    result["state15_distance"] = result["nearest_state15_distance"]
    result["distance_to_state15"] = result["nearest_state15_distance"]
    return result, graph_indices, graph_distances


def assign_distance_bins(distances: pd.Series) -> pd.Series:
    rank = distances.rank(method="first", pct=True).to_numpy(dtype=float)
    labels = np.asarray(DISTANCE_BIN_LABELS, dtype=object)
    indices = np.clip(np.searchsorted(DISTANCE_BIN_EDGES[1:], rank, side="right"), 0, len(labels) - 1)
    return pd.Series(labels[indices], index=distances.index, dtype="string")


def state_composition(table: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for label in DISTANCE_BIN_LABELS:
        sub = table[table["distance_bin"].eq(label)]
        counts = sub["current_state"].replace("", "boundary").value_counts()
        for state, cells in counts.items():
            rows.append(
                {
                    "distance_bin": label,
                    "current_state": str(state),
                    "cells": int(cells),
                    "fraction_in_bin": float(cells / max(len(sub), 1)),
                    "bin_cells": len(sub),
                }
            )
    return pd.DataFrame(rows)


def gradient_table(table: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for label in DISTANCE_BIN_LABELS:
        sub = table[table["distance_bin"].eq(label)]
        row: dict[str, Any] = {
            "distance_bin": label,
            "n_cells": len(sub),
            "nearest_state15_distance_median": float(sub["nearest_state15_distance"].median()) if len(sub) else np.nan,
        }
        for feature in features:
            if feature not in sub:
                continue
            values = pd.to_numeric(sub[feature], errors="coerce").dropna()
            row[f"{feature}_median"] = float(values.median()) if len(values) else np.nan
            row[f"{feature}_mean"] = float(values.mean()) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def prepared_scores(
    prepared_path: Path,
    selected: pd.DataFrame,
    modules: dict[str, list[str]],
    stage18: Any,
    block_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    prepared = ad.read_h5ad(prepared_path, backed="r")
    scores, manifest = stage18.selected_score_table(prepared, selected[["analysis_cell_id"]], modules, block_size)
    prepared.file.close()
    return selected.merge(scores, on="analysis_cell_id", how="left", validate="one_to_one"), manifest


def extract_raw_marker_counts(
    prepared_path: Path,
    selected: pd.DataFrame,
    genes: list[str],
    stage18: Any,
) -> tuple[np.ndarray, list[str]]:
    prepared = ad.read_h5ad(prepared_path, backed="r")
    counts, present = stage18.extract_pseudobulk_counts(prepared, selected[["analysis_cell_id"]], genes)
    prepared.file.close()
    return counts.to_numpy(dtype=float), present


def state16_coexpression(
    scores: pd.DataFrame,
    raw_counts: np.ndarray,
    raw_genes: list[str],
    state16_mask: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    state16 = scores.loc[state16_mask].copy().reset_index(drop=True)
    state16_counts = raw_counts[state16_mask]
    gene_index = {gene: i for i, gene in enumerate(raw_genes)}
    lam_genes = ["PMEL", "MLANA", "VEGFD", "ACTA2"]
    immune_genes = ["NKG7", "GNLY", "GZMB", "CD3D", "CD3E"]
    lam_detected = np.zeros(len(state16), dtype=bool)
    immune_detected = np.zeros(len(state16), dtype=bool)
    for gene in lam_genes:
        if gene in gene_index:
            lam_detected |= state16_counts[:, gene_index[gene]] > 0
    for gene in immune_genes:
        if gene in gene_index:
            immune_detected |= state16_counts[:, gene_index[gene]] > 0
    lam_identity = state16[["LAMCORE_777", "CORE2", "CORE3", "melanocytic"]].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=True)
    immune_score = state16[["T_NK", "macrophage"]].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=True)
    state15 = scores[scores["current_state"].astype(str).eq(f"State_{TARGET_STATE}")]
    lam_threshold = float(
        state15[["LAMCORE_777", "CORE2", "CORE3", "melanocytic"]].apply(pd.to_numeric, errors="coerce").mean(axis=1, skipna=True).median()
    )
    immune_threshold = float(immune_score.quantile(0.75)) if len(immune_score) else np.nan
    lam_high = lam_identity >= lam_threshold
    immune_high = immune_score >= immune_threshold
    category = np.select(
        [lam_high & ~immune_high, lam_high & immune_high, ~lam_high & immune_high],
        ["LAM-high / immune-low", "LAM-high / immune-high", "LAM-low / immune-high"],
        default="LAM-low / immune-low",
    )
    rows: list[dict[str, Any]] = []
    for value in ["LAM-high / immune-low", "LAM-high / immune-high", "LAM-low / immune-high", "LAM-low / immune-low"]:
        n = int((category == value).sum())
        rows.append(
            {
                "row_type": "category",
                "category": value,
                "marker_pair": "",
                "cells": n,
                "fraction_state16": float(n / max(len(state16), 1)),
                "lam_identity_threshold": lam_threshold,
                "immune_score_threshold": immune_threshold,
            }
        )
    for lam_gene in lam_genes:
        for immune_gene in immune_genes:
            if lam_gene not in gene_index or immune_gene not in gene_index:
                continue
            both = (state16_counts[:, gene_index[lam_gene]] > 0) & (state16_counts[:, gene_index[immune_gene]] > 0)
            rows.append(
                {
                    "row_type": "raw_count_marker_pair",
                    "category": "",
                    "marker_pair": f"{lam_gene}+{immune_gene}",
                    "cells": int(both.sum()),
                    "fraction_state16": float(both.mean()) if len(both) else np.nan,
                    "lam_identity_threshold": lam_threshold,
                    "immune_score_threshold": immune_threshold,
                }
            )
    coexpression = pd.DataFrame(rows)
    state16["lam_identity_score"] = lam_identity.to_numpy()
    state16["immune_score"] = immune_score.to_numpy()
    state16["LAM_marker_any_raw_gt0"] = lam_detected
    state16["immune_marker_any_raw_gt0"] = immune_detected
    state16["LAM_immune_raw_coexpressing"] = lam_detected & immune_detected
    state16["LAM_immune_category"] = category
    coexpression.attrs["state16_per_cell"] = state16
    return coexpression, {
        "state16_cells": len(state16),
        "lam_identity_threshold": lam_threshold,
        "immune_score_threshold": immune_threshold,
        "raw_lam_genes": lam_genes,
        "raw_immune_genes": immune_genes,
        "raw_count_coexpressing_cells": int((lam_detected & immune_detected).sum()),
    }


def technical_audit(scores: pd.DataFrame, coexpression_cells: pd.DataFrame) -> pd.DataFrame:
    main = scores.copy()
    main["technical_cohort"] = "other_candidate_or_boundary"
    main.loc[main["current_state"].astype(str).eq(f"State_{TARGET_STATE}"), "technical_cohort"] = "State_15"
    main.loc[main["current_state"].astype(str).eq("State_16"), "technical_cohort"] = "State_16"
    tnk_threshold = float(pd.to_numeric(main["T_NK"], errors="coerce").quantile(0.75))
    tnk_reference = (
        main["T_NK"].astype(float).ge(tnk_threshold)
        & main["T_NK"].astype(float).ge(main["macrophage"].astype(float))
        & ~main["technical_cohort"].isin(["State_15", "State_16"])
    )
    main.loc[tnk_reference, "technical_cohort"] = "T_NK_marker_reference"
    rows: list[dict[str, Any]] = []
    for cohort, sub in main.groupby("technical_cohort", observed=True):
        row: dict[str, Any] = {
            "technical_cohort": str(cohort),
            "n_cells": len(sub),
            "T_NK_reference_threshold": tnk_threshold,
        }
        for feature in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "doublet_score"]:
            if feature not in sub:
                continue
            values = pd.to_numeric(sub[feature], errors="coerce").dropna()
            row[f"{feature}_median"] = float(values.median()) if len(values) else np.nan
            row[f"{feature}_q25"] = float(values.quantile(0.25)) if len(values) else np.nan
            row[f"{feature}_q75"] = float(values.quantile(0.75)) if len(values) else np.nan
        if "doublet_predicted" in sub:
            predicted = bool_series_for_values(sub["doublet_predicted"])
            row["doublet_predicted_fraction"] = float(predicted.mean())
        rows.append(row)
    coexpress_ids = set(coexpression_cells.loc[coexpression_cells["LAM_immune_raw_coexpressing"], "analysis_cell_id"].astype(str))
    state16 = scores[scores["current_state"].astype(str).eq("State_16")].copy()
    state16["coexpression_group"] = np.where(state16["analysis_cell_id"].astype(str).isin(coexpress_ids), "LAM_immune_raw_coexpressing", "other_State_16")
    for cohort, sub in state16.groupby("coexpression_group", observed=True):
        rows.append(
            {
                "technical_cohort": str(cohort),
                "n_cells": len(sub),
                "T_NK_reference_threshold": tnk_threshold,
                "total_counts_median": float(pd.to_numeric(sub["total_counts"], errors="coerce").median()),
                "n_genes_by_counts_median": float(pd.to_numeric(sub["n_genes_by_counts"], errors="coerce").median()),
                "pct_counts_mt_median": float(pd.to_numeric(sub["pct_counts_mt"], errors="coerce").median()),
                "doublet_score_median": float(pd.to_numeric(sub["doublet_score"], errors="coerce").median()),
            }
        )
    return pd.DataFrame(rows)


def bool_series_for_values(values: pd.Series) -> pd.Series:
    return values.map(lambda value: str(value).strip().lower() in {"true", "1", "yes", "y"}).astype(bool)


def patient_or_dataset_gradients(table: pd.DataFrame, features: list[str], group_field: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for group, sub in table.groupby(group_field, observed=True):
        row: dict[str, Any] = {group_field: str(group), "n_cells": len(sub)}
        distances = pd.to_numeric(sub["nearest_state15_distance"], errors="coerce")
        for feature in features:
            values = pd.to_numeric(sub[feature], errors="coerce")
            valid = distances.notna() & values.notna()
            if int(valid.sum()) >= 3:
                rho, pvalue = spearmanr(distances[valid], values[valid])
            else:
                rho, pvalue = np.nan, np.nan
            row[f"distance_spearman_{feature}"] = float(rho) if np.isfinite(rho) else np.nan
            row[f"distance_spearman_p_{feature}"] = float(pvalue) if np.isfinite(pvalue) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def boundary_projection(table: pd.DataFrame, candidate_latent: np.ndarray, candidate_states: np.ndarray, boundary_latent: np.ndarray) -> pd.DataFrame:
    nearest_candidate = NearestNeighbors(n_neighbors=1, metric="euclidean").fit(candidate_latent)
    distances, indices = nearest_candidate.kneighbors(boundary_latent)
    boundary = table[table["analysis_role"].eq("boundary")].copy().reset_index(drop=True)
    boundary["nearest_current_state"] = candidate_states[indices[:, 0]]
    boundary["nearest_candidate_distance"] = distances[:, 0]
    return boundary


def normal_remote_summary(scvi_obs: pd.DataFrame, latent_normal: np.ndarray, reference_latent: np.ndarray) -> pd.DataFrame:
    if len(latent_normal) == 0:
        return pd.DataFrame([{"cohort": "normal_remote", "n_cells": 0}])
    nearest = NearestNeighbors(n_neighbors=1, metric="euclidean").fit(reference_latent)
    distances = nearest.kneighbors(latent_normal, return_distance=True)[0].ravel()
    centroid = reference_latent.mean(axis=0)
    row: dict[str, Any] = {
        "cohort": "normal_remote",
        "n_cells": len(latent_normal),
        "n_datasets": scvi_obs.loc[~scvi_obs["condition"].astype(str).eq("LAM"), "dataset"].astype(str).nunique() if "dataset" in scvi_obs else np.nan,
        "n_patients": scvi_obs.loc[~scvi_obs["condition"].astype(str).eq("LAM"), "patient_id"].astype(str).nunique() if "patient_id" in scvi_obs else np.nan,
        "nearest_state15_distance_median": float(np.median(distances)),
        "nearest_state15_distance_q25": float(np.quantile(distances, 0.25)),
        "nearest_state15_distance_q75": float(np.quantile(distances, 0.75)),
        "nearest_state15_distance_mean": float(np.mean(distances)),
        "state15_centroid_distance_median": float(np.median(np.linalg.norm(latent_normal - centroid, axis=1))),
    }
    return pd.DataFrame([row])


def write_report(
    output_dir: Path,
    manifest: dict[str, Any],
    geometry: pd.DataFrame,
    identity_gradient: pd.DataFrame,
    lineage_gradient: pd.DataFrame,
    patient_gradient: pd.DataFrame,
    dataset_gradient: pd.DataFrame,
    coexpression_manifest: dict[str, Any],
    normal_summary: pd.DataFrame,
) -> None:
    identity_near = identity_gradient.iloc[0] if len(identity_gradient) else pd.Series(dtype=float)
    identity_far = identity_gradient.iloc[-1] if len(identity_gradient) else pd.Series(dtype=float)
    lamcore_near = float(identity_near.get("LAMCORE_777_median", np.nan))
    lamcore_far = float(identity_far.get("LAMCORE_777_median", np.nan))
    patient_lamcore = pd.to_numeric(patient_gradient.get("distance_spearman_LAMCORE_777", pd.Series(dtype=float)), errors="coerce").dropna()
    dataset_lamcore = pd.to_numeric(dataset_gradient.get("distance_spearman_LAMCORE_777", pd.Series(dtype=float)), errors="coerce").dropna()
    report = [
        "# Stage 20：State 15-centered transcriptional manifold",
        "",
        "本阶段固定使用现有 State 15 的 200 个细胞和 `X_scVI`，分析对象为 5,378 个 high-confidence candidates 加 16,883 个 inherited boundary cells。没有重训 scVI、重新 Leiden/consensus clustering 或修改 candidate gate。",
        "",
        "## Frozen input",
        "",
        f"- State 15 cells: {manifest['frozen_input']['state15_cell_count']}",
        f"- State 15 ID SHA-256: `{manifest['frozen_input']['state15_cell_id_sha256']}`",
        f"- Latent artifact: `{manifest['frozen_input']['latent_artifact']}` (`{manifest['frozen_input']['latent_key']}`)",
        f"- Main geometry cells: {manifest['cohort']['main_cell_count']} ({manifest['cohort']['candidate_cell_count']} candidates + {manifest['cohort']['boundary_cell_count']} boundary)",
        "",
        "## Geometry-first distance axis",
        "",
        "距离、State15 邻居比例和 30-NN 图均只由 `X_scVI` 计算；LAMCORE、program 和 lineage score 在几何结果生成后才作为描述性映射。",
        "",
        f"- State 15-nearest distance median: {float(geometry.loc[geometry['current_state'].eq('State_15'), 'nearest_state15_distance'].median()):.4f}",
        f"- State 15 neighbor fraction median: {float(geometry.loc[geometry['current_state'].eq('State_15'), 'state15_neighbor_fraction'].median()):.4f}",
        "",
        "## Distance-bin composition",
        "",
        "详见 `state15_distance_bins.csv`；该表只报告当前 State 1–20 标签和 boundary 的几何组成，不重新定义 state。",
        "",
        "## Identity and lineage gradients",
        "",
        f"- LAMCORE median, nearest bin → farthest bin: {lamcore_near:.4f} → {lamcore_far:.4f}.",
        f"- Patient-level distance~LAMCORE Spearman rho < 0 in {int((patient_lamcore < 0).sum())}/{len(patient_lamcore)} patients with estimable values.",
        f"- Dataset-level distance~LAMCORE Spearman rho < 0 in {int((dataset_lamcore < 0).sum())}/{len(dataset_lamcore)} datasets with estimable values.",
        "",
        "## State 16 audit",
        "",
        f"- State 16 cells audited: {coexpression_manifest.get('state16_cells', 0)}.",
        f"- Raw-count LAM-marker/immune-marker coexpressing State 16 cells: {coexpression_manifest.get('raw_count_coexpressing_cells', 0)}.",
        f"- LAM identity threshold: {coexpression_manifest.get('lam_identity_threshold', np.nan):.4f}; immune threshold: {coexpression_manifest.get('immune_score_threshold', np.nan):.4f}.",
        "详见 `state16_lam_immune_coexpression.csv` 和 `state16_doublet_audit.csv`；阈值只用于诊断分类，不改变 State 16 标签。",
        "",
        "## Boundary and normal scope",
        "",
        "Boundary 仅投射到 State 15-centered geometry，normal 仅作为远端对照，不参与 State 15 邻域图或 state 数量。",
        "",
        normal_summary.to_string(index=False),
        "",
        "## Stage 20 checkpoint",
        "",
        f"- Current geometry checkpoint: `{manifest['checkpoint']}`.",
        f"- Interpretation: {manifest['checkpoint_interpretation']}",
        "",
        "## Outputs",
        "",
        *[f"- {path.name}" for path in sorted(output_dir.iterdir()) if path.is_file()],
    ]
    (output_dir / "stage20_manifold_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/state_modeling.yaml"))
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "results/stage20"))
    parser.add_argument("--block-size", type=int, default=8192)
    args = parser.parse_args()
    config = load_config(Path(args.config).resolve())
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stage18 = load_stage18_module()
    modules, formal_manifest, program_manifest = load_modules(stage18, config)

    consensus_path = PROJECT_ROOT / str(config["outputs"]["consensus_upstream_annotations"])
    prepared_path = PROJECT_ROOT / str(config["outputs"]["prepared_h5ad"])
    scvi_path = PROJECT_ROOT / str(config["outputs"]["scvi_h5ad"])
    if not consensus_path.exists() or not prepared_path.exists() or not scvi_path.exists():
        raise FileNotFoundError(f"Required inputs missing: {consensus_path}, {prepared_path}, {scvi_path}")
    consensus = pd.read_csv(consensus_path)
    consensus["analysis_cell_id"] = consensus["analysis_cell_id"].astype(str)
    consensus["consensus_state"] = consensus["consensus_state"].astype(str)
    if consensus["analysis_cell_id"].duplicated().any():
        raise ValueError("Duplicate analysis_cell_id in consensus annotation")
    if int(consensus["consensus_state"].eq(TARGET_STATE).sum()) != EXPECTED_STATE15_CELLS:
        raise ValueError("Frozen State 15 is not exactly 200 cells")

    main_obs, latent_main, scvi_obs, latent_normal, cohort_manifest = load_scvi_cohort(scvi_path, consensus, stage18)
    frozen_manifest = fixed_state15_manifest(consensus, scvi_path, tuple(latent_main.shape))
    reference_mask = main_obs["current_state"].astype(str).eq(f"State_{TARGET_STATE}").to_numpy()
    geometry, graph_indices, graph_distances = compute_geometry(main_obs, latent_main, reference_mask)
    geometry["distance_bin"] = assign_distance_bins(geometry["nearest_state15_distance"])

    prepared = ad.read_h5ad(prepared_path, backed="r")
    prepared_obs_table = prepared_metadata(stage18, prepared)
    prepared.file.close()
    selected_for_score = main_obs[["analysis_cell_id", "cell_id", "current_state", "patient", "dataset", "analysis_role"]].copy()
    for column in ["total_counts", "n_genes_by_counts", "pct_counts_mt", "doublet_score", "doublet_predicted", "assay", "condition"]:
        if column in main_obs:
            selected_for_score[column] = main_obs[column].to_numpy()
        elif column in prepared_obs_table:
            values = prepared_obs_table.set_index("analysis_cell_id")[column]
            selected_for_score[column] = selected_for_score["analysis_cell_id"].map(values).to_numpy()
    scores, score_manifest = prepared_scores(prepared_path, selected_for_score, modules, stage18, int(args.block_size))
    geometry = geometry.merge(scores.drop(columns=["cell_id", "current_state", "patient", "dataset", "analysis_role"], errors="ignore"), on=["analysis_cell_id"], how="left", validate="one_to_one")
    geometry["current_state"] = geometry["current_state"].astype(str)
    geometry["analysis_role"] = geometry["analysis_role"].astype(str)
    geometry["patient"] = geometry["patient"].astype(str)
    geometry["dataset"] = geometry["dataset"].astype(str)

    composition = state_composition(geometry)
    composition.to_csv(output_dir / "state15_distance_bins.csv", index=False)
    identity_features = [
        "LAMCORE_777", "CORE1", "CORE2", "CORE3", "melanocytic", "lam_support",
        "myogenic", "HOX_PBX", "hormone", "ECM", "protease", "mTOR",
    ]
    lineage_features = [
        "macrophage", "T_NK", "AT2", "endothelial", "lymphatic_endothelial",
        "fibroblast", "pericyte_VSMC", "mesothelial", "ciliated",
    ]
    identity_gradient = gradient_table(geometry, identity_features)
    lineage_gradient = gradient_table(geometry, lineage_features)
    identity_gradient.to_csv(output_dir / "state15_identity_gradient.csv", index=False)
    lineage_gradient.to_csv(output_dir / "state15_lineage_gradient.csv", index=False)

    state16_mask = geometry["current_state"].astype(str).eq("State_16").to_numpy()
    state16_audit = geometry.loc[state16_mask].sort_values("nearest_state15_distance")
    state16_audit.to_csv(output_dir / "state16_cell_audit.csv", index=False)
    raw_marker_genes = [
        "PMEL", "MLANA", "VEGFD", "ACTA2", "NKG7", "GNLY", "GZMB", "CD3D", "CD3E",
    ]
    raw_counts, raw_genes = extract_raw_marker_counts(prepared_path, main_obs, raw_marker_genes, stage18)
    coexpression, coexpression_manifest = state16_coexpression(scores, raw_counts, raw_genes, state16_mask)
    coexpression.to_csv(output_dir / "state16_lam_immune_coexpression.csv", index=False)
    coexpression_cells = coexpression.attrs.get("state16_per_cell", pd.DataFrame())
    if len(coexpression_cells):
        coexpression_cells.to_csv(output_dir / "state16_lam_immune_per_cell.csv", index=False)
    technical = technical_audit(scores, coexpression_cells if len(coexpression_cells) else pd.DataFrame(columns=["analysis_cell_id", "LAM_immune_raw_coexpressing"]))
    technical.to_csv(output_dir / "state16_doublet_audit.csv", index=False)

    patient_features = ["LAMCORE_777", "CORE1", "CORE3", "melanocytic", "lam_support", "myogenic", "HOX_PBX", *lineage_features]
    patient_gradient = patient_or_dataset_gradients(geometry, patient_features, "patient")
    dataset_gradient = patient_or_dataset_gradients(geometry, patient_features, "dataset")
    patient_gradient.to_csv(output_dir / "patient_gradient_consistency.csv", index=False)
    dataset_gradient.to_csv(output_dir / "dataset_gradient_consistency.csv", index=False)

    candidate_mask = main_obs["current_state"].astype(str).ne("").to_numpy()
    candidate_latent = latent_main[candidate_mask]
    candidate_states = main_obs.loc[candidate_mask, "current_state"].astype(str).to_numpy()
    boundary_projection_table = boundary_projection(geometry, candidate_latent, candidate_states, latent_main[~candidate_mask])
    boundary_projection_table.to_csv(output_dir / "boundary_state15_projection.csv", index=False)

    reference_latent = latent_main[reference_mask]
    normal_summary = normal_remote_summary(scvi_obs, latent_normal, reference_latent)
    normal_summary.to_csv(output_dir / "normal_remote_summary.csv", index=False)

    manifold_columns = [
        "cell_id", "analysis_cell_id", "patient", "dataset", "current_state", "analysis_role",
        "nearest_state15_distance", "mean_5_nearest_state15_distance", "mean_15_nearest_state15_distance",
        "state15_centroid_distance", "state15_neighbor_fraction", "state15_neighbor_count_k30",
        "state15_distance", "distance_to_state15", *identity_features, *lineage_features,
    ]
    manifold = geometry[[column for column in manifold_columns if column in geometry.columns]].copy()
    manifold.to_csv(output_dir / "state15_centered_manifold.csv", index=False)
    geometry.to_csv(output_dir / "state15_cell_distances.csv", index=False)

    near = identity_gradient.iloc[0] if len(identity_gradient) else pd.Series(dtype=float)
    far = identity_gradient.iloc[-1] if len(identity_gradient) else pd.Series(dtype=float)
    lamcore_near = float(near.get("LAMCORE_777_median", np.nan))
    lamcore_far = float(far.get("LAMCORE_777_median", np.nan))
    patient_rho = pd.to_numeric(patient_gradient.get("distance_spearman_LAMCORE_777", pd.Series(dtype=float)), errors="coerce").dropna()
    dataset_rho = pd.to_numeric(dataset_gradient.get("distance_spearman_LAMCORE_777", pd.Series(dtype=float)), errors="coerce").dropna()
    lamcore_gradient_support = bool(np.isfinite(lamcore_near) and np.isfinite(lamcore_far) and lamcore_near > lamcore_far)
    patient_direction_support = bool(len(patient_rho) and float((patient_rho < 0).mean()) > 0.5)
    dataset_direction_support = bool(len(dataset_rho) and float((dataset_rho < 0).mean()) > 0.5)
    if lamcore_gradient_support and patient_direction_support and dataset_direction_support:
        checkpoint = "supports_lam_centered_transcriptional_manifold"
        interpretation = "LAM identity signals decline along the State 15 distance axis with majority patient- and dataset-level negative distance correlations."
    elif lamcore_gradient_support:
        checkpoint = "pooled_or_partial_lam_centered_gradient"
        interpretation = "The pooled distance axis retains an identity gradient, but patient/dataset direction is not sufficiently concordant for a uniform manifold claim."
    else:
        checkpoint = "no_monotonic_lam_centered_gradient"
        interpretation = "The current geometry does not show a monotonic pooled LAMCORE decline from State 15; retain a discrete or patient-specific interpretation."
    manifest = {
        "stage": 20,
        "frozen_input": frozen_manifest,
        "cohort": cohort_manifest,
        "score_manifest": score_manifest,
        "formal_signature_manifest": formal_manifest,
        "program_manifest": program_manifest,
        "distance_geometry": {
            "reference": "frozen State 15 cells only",
            "candidate_boundary_graph_k": 30,
            "distance_features_not_used": ["LAMCORE", "CORE1", "CORE2", "CORE3", "marker_scores", "lineage_scores"],
            "distance_bin_labels": DISTANCE_BIN_LABELS,
        },
        "state16_coexpression": coexpression_manifest,
        "normal_remote_summary": normal_summary.to_dict(orient="records"),
        "checkpoint": checkpoint,
        "checkpoint_interpretation": interpretation,
        "no_scvi_training": True,
        "no_reclustering": True,
        "no_candidate_gate_change": True,
    }
    (output_dir / "stage20_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    write_report(output_dir, manifest, geometry, identity_gradient, lineage_gradient, patient_gradient, dataset_gradient, coexpression_manifest, normal_summary)
    print(f"State 15-centered geometry: {len(geometry)} cells")
    print(f"State 16 audited: {int(state16_mask.sum())} cells")
    print(f"Checkpoint: {checkpoint}")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
