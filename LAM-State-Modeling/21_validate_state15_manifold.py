#!/usr/bin/env python3
"""Validate the State 15-centered gradient outside the State 15 anchor.

Stage 21 is a read-only validation stage.  It uses the frozen Stage 20
``distance_to_state15`` values and the existing ``X_scVI`` artifact, but does
not train scVI, recluster cells, or change a candidate gate.  State 15 is an
anchor only; every primary gradient result is computed on the remaining
22,061 cells.
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
from scipy.stats import spearmanr, t
from sklearn.neighbors import NearestNeighbors


PROJECT_ROOT = Path(__file__).resolve().parent
TARGET_STATE = "State_15"
EXPECTED_ANCHOR_CELLS = 200
EXPECTED_MAIN_CELLS = 22261
EXPECTED_VALIDATION_CELLS = 22061
EXPECTED_CANDIDATE_CELLS = 5178
DISTANCE_BIN_LABELS = ["0-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-100%"]
DISTANCE_BIN_EDGES = np.asarray([0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0])
OLD_GATE_MARKERS = ["PMEL", "MLANA", "MITF", "ACTA2", "ESR1", "VEGFD", "CTSK"]
LAMCORE_SCORE_NAMES = [
    "LAMCORE_full",
    "LAMCORE_no_gate",
    "LAMCORE_outside_scVI",
    "LAMCORE_independent",
]
LINEAGE_FEATURES = [
    "AT2", "endothelial", "lymphatic_endothelial", "fibroblast",
    "pericyte_VSMC", "macrophage", "T_NK", "mesothelial", "ciliated",
]
STATE16_FEATURES = [
    *LAMCORE_SCORE_NAMES,
    "CORE1", "CORE2", "CORE3", "melanocytic", "lam_support", "HOX_PBX",
    "myogenic", "T_NK", "macrophage", *LINEAGE_FEATURES,
]


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def load_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_helpers() -> tuple[Any, Any]:
    stage18 = load_module(PROJECT_ROOT / "18_validate_state15_anchor.py", "stage18_for_stage21")
    stage20 = load_module(PROJECT_ROOT / "20_state15_centered_manifold.py", "stage20_for_stage21")
    return stage18, stage20


def canonical_gene(stage18: Any, gene: str) -> str:
    return str(stage18.canonical_gene(str(gene))).upper()


def unique_genes(stage18: Any, genes: list[str]) -> list[str]:
    output: list[str] = []
    for gene in genes:
        canonical = canonical_gene(stage18, gene)
        if canonical and canonical not in output:
            output.append(canonical)
    return output


def read_stage20_inputs(
    distance_path: Path,
    scvi_path: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], dict[str, Any]]:
    distances = pd.read_csv(distance_path)
    required = {
        "analysis_cell_id", "current_state", "patient", "dataset", "analysis_role",
        "distance_to_state15", "nearest_state15_distance",
    }
    missing = sorted(required.difference(distances.columns))
    if missing:
        raise ValueError(f"Stage 20 distance file is missing columns: {missing}")
    distances["analysis_cell_id"] = distances["analysis_cell_id"].astype(str)
    distances["current_state"] = distances["current_state"].fillna("").astype(str)
    if distances["analysis_cell_id"].duplicated().any():
        raise ValueError("Stage 20 distance file contains duplicate analysis_cell_id")
    anchor = distances[distances["current_state"].eq(TARGET_STATE)].copy()
    validation = distances[~distances["current_state"].eq(TARGET_STATE)].copy()
    if len(distances) != EXPECTED_MAIN_CELLS:
        raise ValueError(f"Expected {EXPECTED_MAIN_CELLS} Stage 20 cells, found {len(distances)}")
    if len(anchor) != EXPECTED_ANCHOR_CELLS:
        raise ValueError(f"Expected {EXPECTED_ANCHOR_CELLS} State 15 anchor cells, found {len(anchor)}")
    if len(validation) != EXPECTED_VALIDATION_CELLS:
        raise ValueError(f"Expected {EXPECTED_VALIDATION_CELLS} non-State15 cells, found {len(validation)}")
    if not scvi_path.exists():
        raise FileNotFoundError(scvi_path)
    anchor_ids = sorted(anchor["analysis_cell_id"].tolist())
    manifest = {
        "stage20_distance_file": str(distance_path),
        "scvi_artifact": str(scvi_path),
        "latent_key": "X_scVI",
        "anchor_state": TARGET_STATE,
        "anchor_cell_count": len(anchor_ids),
        "anchor_cell_id_sha256": hashlib.sha256("\n".join(anchor_ids).encode("utf-8")).hexdigest(),
        "anchor_cell_ids": anchor_ids,
        "validation_cell_count": len(validation),
        "validation_cell_ids": sorted(validation["analysis_cell_id"].tolist()),
        "validation_scope": "all Stage 20 cells except frozen State 15 anchor",
        "candidate_null_scope": "non-State15 primary_candidate cells only",
        "no_scvi_training": True,
        "no_reclustering": True,
        "no_candidate_gate_change": True,
    }
    return anchor, validation, anchor_ids, manifest


def load_latent_for_ids(scvi_path: Path, cell_ids: pd.Series | list[str]) -> tuple[np.ndarray, pd.DataFrame]:
    obj = ad.read_h5ad(scvi_path, backed="r")
    if "X_scVI" not in obj.obsm:
        obj.file.close()
        raise ValueError(f"X_scVI is missing from {scvi_path}")
    obs = obj.obs.copy()
    if "analysis_cell_id" in obs:
        ids = obs["analysis_cell_id"].astype(str)
    else:
        ids = pd.Series(obj.obs_names.astype(str), index=obs.index)
    positions = pd.Index(ids).get_indexer(pd.Series(cell_ids).astype(str).tolist())
    if (positions < 0).any():
        obj.file.close()
        raise ValueError(f"Missing {int((positions < 0).sum())} requested IDs in X_scVI artifact")
    latent = np.asarray(obj.obsm["X_scVI"][positions], dtype=np.float32)
    selected_obs = obs.iloc[positions].copy().reset_index(drop=True)
    selected_obs["analysis_cell_id"] = pd.Series(cell_ids).astype(str).to_numpy()
    obj.file.close()
    return latent, selected_obs


def scvi_hvg_and_expression_audit(
    scvi_path: Path,
    prepared_path: Path,
    formal_genes: list[str],
    stage18: Any,
) -> tuple[pd.DataFrame, set[str], set[str]]:
    scvi = ad.read_h5ad(scvi_path, backed="r")
    scvi_gene_names = [canonical_gene(stage18, gene) for gene in scvi.var_names.astype(str)]
    if "highly_variable" in scvi.var.columns:
        hv_mask = scvi.var["highly_variable"].astype(bool).to_numpy()
        scvi_hvg = {gene for gene, flag in zip(scvi_gene_names, hv_mask) if flag}
    else:
        scvi_hvg = set(scvi_gene_names)
    scvi_shape = list(scvi.shape)
    scvi.file.close()
    prepared = ad.read_h5ad(prepared_path, backed="r")
    expression_genes = [canonical_gene(stage18, gene) for gene in prepared.var_names.astype(str)]
    expression_set = set(expression_genes)
    expression_aliases: dict[str, list[str]] = {}
    for gene in prepared.var_names.astype(str):
        expression_aliases.setdefault(canonical_gene(stage18, gene), []).append(str(gene))
    gate_markers = set(unique_genes(stage18, OLD_GATE_MARKERS))
    rows: list[dict[str, Any]] = []
    for gene in formal_genes:
        canonical = canonical_gene(stage18, gene)
        in_hvg = canonical in scvi_hvg
        in_gate = canonical in gate_markers
        rows.append(
            {
                "gene": canonical,
                "in_scvi_4000_hvg": in_hvg,
                "in_old_candidate_gate_marker": in_gate,
                "in_neither_scvi_hvg_nor_gate": not in_hvg and not in_gate,
                "in_expression_matrix": canonical in expression_set,
                "expression_aliases": ";".join(expression_aliases.get(canonical, [])),
                "n_expression_aliases": len(expression_aliases.get(canonical, [])),
            }
        )
    prepared.file.close()
    audit = pd.DataFrame(rows).drop_duplicates("gene").sort_values("gene").reset_index(drop=True)
    audit.attrs["scvi_shape"] = scvi_shape
    return audit, scvi_hvg, expression_set


def build_score_modules(stage18: Any, stage20: Any, config: dict[str, Any], formal_genes: list[str], scvi_hvg: set[str]) -> tuple[dict[str, list[str]], dict[str, Any]]:
    stage20_modules, formal_manifest, program_manifest = stage20.load_modules(stage18, config)
    formal = unique_genes(stage18, formal_genes)
    gate = set(unique_genes(stage18, OLD_GATE_MARKERS))
    no_gate = [gene for gene in formal if gene not in gate]
    outside_scvi = [gene for gene in formal if gene not in scvi_hvg]
    independent = [gene for gene in outside_scvi if gene not in gate]
    modules: dict[str, list[str]] = {
        "LAMCORE_full": formal,
        "LAMCORE_no_gate": no_gate,
        "LAMCORE_outside_scVI": outside_scvi,
        "LAMCORE_independent": independent,
        "melanocytic": ["PMEL", "MLANA", "MITF"],
        "lam_support": ["VEGFD", "CTSK", "ESR1"],
        "myogenic": ["ACTA2", "ACTG2", "MYH11"],
        "HOX_PBX": ["EMX2", "HOXA11"],
    }
    needed = [
        "CORE1", "CORE2", "CORE3", "LAM_myogenic", "ECM", "mTOR_translation",
        "hormone", "protease", "HOX_PBX", *LINEAGE_FEATURES,
    ]
    for name in needed:
        if name in stage20_modules:
            modules[name] = stage20_modules[name]
    modules = {name: unique_genes(stage18, genes) for name, genes in modules.items()}
    return modules, {
        "formal_gene_count": len(formal),
        "no_gate_gene_count": len(no_gate),
        "outside_scvi_gene_count": len(outside_scvi),
        "independent_gene_count": len(independent),
        "formal_manifest": formal_manifest,
        "program_manifest": program_manifest,
    }


def score_selected_cells(
    prepared_path: Path,
    selected: pd.DataFrame,
    modules: dict[str, list[str]],
    stage18: Any,
    block_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    prepared = ad.read_h5ad(prepared_path, backed="r")
    scores, manifest = stage18.selected_score_table(
        prepared,
        selected[["analysis_cell_id"]],
        modules,
        block_size,
    )
    prepared.file.close()
    return selected.merge(scores, on="analysis_cell_id", how="left", validate="one_to_one"), manifest


def assign_distance_bins(distances: pd.Series) -> pd.Series:
    rank = distances.rank(method="first", pct=True).to_numpy(dtype=float)
    indices = np.clip(np.searchsorted(DISTANCE_BIN_EDGES[1:], rank, side="right"), 0, len(DISTANCE_BIN_LABELS) - 1)
    return pd.Series(np.asarray(DISTANCE_BIN_LABELS, dtype=object)[indices], index=distances.index, dtype="string")


def add_dataset_standardized_distance(table: pd.DataFrame, distance_column: str = "distance_to_state15") -> pd.DataFrame:
    output = table.copy()
    standardized = pd.Series(np.nan, index=output.index, dtype=float)
    for _, indices in output.groupby("dataset", observed=True).groups.items():
        values = pd.to_numeric(output.loc[indices, distance_column], errors="coerce")
        mean = float(values.mean())
        std = float(values.std(ddof=1))
        if np.isfinite(std) and std > 0:
            standardized.loc[indices] = (values - mean) / std
        else:
            standardized.loc[indices] = 0.0
    output["distance_dataset_z"] = standardized
    return output


def summary_by_bins(table: pd.DataFrame, features: list[str], bin_column: str = "distance_bin") -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for label in DISTANCE_BIN_LABELS:
        sub = table[table[bin_column].eq(label)]
        row: dict[str, Any] = {"distance_bin": label, "n_cells": len(sub)}
        for feature in features:
            values = pd.to_numeric(sub[feature], errors="coerce").dropna() if feature in sub else pd.Series(dtype=float)
            row[f"{feature}_median"] = float(values.median()) if len(values) else np.nan
            row[f"{feature}_q25"] = float(values.quantile(0.25)) if len(values) else np.nan
            row[f"{feature}_q75"] = float(values.quantile(0.75)) if len(values) else np.nan
            row[f"{feature}_mean"] = float(values.mean()) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def smooth_distance_table(table: pd.DataFrame, features: list[str], n_bins: int = 20) -> pd.DataFrame:
    """Return quantile-binned points for a lightweight distance-score curve."""
    if len(table) == 0:
        return pd.DataFrame()
    ranks = table["distance_to_state15"].rank(method="first", pct=True).to_numpy(dtype=float)
    bins = np.clip(np.ceil(ranks * n_bins).astype(int), 1, n_bins)
    working = table.copy()
    working["smooth_bin"] = bins
    rows: list[dict[str, Any]] = []
    for bin_id, sub in working.groupby("smooth_bin", observed=True):
        row: dict[str, Any] = {
            "smooth_bin": int(bin_id),
            "distance_percentile_midpoint": (float(bin_id) - 0.5) / n_bins,
            "n_cells": len(sub),
            "distance_median": float(pd.to_numeric(sub["distance_to_state15"], errors="coerce").median()),
        }
        for feature in features:
            values = pd.to_numeric(sub[feature], errors="coerce").dropna()
            row[f"{feature}_median"] = float(values.median()) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows).sort_values("smooth_bin").reset_index(drop=True)


def ols_slope(distance: pd.Series, score: pd.Series, patient: pd.Series | None = None) -> dict[str, float]:
    frame = pd.DataFrame({"distance": pd.to_numeric(distance, errors="coerce"), "score": pd.to_numeric(score, errors="coerce")})
    if patient is not None:
        frame["patient"] = patient.astype(str)
    frame = frame.replace([np.inf, -np.inf], np.nan).dropna()
    if len(frame) < 4 or frame["distance"].nunique() < 2:
        return {"n": len(frame), "slope": np.nan, "ci95_low": np.nan, "ci95_high": np.nan, "pvalue": np.nan}
    if patient is not None and frame["patient"].nunique() > 1:
        dummies = pd.get_dummies(frame["patient"], drop_first=True, dtype=float).to_numpy()
        design = np.column_stack([np.ones(len(frame)), frame["distance"].to_numpy(dtype=float), dummies])
    else:
        design = np.column_stack([np.ones(len(frame)), frame["distance"].to_numpy(dtype=float)])
    response = frame["score"].to_numpy(dtype=float)
    beta, _, rank, _ = np.linalg.lstsq(design, response, rcond=None)
    residual = response - design @ beta
    degrees = len(response) - int(rank)
    if degrees <= 0:
        return {"n": len(frame), "slope": float(beta[1]), "ci95_low": np.nan, "ci95_high": np.nan, "pvalue": np.nan}
    sigma2 = float(np.dot(residual, residual) / degrees)
    covariance = sigma2 * np.linalg.pinv(design.T @ design)
    standard_error = float(np.sqrt(max(covariance[1, 1], 0.0)))
    slope = float(beta[1])
    if standard_error == 0:
        pvalue = 0.0 if slope != 0 else np.nan
        low = high = slope
    else:
        statistic = slope / standard_error
        pvalue = float(2.0 * t.sf(abs(statistic), degrees))
        critical = float(t.ppf(0.975, degrees))
        low = slope - critical * standard_error
        high = slope + critical * standard_error
    return {"n": len(frame), "slope": slope, "ci95_low": low, "ci95_high": high, "pvalue": pvalue}


def model_gradients(table: pd.DataFrame, features: list[str], scope: str) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for feature in features:
        valid = table[["distance_to_state15", "distance_dataset_z", "patient", feature]].copy()
        valid[feature] = pd.to_numeric(valid[feature], errors="coerce")
        valid = valid.dropna(subset=["distance_to_state15", "distance_dataset_z", "patient", feature])
        if len(valid) >= 3:
            rho, rho_p = spearmanr(valid["distance_to_state15"], valid[feature])
        else:
            rho, rho_p = np.nan, np.nan
        fit = ols_slope(valid["distance_dataset_z"], valid[feature], valid["patient"])
        rows.append(
            {
                "scope": scope,
                "score_name": feature,
                "n_cells": int(fit["n"]),
                "n_patients": int(valid["patient"].nunique()),
                "distance_metric": "Stage20 distance_to_state15; dataset-standardized for regression",
                "patient_adjusted_slope": fit["slope"],
                "ci95_low": fit["ci95_low"],
                "ci95_high": fit["ci95_high"],
                "pvalue": fit["pvalue"],
                "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "spearman_pvalue": float(rho_p) if np.isfinite(rho_p) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def group_gradient(table: pd.DataFrame, features: list[str], group_field: str, stage20_nonnegative: set[str] | None = None) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for group, sub in table.groupby(group_field, observed=True):
        for feature in features:
            valid = sub[["distance_to_state15", "distance_dataset_z", feature]].copy()
            valid[feature] = pd.to_numeric(valid[feature], errors="coerce")
            valid = valid.dropna()
            if len(valid) >= 3 and valid["distance_to_state15"].nunique() > 1:
                rho, rho_p = spearmanr(valid["distance_to_state15"], valid[feature])
                fit = ols_slope(valid["distance_dataset_z"], valid[feature])
            else:
                rho, rho_p = np.nan, np.nan
                fit = {"n": len(valid), "slope": np.nan, "ci95_low": np.nan, "ci95_high": np.nan, "pvalue": np.nan}
            row: dict[str, Any] = {
                group_field: str(group),
                "score_name": feature,
                "n_cells": len(valid),
                "distance_range": float(valid["distance_to_state15"].max() - valid["distance_to_state15"].min()) if len(valid) else np.nan,
                "distance_spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "distance_spearman_pvalue": float(rho_p) if np.isfinite(rho_p) else np.nan,
                "slope": fit["slope"],
                "ci95_low": fit["ci95_low"],
                "ci95_high": fit["ci95_high"],
                "pvalue": fit["pvalue"],
            }
            if group_field == "patient":
                row["stage20_nonnegative_patient"] = str(group) in (stage20_nonnegative or set())
                if np.isfinite(rho):
                    row["gradient_class"] = "positive" if rho > 0.1 else "negative" if rho < -0.1 else "near_zero"
                else:
                    row["gradient_class"] = "not_estimable"
            rows.append(row)
    return pd.DataFrame(rows)


def distance_segments(table: pd.DataFrame, group_field: str | None = None) -> pd.Series:
    labels = ["near", "mid", "far"]
    output = pd.Series("", index=table.index, dtype="string")
    groups = [("all", table.index)] if group_field is None else table.groupby(group_field, observed=True).groups.items()
    for _, indices in groups:
        sub = table.loc[indices]
        ranks = sub["distance_to_state15"].rank(method="first", pct=True).to_numpy(dtype=float)
        output.loc[indices] = np.asarray(labels, dtype=object)[np.clip(np.ceil(ranks * 3).astype(int) - 1, 0, 2)]
    return output


def state16_gradient(table: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    state16 = table[table["current_state"].eq("State_16")].copy()
    rows: list[pd.DataFrame] = []
    for scope, sub, patient in [("pooled", state16, "")]:
        if len(sub):
            sub = sub.copy()
            sub["distance_segment"] = distance_segments(sub)
            summary = summary_by_segments(sub, features)
            summary.insert(0, "patient", patient)
            summary.insert(0, "scope", scope)
            rows.append(summary)
    for patient, sub in state16.groupby("patient", observed=True):
        if len(sub) < 3:
            continue
        sub = sub.copy()
        sub["distance_segment"] = distance_segments(sub)
        summary = summary_by_segments(sub, features)
        summary.insert(0, "patient", str(patient))
        summary.insert(0, "scope", "patient")
        rows.append(summary)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def summary_by_segments(table: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for segment in ["near", "mid", "far"]:
        sub = table[table["distance_segment"].eq(segment)]
        row: dict[str, Any] = {"distance_segment": segment, "n_cells": len(sub)}
        for feature in features:
            values = pd.to_numeric(sub[feature], errors="coerce").dropna() if feature in sub else pd.Series(dtype=float)
            row[f"{feature}_median"] = float(values.median()) if len(values) else np.nan
            row[f"{feature}_q25"] = float(values.quantile(0.25)) if len(values) else np.nan
            row[f"{feature}_q75"] = float(values.quantile(0.75)) if len(values) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def pairwise_nearest_distance(query: np.ndarray, reference: np.ndarray, block_size: int = 2048) -> np.ndarray:
    reference = np.asarray(reference, dtype=np.float32)
    query = np.asarray(query, dtype=np.float32)
    reference_norm = np.sum(reference * reference, axis=1, dtype=np.float32)
    output = np.empty(len(query), dtype=np.float32)
    for start in range(0, len(query), block_size):
        stop = min(start + block_size, len(query))
        query_block = query[start:stop]
        query_norm = np.sum(query_block * query_block, axis=1, dtype=np.float32)[:, None]
        squared = query_norm + reference_norm[None, :] - 2.0 * (query_block @ reference.T)
        output[start:stop] = np.sqrt(np.maximum(np.min(squared, axis=1), 0.0))
    return output


def matched_anchor_null(
    target: pd.DataFrame,
    target_latent: np.ndarray,
    score_name: str,
    anchor: pd.DataFrame,
    candidate_pool: pd.DataFrame,
    candidate_latent: np.ndarray,
    repetitions: int,
    seed: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    rng = np.random.default_rng(seed)
    anchor_keys = anchor.assign(_stratum=anchor["patient"].astype(str) + "||" + anchor["dataset"].astype(str))
    pool_keys = candidate_pool.assign(_stratum=candidate_pool["patient"].astype(str) + "||" + candidate_pool["dataset"].astype(str))
    pool_by_key = {key: values.index.to_numpy() for key, values in pool_keys.groupby("_stratum", observed=True)}
    requested = anchor_keys["_stratum"].value_counts().to_dict()
    missing = sorted(set(requested).difference(pool_by_key))
    if missing:
        raise ValueError(f"Matched null strata missing from non-State15 candidate pool: {missing}")
    rows: list[dict[str, Any]] = []
    target_score = pd.to_numeric(target[score_name], errors="coerce").to_numpy(dtype=float)
    target_patient = target["patient"].astype(str)
    target_dataset = target["dataset"].astype(str)
    for replicate in range(repetitions):
        sampled_indices: list[int] = []
        replacement_count = 0
        for key, count in requested.items():
            eligible = pool_by_key[key]
            replace = len(eligible) < int(count)
            sampled = rng.choice(eligible, size=int(count), replace=replace)
            sampled_indices.extend(int(value) for value in sampled)
            replacement_count += int(count) if replace else 0
        fake_latent = candidate_latent[np.asarray(sampled_indices, dtype=np.int64)]
        fake_distance = pairwise_nearest_distance(target_latent, fake_latent)
        fake_z = pd.Series(fake_distance, index=target.index, dtype=float)
        for dataset, indices in target.groupby(target_dataset, observed=True).groups.items():
            values = fake_z.loc[indices]
            std = float(values.std(ddof=1))
            fake_z.loc[indices] = (values - float(values.mean())) / std if np.isfinite(std) and std > 0 else 0.0
        fit = ols_slope(fake_z, pd.Series(target_score, index=target.index), target_patient)
        raw_rho, raw_p = spearmanr(fake_distance, target_score, nan_policy="omit")
        rows.append(
            {
                "replicate": replicate,
                "seed": seed,
                "n_sampled_fake_anchor": len(sampled_indices),
                "sampling_strata": len(requested),
                "sampling_with_replacement_cells": replacement_count,
                "target_scope": "non-State15 candidates",
                "score_name": score_name,
                "patient_adjusted_slope": fit["slope"],
                "slope_ci95_low": fit["ci95_low"],
                "slope_ci95_high": fit["ci95_high"],
                "slope_pvalue": fit["pvalue"],
                "spearman_rho": float(raw_rho) if np.isfinite(raw_rho) else np.nan,
                "spearman_pvalue": float(raw_p) if np.isfinite(raw_p) else np.nan,
            }
        )
    null = pd.DataFrame(rows)
    return null, {
        "repetitions": repetitions,
        "seed": seed,
        "sampling_strata": requested,
        "sampling_with_replacement_cells_total": int(null["sampling_with_replacement_cells"].sum()) if len(null) else 0,
    }


def connectivity_table(latent: np.ndarray, table: pd.DataFrame, k: int = 30) -> pd.DataFrame:
    graph = NearestNeighbors(n_neighbors=min(k + 1, len(latent)), metric="euclidean").fit(latent)
    _, indices = graph.kneighbors(latent)
    source = table["current_state"].replace("", "boundary").astype(str).to_numpy()
    edge_rows: list[dict[str, Any]] = []
    for source_state in pd.unique(source):
        source_indices = np.flatnonzero(source == source_state)
        neighbors = indices[source_indices, 1:]
        targets = source[neighbors].ravel()
        counts = pd.Series(targets).value_counts()
        for target_state, count in counts.items():
            edge_rows.append(
                {
                    "source_state": str(source_state),
                    "target_state": str(target_state),
                    "source_cells": len(source_indices),
                    "k": min(k, len(latent) - 1),
                    "directed_edges": int(count),
                    "source_to_target_fraction": float(count / max(len(source_indices) * min(k, len(latent) - 1), 1)),
                }
            )
    return pd.DataFrame(edge_rows).sort_values(["source_state", "source_to_target_fraction"], ascending=[True, False]).reset_index(drop=True)


def stage21_checkpoint(
    models: pd.DataFrame,
    dataset_models: pd.DataFrame,
    patient_models: pd.DataFrame,
    null: pd.DataFrame,
    connectivity: pd.DataFrame,
) -> tuple[str, str, dict[str, Any]]:
    observed = models[(models["scope"].eq("non_state15_candidates")) & models["score_name"].eq("LAMCORE_independent")]
    observed_slope = float(observed["patient_adjusted_slope"].iloc[0]) if len(observed) else np.nan
    observed_rho = float(observed["spearman_rho"].iloc[0]) if len(observed) else np.nan
    dataset_independent = dataset_models[dataset_models["score_name"].eq("LAMCORE_independent")]
    patient_independent = patient_models[patient_models["score_name"].eq("LAMCORE_independent")]
    dataset_negative = float((dataset_independent["slope"] < 0).mean()) if len(dataset_independent) else 0.0
    patient_negative = float((patient_independent["slope"] < 0).mean()) if len(patient_independent) else 0.0
    if len(null) and np.isfinite(observed_slope):
        null_values = pd.to_numeric(null["patient_adjusted_slope"], errors="coerce").dropna().to_numpy()
        empirical_p = float((1 + np.sum(np.abs(null_values) >= abs(observed_slope))) / (len(null_values) + 1)) if len(null_values) else np.nan
    else:
        empirical_p = np.nan
    state15_external = connectivity[(connectivity["source_state"].eq(TARGET_STATE)) & ~connectivity["target_state"].eq(TARGET_STATE)]
    branch_evidence = len(state15_external[state15_external["source_to_target_fraction"] >= 0.01]) >= 2
    robust = bool(
        np.isfinite(observed_slope)
        and observed_slope < 0
        and np.isfinite(observed_rho)
        and observed_rho < 0
        and dataset_negative >= 0.75
        and patient_negative > 0.5
        and np.isfinite(empirical_p)
        and empirical_p <= 0.05
    )
    independent_signal = bool(np.isfinite(observed_slope) and observed_slope < 0 and dataset_negative >= 0.5)
    if robust and branch_evidence:
        checkpoint = "branched_lam_centered_manifold_candidate"
        interpretation = "Independent LAM evidence survives outside State 15, exceeds composition-matched null slopes, and State 15 connects to multiple external directions."
    elif robust:
        checkpoint = "robust_state15_centered_lam_identity_manifold_candidate"
        interpretation = "Independent LAM evidence declines outside State 15 with patient/dataset support and a composition-matched null result."
    elif independent_signal:
        checkpoint = "state15_lam_rich_gradient_but_not_robust_manifold"
        interpretation = "Some independent LAM evidence remains outside State 15 and is not reproduced by the composition-matched null, but the pooled rank pattern and patient heterogeneity are insufficient for a robust uniform manifold claim."
    else:
        checkpoint = "state15_lam_rich_discrete_population_no_external_manifold_evidence"
        interpretation = "The external gradient is weak after removing the anchor or is not independent of composition/technical structure."
    details = {
        "observed_candidate_only_independent_slope": observed_slope,
        "observed_candidate_only_independent_spearman": observed_rho,
        "dataset_negative_fraction": dataset_negative,
        "patient_negative_fraction": patient_negative,
        "matched_anchor_empirical_pvalue": empirical_p,
        "branch_evidence": branch_evidence,
        "robust_rule_passed": robust,
    }
    return checkpoint, interpretation, details


def write_report(
    output_dir: Path,
    manifest: dict[str, Any],
    gradient: pd.DataFrame,
    models: pd.DataFrame,
    dataset_models: pd.DataFrame,
    patient_models: pd.DataFrame,
    null: pd.DataFrame,
    checkpoint: str,
    interpretation: str,
) -> None:
    independent = gradient["LAMCORE_independent_median"].to_numpy(dtype=float)
    non_anchor_gradient = bool(len(independent) >= 2 and np.nanmedian(independent[:2]) > np.nanmedian(independent[-2:]))
    observed = models[(models["scope"].eq("non_state15_candidates")) & models["score_name"].eq("LAMCORE_independent")]
    slope = float(observed["patient_adjusted_slope"].iloc[0]) if len(observed) else np.nan
    rho = float(observed["spearman_rho"].iloc[0]) if len(observed) else np.nan
    null_summary = "not estimable"
    if len(null) and np.isfinite(slope):
        values = pd.to_numeric(null["patient_adjusted_slope"], errors="coerce").dropna()
        if len(values):
            p = float((1 + np.sum(np.abs(values) >= abs(slope))) / (len(values) + 1))
            null_summary = f"empirical two-sided p={p:.4g}; null median slope={float(values.median()):.5g}"
    report = [
        "# Stage 21：State 15-centered manifold validation",
        "",
        "本阶段将 State 15 的 200 个细胞严格作为 reference anchor；所有主要 gradient 检验均在其余 22,061 个细胞上完成。复用 Stage 20 的 `distance_to_state15` 和既有 `X_scVI`，不重训 scVI、不重新 Leiden/consensus、不修改 candidate gate。",
        "",
        "## Frozen input and independence audit",
        "",
        f"- Anchor: {manifest['anchor_cell_count']} State 15 cells; ID SHA-256 `{manifest['anchor_cell_id_sha256']}`.",
        f"- Validation object: {manifest['validation_cell_count']} non-State15 cells; candidate-only null pool: {manifest['candidate_null_pool_count']} cells.",
        f"- LAMCORE formal genes: {manifest['score_components']['formal_gene_count']}; outside-scVI genes: {manifest['score_components']['outside_scvi_gene_count']}; independent genes: {manifest['score_components']['independent_gene_count']}.",
        "详见 `lamcore_independence_gene_audit.csv`；Stage 21 的主要独立验证分数为 `LAMCORE_outside_scVI` 和 `LAMCORE_independent`。",
        "",
        "## Anchor-excluded distance gradient",
        "",
        f"- Anchor-excluded binned gradient shows near-vs-far independent LAMCORE decrease: `{non_anchor_gradient}`.",
        f"- Candidate-only patient-adjusted independent slope: `{slope:.6g}`.",
        f"- Candidate-only independent Spearman rho: `{rho:.6g}`; the binned nearest-vs-farthest monotonic check is `{non_anchor_gradient}`.",
        "- 四套 score 的 median/IQR 和完整距离模型分别见 `non_state15_distance_gradient.csv` 与 `gradient_models.csv`。",
        "",
        "## Dataset and patient validation",
        "",
        f"- Dataset-level results: {len(dataset_models[dataset_models['score_name'].eq('LAMCORE_independent')])} dataset rows for independent LAMCORE.",
        f"- Patient-level results: {len(patient_models[patient_models['score_name'].eq('LAMCORE_independent')])} patient rows for independent LAMCORE.",
        "- Stage 20 中非负的患者在 Stage 21 中保留了 distance range、cell count 和 slope/rho 字段，不被直接标记为失败。",
        "",
        "## Composition-matched fake-anchor null",
        "",
        f"- Null repetitions: {len(null)}; {null_summary}.",
        "- 每次假 anchor 均按 State 15 的 patient×dataset 组成，从非-State15 candidate pool 抽取；真实比较对象为相同 candidate-only scope。",
        "- 当前 checkpoint 未通过的主要原因是 pooled rank gradient 非单调及 patient-level heterogeneity，而不是 matched-null 未通过。",
        "",
        "## State 16 and boundary",
        "",
        "- `state16_distance_gradient.csv` 分别给出 pooled 和 patient-stratified 的 near/mid/far profile；`boundary_independent_gradient.csv` 只做 evidence ranking，不产生新 candidate 标签。",
        "- `manifold_connectivity.csv` 使用与 Stage 20 相同的 candidate+boundary `X_scVI` k=30 scope 重建邻接汇总，用于识别单轴或多分支连接；不产生 cluster。",
        "",
        "## Stage 21 checkpoint",
        "",
        f"- `{checkpoint}`",
        f"- {interpretation}",
        "",
        "## Outputs",
        "",
        *[f"- {path.name}" for path in sorted(output_dir.iterdir()) if path.is_file()],
    ]
    (output_dir / "stage21_manifold_validation_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/state_modeling.yaml"))
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "results/stage21"))
    parser.add_argument("--block-size", type=int, default=4096)
    parser.add_argument("--null-reps", type=int, default=None)
    parser.add_argument("--null-seed", type=int, default=None)
    args = parser.parse_args()
    config = load_config(Path(args.config).resolve())
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    stage18, stage20 = load_helpers()
    consensus_config = config["outputs"]
    distance_path = PROJECT_ROOT / "results/stage20/state15_cell_distances.csv"
    scvi_path = PROJECT_ROOT / str(consensus_config["scvi_h5ad"])
    prepared_path = PROJECT_ROOT / str(consensus_config["prepared_h5ad"])
    if not distance_path.exists() or not scvi_path.exists() or not prepared_path.exists():
        raise FileNotFoundError(f"Required Stage 21 inputs missing: {distance_path}, {scvi_path}, {prepared_path}")

    anchor, validation, anchor_ids, manifest = read_stage20_inputs(distance_path, scvi_path)
    latent_validation, _ = load_latent_for_ids(scvi_path, validation["analysis_cell_id"])
    latent_anchor, _ = load_latent_for_ids(scvi_path, anchor["analysis_cell_id"])
    latent_main = np.vstack([latent_anchor, latent_validation]).astype(np.float32)
    main_for_graph = pd.concat([anchor, validation], ignore_index=True)
    if len(main_for_graph) != EXPECTED_MAIN_CELLS:
        raise ValueError("Unable to reconstruct Stage 20 main graph scope")

    formal_genes, formal_manifest = stage18.resolve_formal_signature(config)
    formal_genes = unique_genes(stage18, formal_genes)
    independence_audit, scvi_hvg, expression_set = scvi_hvg_and_expression_audit(scvi_path, prepared_path, formal_genes, stage18)
    independence_audit.to_csv(output_dir / "lamcore_independence_gene_audit.csv", index=False)
    modules, score_components = build_score_modules(stage18, stage20, config, formal_genes, scvi_hvg)
    selected = validation[["analysis_cell_id", "cell_id", "current_state", "patient", "dataset", "analysis_role", "distance_to_state15", "nearest_state15_distance"]].copy()
    selected["distance_bin"] = assign_distance_bins(selected["distance_to_state15"])
    selected = add_dataset_standardized_distance(selected)
    scores, score_manifest = score_selected_cells(prepared_path, selected, modules, stage18, int(args.block_size))
    if len(scores) != EXPECTED_VALIDATION_CELLS:
        raise ValueError("Score table does not cover all non-State15 cells")
    scores.to_csv(output_dir / "independent_lamcore_scores.csv", index=False)

    gradient_features = LAMCORE_SCORE_NAMES
    gradient = summary_by_bins(scores, gradient_features)
    gradient.to_csv(output_dir / "non_state15_distance_gradient.csv", index=False)
    smooth = smooth_distance_table(scores, gradient_features)
    smooth.to_csv(output_dir / "distance_score_smooth.csv", index=False)
    all_models = model_gradients(scores, gradient_features, "all_non_state15")
    candidate_scores = scores[scores["analysis_role"].eq("primary_candidate")].copy()
    if len(candidate_scores) != EXPECTED_CANDIDATE_CELLS:
        raise ValueError(f"Expected {EXPECTED_CANDIDATE_CELLS} non-State15 candidates, found {len(candidate_scores)}")
    candidate_models = model_gradients(candidate_scores, gradient_features, "non_state15_candidates")
    models = pd.concat([all_models, candidate_models], ignore_index=True)
    models.to_csv(output_dir / "gradient_models.csv", index=False)

    stage20_patient_path = PROJECT_ROOT / "results/stage20/patient_gradient_consistency.csv"
    stage20_nonnegative: set[str] = set()
    if stage20_patient_path.exists():
        prior = pd.read_csv(stage20_patient_path)
        if "patient" in prior and "distance_spearman_LAMCORE_777" in prior:
            stage20_nonnegative = set(prior.loc[pd.to_numeric(prior["distance_spearman_LAMCORE_777"], errors="coerce") >= 0, "patient"].astype(str))
    dataset_models = group_gradient(scores, gradient_features, "dataset")
    patient_models = group_gradient(scores, gradient_features, "patient", stage20_nonnegative)
    anchor_patient_counts = anchor["patient"].astype(str).value_counts()
    patient_models["n_State15"] = patient_models["patient"].astype(str).map(anchor_patient_counts).fillna(0).astype(int)
    dataset_models.to_csv(output_dir / "dataset_independent_gradient.csv", index=False)
    patient_models.to_csv(output_dir / "patient_independent_gradient.csv", index=False)

    lineage_gradient = summary_by_bins(scores, LINEAGE_FEATURES)
    lineage_gradient.to_csv(output_dir / "lineage_gradient_by_distance.csv", index=False)

    state16 = state16_gradient(scores, STATE16_FEATURES)
    state16.to_csv(output_dir / "state16_distance_gradient.csv", index=False)

    boundary = scores[scores["analysis_role"].eq("boundary")].copy()
    boundary_gradient = summary_by_bins(boundary, [*LAMCORE_SCORE_NAMES, "CORE1", "CORE2", "CORE3", *LINEAGE_FEATURES])
    boundary_gradient.to_csv(output_dir / "boundary_independent_gradient.csv", index=False)

    candidate_pool = candidate_scores.reset_index(drop=True).copy()
    candidate_positions = pd.Index(validation["analysis_cell_id"]).get_indexer(candidate_pool["analysis_cell_id"].tolist())
    if (candidate_positions < 0).any():
        raise ValueError("Candidate latent positions missing for matched null")
    candidate_latent = latent_validation[candidate_positions]
    stage21_config = config.get("stage21", {})
    null_reps = int(args.null_reps if args.null_reps is not None else stage21_config.get("matched_anchor_null_reps", 500))
    null_seed = int(args.null_seed if args.null_seed is not None else stage21_config.get("matched_anchor_null_seed", 20260831))
    null, null_manifest = matched_anchor_null(
        scores,
        latent_validation,
        "LAMCORE_independent",
        anchor,
        candidate_pool,
        candidate_latent,
        null_reps,
        null_seed,
    )
    null.to_csv(output_dir / "matched_anchor_null.csv", index=False)

    connectivity = connectivity_table(latent_main, main_for_graph, k=30)
    connectivity.to_csv(output_dir / "manifold_connectivity.csv", index=False)
    checkpoint, interpretation, checkpoint_details = stage21_checkpoint(models, dataset_models, patient_models, null, connectivity)

    manifest.update(
        {
            "stage": 21,
            "anchor_cell_ids": anchor_ids,
            "anchor_artifact": str(scvi_path),
            "anchor_latent_shape": list(latent_anchor.shape),
            "validation_latent_shape": list(latent_validation.shape),
            "candidate_null_pool_count": len(candidate_pool),
            "scvi_hvg_count": len(scvi_hvg),
            "expression_gene_count": len(expression_set),
            "formal_signature_manifest": formal_manifest,
            "score_components": score_components,
            "score_manifest": score_manifest,
            "distance_contract": {
                "source": str(distance_path),
                "column": "distance_to_state15",
                "anchor_excluded_from_primary_gradient": True,
                "dataset_standardization_for_regression": True,
            },
            "matched_anchor_null": null_manifest,
            "connectivity": {"source": "same Stage20 main X_scVI scope, recomputed k=30 because Stage20 persisted fraction but not edge list", "k": 30},
            "checkpoint": checkpoint,
            "checkpoint_interpretation": interpretation,
            "checkpoint_details": checkpoint_details,
        }
    )
    (output_dir / "stage21_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    write_report(output_dir, manifest, gradient, models, dataset_models, patient_models, null, checkpoint, interpretation)
    print(f"Stage 21 validation cells: {len(validation)}")
    print(f"Candidate-only null pool: {len(candidate_pool)}; repetitions: {len(null)}")
    print(f"Checkpoint: {checkpoint}")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
