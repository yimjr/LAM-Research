#!/usr/bin/env python3
"""Audit cross-dataset identity calibration without changing the Stage 16 gate.

This script reads Stage 16 tables, extracts raw counts only for upstream
positive-reference cells, and decomposes why references were recovered or
missed. It never writes a Stage 16 artifact, retrains scVI, or reclusters.

If the formal LAMCORE signature was added after Stage 16 was run, its score is
reported as supplemental only; it is not substituted into the historical
Stage 16 score or assignment.
"""

from __future__ import annotations

import argparse
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


PROJECT_ROOT = Path(__file__).resolve().parent
DATASETS = ["GSE135851", "GSE190260", "GSE217108", "GSE302356"]
FOCUS_DATASET = "GSE190260"
ALIASES = {"FIGF": "VEGFD"}
MARKER_GENES = ["PMEL", "MLANA", "MITF", "ACTA2", "ESR1", "VEGFD", "CTSK"]


def load_stage16_module() -> Any:
    path = PROJECT_ROOT / "16_rebuild_lam_identity_gate.py"
    spec = importlib.util.spec_from_file_location("stage16_identity_gate", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Stage 16 module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def bool_series(values: pd.Series) -> pd.Series:
    return values.map(as_bool).astype(bool)


def canonical_gene(gene: str) -> str:
    upper = str(gene).strip().upper()
    return ALIASES.get(upper, upper)


def unique_genes(genes: list[str]) -> list[str]:
    result: list[str] = []
    for gene in genes:
        canonical = canonical_gene(gene)
        if canonical and canonical not in result:
            result.append(canonical)
    return result


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def configured_roots(config: dict[str, Any]) -> list[Path]:
    roots = [PROJECT_ROOT / str(root) for root in config.get("input_roots", [])]
    roots.extend(
        [
            Path("/mnt/e/LAM-Research/LAM-Cell-Research"),
            Path("/mnt/e/LAM-Research/data-temp"),
            PROJECT_ROOT / "data/upstream",
        ]
    )
    result: list[Path] = []
    seen: set[str] = set()
    for root in roots:
        key = str(root.resolve())
        if key not in seen:
            result.append(root)
            seen.add(key)
    return result


def resolve_prepared_path(config: dict[str, Any]) -> Path:
    path = PROJECT_ROOT / str(config["outputs"]["prepared_h5ad"])
    if not path.exists():
        raise FileNotFoundError(f"Prepared AnnData not found: {path}")
    return path


def resolve_formal_signature(config: dict[str, Any]) -> tuple[list[str], dict[str, Any]]:
    candidates: list[Path] = []
    for root in configured_roots(config):
        candidates.extend(
            [
                root / "data/raw/reference/LAM_core_signature_genes.csv",
                root / "LAM_core_signature_genes.csv",
            ]
        )
    for path in candidates:
        if path.exists():
            table = pd.read_csv(path)
            if table.empty:
                return [], {"status": "empty", "path": str(path)}
            table.columns = [str(column).strip() for column in table.columns]
            gene_column = next(
                (column for column in ["Gene", "gene", "gene_symbol", "symbol"] if column in table.columns),
                str(table.columns[0]),
            )
            genes = unique_genes(table[gene_column].dropna().astype(str).tolist())
            return genes, {
                "status": "available",
                "path": str(path),
                "gene_column": gene_column,
                "n_genes": len(genes),
            }
    return [], {"status": "not_available", "path": ""}


def read_stage16_inputs(config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    stage16_dir = PROJECT_ROOT / "results/stage16"
    evidence_path = stage16_dir / "cell_identity_evidence.csv"
    calibration_path = stage16_dir / "reference_calibration.csv"
    lodo_path = stage16_dir / "leave_one_dataset_out_validation.csv"
    manifest_path = stage16_dir / "stage16_manifest.json"
    required = [evidence_path, calibration_path, lodo_path]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing Stage 16 input(s): " + ", ".join(missing))
    evidence = pd.read_csv(evidence_path)
    calibration = pd.read_csv(calibration_path)
    lodo = pd.read_csv(lodo_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8")) if manifest_path.exists() else {}
    required_evidence = {
        "analysis_cell_id", "cell_id", "dataset", "patient_id", "positive_reference",
        "positive_reference_type", "LAM_identity_score", "identity_anchor_score", "support_score",
        "competing_lineage_penalty", "LAM_core_candidate", "LAM_boundary_candidate", "non_LAM_like",
    }
    missing_columns = sorted(required_evidence - set(evidence.columns))
    if missing_columns:
        raise ValueError("Stage 16 evidence is missing columns: " + ", ".join(missing_columns))
    evidence["dataset"] = evidence["dataset"].astype(str)
    return evidence, calibration, lodo, manifest


def dense(matrix: Any) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def extract_raw_and_log_scores(
    prepared: ad.AnnData,
    positive: pd.DataFrame,
    modules: dict[str, list[str]],
    block_size: int = 2048,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Extract counts for positive references and calculate log-normalized scores."""

    if "counts" not in prepared.layers:
        raise RuntimeError("Prepared AnnData has no layers['counts']; raw-count audit cannot run")
    obs_names = pd.Index(prepared.obs_names.astype(str))
    requested_ids = positive["analysis_cell_id"].astype(str).tolist()
    positions = obs_names.get_indexer(requested_ids)
    missing = [requested_ids[index] for index, value in enumerate(positions) if value < 0]
    if missing:
        raise ValueError(f"Positive-reference IDs missing from prepared AnnData: {len(missing)}")

    var_names = pd.Index(prepared.var_names.astype(str))
    actual_by_canonical: dict[str, list[str]] = {}
    for actual in var_names:
        actual_by_canonical.setdefault(canonical_gene(actual), []).append(str(actual))
    canonical_union: list[str] = []
    for genes in modules.values():
        for gene in unique_genes(genes):
            if gene in actual_by_canonical and gene not in canonical_union:
                canonical_union.append(gene)
    actual_union: list[str] = []
    canonical_to_actual_indices: dict[str, list[int]] = {}
    for gene in canonical_union:
        canonical_to_actual_indices[gene] = []
        for actual in actual_by_canonical[gene]:
            canonical_to_actual_indices[gene].append(len(actual_union))
            actual_union.append(actual)

    n_cells = len(positive)
    counts_by_gene = {gene: np.zeros(n_cells, dtype=np.float32) for gene in canonical_union}
    library_size = pd.to_numeric(positive.get("total_counts", np.nan), errors="coerce").to_numpy(dtype=np.float32)
    for start in range(0, n_cells, block_size):
        stop = min(start + block_size, n_cells)
        row_positions = np.asarray(positions[start:stop], dtype=np.int64)
        subset = prepared[row_positions, actual_union].to_memory()
        count_values = dense(subset.layers["counts"])
        if not np.isfinite(library_size[start:stop]).all():
            inferred = count_values.sum(axis=1)
            library_size[start:stop] = np.where(library_size[start:stop] > 0, library_size[start:stop], inferred)
        for gene, indices in canonical_to_actual_indices.items():
            counts_by_gene[gene][start:stop] = count_values[:, indices].sum(axis=1)

    safe_library = np.maximum(library_size, 1.0)
    log_by_gene = {
        gene: np.log1p(values / safe_library * 10000.0).astype(np.float32)
        for gene, values in counts_by_gene.items()
    }
    scores: dict[str, np.ndarray] = {}
    detections: dict[str, np.ndarray] = {}
    for gene, values in log_by_gene.items():
        scores[gene] = values
        detections[gene] = counts_by_gene[gene] > 0
    for name, genes in modules.items():
        present = [gene for gene in unique_genes(genes) if gene in log_by_gene]
        if present:
            scores[f"{name}_score_recomputed"] = np.column_stack([log_by_gene[gene] for gene in present]).mean(axis=1).astype(np.float32)
            detections[f"{name}_detected"] = np.column_stack([detections[gene] for gene in present]).any(axis=1)
        else:
            scores[f"{name}_score_recomputed"] = np.full(n_cells, np.nan, dtype=np.float32)
            detections[f"{name}_detected"] = np.zeros(n_cells, dtype=bool)

    output = positive[["analysis_cell_id"]].copy().reset_index(drop=True)
    for gene in MARKER_GENES:
        output[gene] = scores.get(gene, np.full(n_cells, np.nan, dtype=np.float32))
        output[f"{gene}_raw_count"] = counts_by_gene.get(gene, np.zeros(n_cells, dtype=np.float32))
        output[f"{gene}_detected"] = detections.get(gene, np.zeros(n_cells, dtype=bool))
    for name, values in scores.items():
        if name.endswith("_score_recomputed"):
            output[name] = values
    for name, values in detections.items():
        if name.endswith("_detected") and name not in output:
            output[name] = values
    output["total_umi_raw_audit"] = library_size
    output["identity_markers_detected"] = output[[f"{gene}_detected" for gene in ["PMEL", "MLANA", "MITF"]]].sum(axis=1).astype(int)
    output["support_markers_detected"] = output[[f"{gene}_detected" for gene in ["ACTA2", "ESR1", "VEGFD", "CTSK"]]].sum(axis=1).astype(int)
    return output, {
        "n_positive_cells": n_cells,
        "n_resolved_canonical_genes": len(canonical_union),
        "resolved_canonical_genes": canonical_union,
        "alias_map": ALIASES,
        "expression_score_source": "layers['counts'] -> library-size normalization target_sum=10000 -> log1p",
        "block_size": block_size,
    }


def apply_manifest_z(values: pd.Series, datasets: pd.Series, parameters: dict[str, Any], module_name: str) -> np.ndarray:
    result = np.full(len(values), np.nan, dtype=np.float32)
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    dataset_values = datasets.astype(str).to_numpy()
    for dataset in sorted(set(dataset_values)):
        mask = dataset_values == dataset
        params = parameters.get(module_name, {}).get(dataset, {})
        median = float(params.get("median", 0.0))
        scale = float(params.get("scale", 1.0))
        scale = scale if np.isfinite(scale) and scale > 1e-8 else 1.0
        result[mask] = ((numeric[mask] - median) / scale).astype(np.float32)
    return result


def get_full_thresholds(manifest: dict[str, Any]) -> dict[str, Any]:
    thresholds = manifest.get("full_thresholds", {})
    required = [
        "LAM_identity_score_core_threshold", "identity_anchor_core_threshold", "support_core_threshold",
        "LAM_identity_score_boundary_threshold", "identity_anchor_boundary_threshold", "support_boundary_threshold",
        "competing_lineage_95pct_cutoffs", "pericyte_VSMC_95pct_cutoff",
    ]
    missing = [name for name in required if name not in thresholds]
    if missing:
        raise ValueError("Stage 16 manifest lacks full thresholds: " + ", ".join(missing))
    return thresholds


def classify_failure(row: pd.Series, thresholds: dict[str, Any]) -> tuple[str, str, bool, bool]:
    missed = not as_bool(row["LAM_core_candidate"]) and not as_bool(row["LAM_boundary_candidate"])
    if not missed:
        return "not_missed", "not_missed", False, False
    identity_fail = (
        float(row["identity_anchor_score"]) < float(thresholds["identity_anchor_boundary_threshold"])
        or float(row["support_score"]) < float(thresholds["support_boundary_threshold"])
        or float(row["LAM_identity_score"]) < float(thresholds["LAM_identity_score_boundary_threshold"])
    )
    competing_fail = as_bool(row["competing_lineage_exclusion"])
    composite_high = float(row["LAM_identity_score"]) >= float(thresholds["LAM_identity_score_core_threshold"])
    if identity_fail and competing_fail:
        detailed = "both"
        root = "both"
    elif competing_fail:
        detailed = "excessive_competing_penalty"
        root = "competing_penalty_only"
    elif composite_high:
        detailed = "score_above_expected_but_threshold_failed"
        root = "calibration_threshold_only"
    else:
        detailed = "insufficient_positive_evidence"
        root = "low_identity_evidence_only"
    return detailed, root, identity_fail, competing_fail


def best_threshold(positive: np.ndarray, negative: np.ndarray) -> tuple[float, dict[str, float]]:
    positive = positive[np.isfinite(positive)]
    negative = negative[np.isfinite(negative)]
    if positive.size == 0 or negative.size == 0:
        fallback = float(np.nanmedian(positive)) if positive.size else 0.0
        return fallback, {"balanced_accuracy": np.nan, "tpr": np.nan, "tnr": np.nan, "youden_j": np.nan}
    candidates = np.unique(np.concatenate([positive, negative]).astype(float))
    best: tuple[float, float, float, float, float] | None = None
    for threshold in candidates:
        tpr = float(np.mean(positive >= threshold))
        tnr = float(np.mean(negative < threshold))
        youden = tpr + tnr - 1.0
        balanced = 0.5 * (tpr + tnr)
        candidate = (youden, tnr, threshold, tpr, balanced)
        if best is None or candidate > best:
            best = candidate
    assert best is not None
    return float(best[2]), {"balanced_accuracy": best[4], "tpr": best[3], "tnr": best[1], "youden_j": best[0]}


def q10(values: np.ndarray, fallback: float) -> float:
    values = values[np.isfinite(values)]
    return float(np.quantile(values, 0.1)) if values.size else fallback


def lodo_recovery(
    evidence: pd.DataFrame,
    calibration: pd.DataFrame,
    stage16_lodo: pd.DataFrame,
) -> pd.DataFrame:
    positive = evidence[bool_series(evidence["positive_reference"])].copy()
    positive["dataset"] = positive["dataset"].astype(str)
    rows: list[dict[str, Any]] = []
    for held_out in DATASETS:
        label = f"leave_out_{held_out}"
        sub_cal = calibration[calibration["calibration"].astype(str).eq(label)].copy()
        threshold_by_metric = dict(zip(sub_cal["metric"].astype(str), pd.to_numeric(sub_cal["threshold"], errors="coerce")))
        held = positive[positive["dataset"].eq(held_out)].copy()
        train = positive[~positive["dataset"].eq(held_out)].copy()
        n_positive = len(held)
        if n_positive == 0:
            rows.append({
                "held_out_dataset": held_out, "n_positive_reference": 0, "core_recovered": 0,
                "boundary_recovered": 0, "missed": 0, "core_recovery_rate": np.nan,
                "core_or_boundary_recovery_rate": np.nan,
                "full_stage16_core_recovered": np.nan,
                "full_stage16_boundary_or_core_recovered": np.nan,
                "full_stage16_missed": np.nan,
            })
            continue
        core_score = float(threshold_by_metric.get("LAM_identity_score", 0.0))
        core_identity = float(threshold_by_metric.get("identity_anchor_score", 0.0))
        core_support = float(threshold_by_metric.get("support_score", 0.0))
        boundary_score = min(core_score, q10(pd.to_numeric(train["LAM_identity_score"], errors="coerce").to_numpy(dtype=float), core_score))
        boundary_identity = min(core_identity, q10(pd.to_numeric(train["identity_anchor_score"], errors="coerce").to_numpy(dtype=float), core_identity))
        boundary_support = q10(pd.to_numeric(train["support_score"], errors="coerce").to_numpy(dtype=float), core_support)
        non_pericyte_cutoffs = {
            metric.removesuffix("_score_z_negative_q95"): float(value)
            for metric, value in threshold_by_metric.items()
            if metric.endswith("_score_z_negative_q95") and metric != "pericyte_VSMC_score_z_negative_q95"
        }
        pericyte_cutoff = float(threshold_by_metric.get("pericyte_VSMC_score_z_negative_q95", 2.0))
        non_pericyte = np.zeros(len(held), dtype=bool)
        for lineage, cutoff in non_pericyte_cutoffs.items():
            column = f"{lineage}_score_z"
            if column in held:
                non_pericyte |= pd.to_numeric(held[column], errors="coerce").to_numpy(dtype=float) >= cutoff
        pericyte_values = pd.to_numeric(held.get("pericyte_VSMC_score_z", np.nan), errors="coerce").to_numpy(dtype=float)
        identity_values = pd.to_numeric(held["identity_anchor_score"], errors="coerce").to_numpy(dtype=float)
        pericyte_only = (pericyte_values >= pericyte_cutoff) & (identity_values < boundary_identity)
        competing = non_pericyte | pericyte_only
        score_values = pd.to_numeric(held["LAM_identity_score"], errors="coerce").to_numpy(dtype=float)
        support_values = pd.to_numeric(held["support_score"], errors="coerce").to_numpy(dtype=float)
        core = (score_values >= core_score) & (identity_values >= core_identity) & (support_values >= core_support) & ~competing
        boundary = (~core & (score_values >= boundary_score) & (identity_values >= boundary_identity) & (support_values >= boundary_support) & ~competing)
        stage16_row = stage16_lodo[stage16_lodo["held_out_dataset"].astype(str).eq(held_out)]
        reported_core = float(stage16_row["positive_reference_core_recall"].iloc[0]) * n_positive if not stage16_row.empty and pd.notna(stage16_row["positive_reference_core_recall"].iloc[0]) else np.nan
        reported_both = float(stage16_row["positive_reference_boundary_or_core_recall"].iloc[0]) * n_positive if not stage16_row.empty and pd.notna(stage16_row["positive_reference_boundary_or_core_recall"].iloc[0]) else np.nan
        rows.append({
            "held_out_dataset": held_out,
            "n_positive_reference": n_positive,
            "core_recovered": int(core.sum()),
            "boundary_recovered": int((boundary & ~core).sum()),
            "missed": int((~core & ~boundary).sum()),
            "core_recovery_rate": float(core.mean()),
            "core_or_boundary_recovery_rate": float((core | boundary).mean()),
            "full_stage16_core_recovered": reported_core,
            "full_stage16_boundary_or_core_recovered": reported_both,
            "full_stage16_missed": float(n_positive - reported_both) if np.isfinite(reported_both) else np.nan,
            "lodo_core_threshold": core_score,
            "lodo_boundary_threshold": boundary_score,
            "lodo_identity_threshold": core_identity,
            "lodo_support_threshold": core_support,
        })
    return pd.DataFrame(rows)


def counterfactual_calibration(evidence: pd.DataFrame, manifest: dict[str, Any]) -> pd.DataFrame:
    frame = evidence.copy()
    frame["positive_reference_bool"] = bool_series(frame["positive_reference"])
    frame["strong_competing_bool"] = bool_series(frame["strong_competing_reference"]) if "strong_competing_reference" in frame else False
    frame["dataset"] = frame["dataset"].astype(str)
    frame["raw_stage16_score"] = pd.to_numeric(frame["LAM_identity_score"], errors="coerce")
    frame["dataset_z_score"] = frame.groupby("dataset", observed=True)["raw_stage16_score"].transform(
        lambda values: (values - values.median()) / max(float(values.quantile(0.75) - values.quantile(0.25)), 1e-8)
    )
    frame["dataset_percentile_score"] = frame.groupby("dataset", observed=True)["raw_stage16_score"].rank(method="average", pct=True)
    positive = frame[frame["positive_reference_bool"]]
    negative = frame[frame["strong_competing_bool"] & ~frame["positive_reference_bool"]]
    methods = ["raw_stage16_score", "dataset_z_score", "dataset_percentile_score"]
    full_thresholds = get_full_thresholds(manifest)
    rows: list[dict[str, Any]] = []
    for method in methods:
        if method == "raw_stage16_score":
            core_threshold = float(full_thresholds["LAM_identity_score_core_threshold"])
            boundary_threshold = float(full_thresholds["LAM_identity_score_boundary_threshold"])
            calibration_mode = "Stage16 full threshold, score-only counterfactual"
        else:
            core_threshold, _ = best_threshold(positive[method].to_numpy(dtype=float), negative[method].to_numpy(dtype=float))
            boundary_threshold = min(core_threshold, q10(positive[method].to_numpy(dtype=float), core_threshold))
            calibration_mode = "audit score-specific full-reference threshold"
        for dataset in DATASETS:
            pos = positive[positive["dataset"].eq(dataset)]
            neg = negative[negative["dataset"].eq(dataset)]
            core = pos[method] >= core_threshold
            boundary = (pos[method] >= boundary_threshold) & ~core
            neg_non_lam = neg[method] < boundary_threshold
            rows.append({
                "method": method,
                "calibration_mode": calibration_mode,
                "dataset": dataset,
                "n_positive_reference": int(len(pos)),
                "core_threshold": core_threshold,
                "boundary_threshold": boundary_threshold,
                "core_recovered": int(core.sum()),
                "core_or_boundary_recovered": int((core | boundary).sum()),
                "missed": int((~core & ~boundary).sum()),
                "core_recovery_rate": float(core.mean()) if len(pos) else np.nan,
                "core_or_boundary_recovery_rate": float((core | boundary).mean()) if len(pos) else np.nan,
                "n_negative_reference_lam_competing": int(len(neg)),
                "non_LAM_like_retention": float(neg_non_lam.mean()) if len(neg) else np.nan,
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/state_modeling.yaml"))
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "results/stage17"))
    parser.add_argument("--block-size", type=int, default=2048)
    args = parser.parse_args()
    config = load_config(Path(args.config).resolve())
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    evidence, calibration, stage16_lodo, stage16_manifest = read_stage16_inputs(config)
    positive = evidence[bool_series(evidence["positive_reference"])].copy().reset_index(drop=True)
    positive["dataset"] = positive["dataset"].astype(str)
    positive["focus_dataset"] = positive["dataset"].eq(FOCUS_DATASET)
    thresholds = get_full_thresholds(stage16_manifest)

    stage16_module = load_stage16_module()
    lamcore_programs = stage16_manifest.get("lamcore_program_manifest", {}).get("programs", {})
    core2_genes = unique_genes(lamcore_programs.get("CORE2", ["PMEL", "MLANA", "MITF", "TYR", "DCT", "GPNMB", "CTSK", "VEGFD"]))
    core3_genes = unique_genes(lamcore_programs.get("CORE3_identity", ["ACTA2", "PMEL", "MLANA", "MITF", "CTSK", "VEGFD", "ESR1"]))
    formal_genes, formal_manifest = resolve_formal_signature(config)
    modules = {
        "melanocytic": list(stage16_module.MELANOCYTIC_GENES),
        "support": list(stage16_module.SUPPORT_GENES),
        "CORE2": core2_genes,
        "CORE3_identity": core3_genes,
    }
    if formal_genes:
        modules["LAMCORE_777"] = formal_genes
    modules.update({name: list(genes) for name, genes in stage16_module.COMPETING_GENES.items()})

    prepared_path = resolve_prepared_path(config)
    prepared = ad.read_h5ad(prepared_path, backed="r")
    raw_scores, raw_manifest = extract_raw_and_log_scores(prepared, positive, modules, block_size=int(args.block_size))
    prepared.file.close()
    component = positive.merge(raw_scores, on="analysis_cell_id", how="left", validate="one_to_one")

    z_parameters = stage16_manifest.get("robust_z_parameters", {})
    component["CORE2_identity_score_z_recomputed"] = apply_manifest_z(component["CORE2_score_recomputed"], component["dataset"], z_parameters, "lamcore_CORE2")
    component["CORE3_identity_score_z_recomputed"] = apply_manifest_z(component["CORE3_identity_score_recomputed"], component["dataset"], z_parameters, "CORE3_identity")
    component["support_score_z_recomputed"] = apply_manifest_z(component["support_score_recomputed"], component["dataset"], z_parameters, "support")
    component["CORE3_related_score"] = pd.to_numeric(component["lamcore_identity_score"], errors="coerce")
    component["CORE3_related_score_recomputed_from_components"] = component[["CORE2_identity_score_z_recomputed", "CORE3_identity_score_z_recomputed"]].mean(axis=1, skipna=True)
    if "core3_identity_upstream_z" in component:
        component["CORE3_related_score_recomputed_from_components"] = component[["CORE2_identity_score_z_recomputed", "CORE3_identity_score_z_recomputed", "core3_identity_upstream_z"]].mean(axis=1, skipna=True)
    component["positive_component_total"] = 0.60 * pd.to_numeric(component["identity_anchor_score"], errors="coerce") + 0.25 * pd.to_numeric(component["support_score"], errors="coerce")
    component["weighted_competing_penalty"] = 0.20 * pd.to_numeric(component["non_pericyte_competing_score"], errors="coerce") + 0.10 * pd.to_numeric(component["pericyte_VSMC_conditional_penalty"], errors="coerce")
    component["final_LAM_identity_score"] = pd.to_numeric(component["LAM_identity_score"], errors="coerce")
    component["final_score_formula_recomputed"] = component["positive_component_total"] - component["weighted_competing_penalty"]
    component["final_score_formula_abs_diff"] = (component["final_LAM_identity_score"] - component["final_score_formula_recomputed"]).abs()

    failure_details = component.apply(lambda row: classify_failure(row, thresholds), axis=1, result_type="expand")
    failure_details.columns = ["failure_reason", "root_cause_category", "positive_evidence_failed", "competing_penalty_failed"]
    component = pd.concat([component, failure_details], axis=1)
    component["final_assignment"] = np.select(
        [bool_series(component["LAM_core_candidate"]), bool_series(component["LAM_boundary_candidate"])],
        ["LAM_core_candidate", "LAM_boundary_candidate"],
        default="non_LAM_like",
    )
    component["recovered_core"] = component["final_assignment"].eq("LAM_core_candidate")
    component["recovered_boundary"] = component["final_assignment"].eq("LAM_boundary_candidate")
    component["missed"] = component["final_assignment"].eq("non_LAM_like")

    lineage_names = list(stage16_module.COMPETING_GENES)
    lineage_scores = component[[f"{name}_score_z" for name in lineage_names]].apply(pd.to_numeric, errors="coerce")
    component["highest_competing_lineage"] = lineage_scores.idxmax(axis=1).str.removesuffix("_score_z")
    component["highest_competing_score"] = lineage_scores.max(axis=1)
    lineage_array = lineage_scores.to_numpy(dtype=float)
    sorted_indices = np.argsort(lineage_array, axis=1)
    component["second_highest_competing_lineage"] = [lineage_names[index] for index in sorted_indices[:, -2]]
    component["second_highest_competing_score"] = np.sort(lineage_array, axis=1)[:, -2]
    component["focus_dataset"] = component["dataset"].eq(FOCUS_DATASET)

    failure_columns = [
        "analysis_cell_id", "cell_id", "dataset", "focus_dataset", "patient_id", "positive_reference_type",
        "final_assignment", "recovered_core", "recovered_boundary", "missed", "failure_reason", "root_cause_category",
        "positive_evidence_failed", "competing_penalty_failed", *MARKER_GENES,
        *[f"{gene}_raw_count" for gene in MARKER_GENES],
        *[f"{gene}_detected" for gene in MARKER_GENES],
        "identity_markers_detected", "support_markers_detected", "total_umi_raw_audit", "n_genes_by_counts",
        "melanocytic_identity_score", "support_score", "CORE3_related_score",
        "CORE3_related_score_recomputed_from_components", "CORE2_score_recomputed",
        "CORE3_identity_score_recomputed", "LAMCORE_777_score_recomputed", "core3_identity_upstream",
        "core3_identity_upstream_z", "identity_anchor_score", "positive_component_total",
        "non_pericyte_competing_score", "pericyte_VSMC_conditional_penalty", "weighted_competing_penalty",
        "competing_lineage_penalty", "final_LAM_identity_score", "final_score_formula_recomputed",
        "final_score_formula_abs_diff", "highest_competing_lineage", "highest_competing_score",
        "second_highest_competing_lineage", "second_highest_competing_score",
    ]
    failure_columns = list(dict.fromkeys(column for column in failure_columns if column in component.columns))
    component[failure_columns].to_csv(output_dir / "positive_reference_failures.csv", index=False)

    component_score_names = ["melanocytic_identity_score", "support_score", "CORE3_related_score", "competing_lineage_penalty", "final_LAM_identity_score"]
    score_rows: list[dict[str, Any]] = []
    for dataset in DATASETS:
        subset = component[component["dataset"].eq(dataset)]
        for score_name in component_score_names:
            values = pd.to_numeric(subset[score_name], errors="coerce").to_numpy(dtype=float) if score_name in subset else np.array([], dtype=float)
            values = values[np.isfinite(values)]
            score_rows.append({
                "dataset": dataset, "focus_dataset": dataset == FOCUS_DATASET,
                "n_positive_reference": int(len(subset)), "score": score_name, "n": int(len(values)),
                "mean": float(np.mean(values)) if len(values) else np.nan,
                "median": float(np.median(values)) if len(values) else np.nan,
                "q25": float(np.quantile(values, 0.25)) if len(values) else np.nan,
                "q75": float(np.quantile(values, 0.75)) if len(values) else np.nan,
                "detection_rate_score_gt_0": float(np.mean(values > 0)) if len(values) else np.nan,
            })
    pd.DataFrame(score_rows).to_csv(output_dir / "component_scores_by_dataset.csv", index=False)

    marker_rows: list[dict[str, Any]] = []
    for dataset in DATASETS:
        subset = component[component["dataset"].eq(dataset)]
        row: dict[str, Any] = {
            "dataset": dataset, "focus_dataset": dataset == FOCUS_DATASET, "n_positive_reference": int(len(subset)),
            "median_total_umi": float(pd.to_numeric(subset["total_umi_raw_audit"], errors="coerce").median()) if len(subset) else np.nan,
            "median_detected_genes": float(pd.to_numeric(subset["n_genes_by_counts"], errors="coerce").median()) if len(subset) else np.nan,
        }
        for gene in MARKER_GENES:
            row[f"{gene}_detection_rate"] = float(bool_series(subset[f"{gene}_detected"]).mean()) if len(subset) else np.nan
        for label, column in [("identity", "identity_markers_detected"), ("support", "support_markers_detected")]:
            counts = pd.to_numeric(subset[column], errors="coerce")
            for category, mask in [("0", counts.eq(0)), ("1", counts.eq(1)), ("2", counts.eq(2)), ("ge3", counts.ge(3))]:
                row[f"{label}_markers_detected_{category}_fraction"] = float(mask.mean()) if len(subset) else np.nan
        marker_rows.append(row)
    pd.DataFrame(marker_rows).to_csv(output_dir / "marker_detection_by_dataset.csv", index=False)

    penalty_rows: list[dict[str, Any]] = []
    for dataset in DATASETS:
        subset = component[component["dataset"].eq(dataset)]
        for lineage in lineage_names:
            values = pd.to_numeric(subset[f"{lineage}_score_z"], errors="coerce") if len(subset) else pd.Series(dtype=float)
            highest = subset["highest_competing_lineage"].eq(lineage) if len(subset) else pd.Series(dtype=bool)
            penalty_rows.append({
                "dataset": dataset, "focus_dataset": dataset == FOCUS_DATASET, "lineage": lineage,
                "n_positive_reference": int(len(subset)),
                "n_highest_competing": int(highest.sum()) if len(subset) else 0,
                "fraction_highest_competing": float(highest.mean()) if len(subset) else np.nan,
                "n_score_z_ge_2": int((values >= 2).sum()) if len(subset) else 0,
                "fraction_score_z_ge_2": float((values >= 2).mean()) if len(subset) else np.nan,
                "median_score_z": float(values.median()) if len(subset) else np.nan,
                "q95_score_z": float(values.quantile(0.95)) if len(subset) else np.nan,
            })
    pd.DataFrame(penalty_rows).to_csv(output_dir / "competing_lineage_penalty_audit.csv", index=False)

    counterfactual = counterfactual_calibration(evidence, stage16_manifest)
    counterfactual.to_csv(output_dir / "counterfactual_calibration.csv", index=False)
    lodo = lodo_recovery(evidence, calibration, stage16_lodo)
    lodo.to_csv(output_dir / "lodo_recovery.csv", index=False)

    root_rows: list[dict[str, Any]] = []
    lodo_by_dataset = lodo.rename(columns={"held_out_dataset": "dataset"})
    for dataset in DATASETS:
        sub = component[component["dataset"].eq(dataset)]
        n = len(sub)
        category_counts = sub["root_cause_category"].value_counts().to_dict()
        median_score = float(pd.to_numeric(sub["final_LAM_identity_score"], errors="coerce").median()) if n else np.nan
        other = component[~component["dataset"].eq(dataset)]
        other_median = float(pd.to_numeric(other["final_LAM_identity_score"], errors="coerce").median()) if len(other) else np.nan
        largest_category = max(
            ["low_identity_evidence_only", "competing_penalty_only", "both", "calibration_threshold_only"],
            key=lambda name: int(category_counts.get(name, 0)),
            default="none",
        )
        lodo_row = lodo_by_dataset[lodo_by_dataset["dataset"].eq(dataset)]
        lodo_core_both = float(lodo_row["core_or_boundary_recovery_rate"].iloc[0]) if not lodo_row.empty else np.nan
        full_core_both = float((sub["recovered_core"] | sub["recovered_boundary"]).mean()) if n else np.nan
        root_rows.append({
            "dataset": dataset, "focus_dataset": dataset == FOCUS_DATASET, "n_positive_reference": n,
            "core_recovered": int(sub["recovered_core"].sum()) if n else 0,
            "boundary_recovered": int(sub["recovered_boundary"].sum()) if n else 0,
            "missed": int(sub["missed"].sum()) if n else 0,
            "positive_reference_core_recovery": float(sub["recovered_core"].mean()) if n else np.nan,
            "positive_reference_core_or_boundary_recovery": full_core_both,
            "identity_marker_dropout_fraction_lt2": float((pd.to_numeric(sub["identity_markers_detected"], errors="coerce") < 2).mean()) if n else np.nan,
            "support_marker_dropout_fraction_lt2": float((pd.to_numeric(sub["support_markers_detected"], errors="coerce") < 2).mean()) if n else np.nan,
            "median_total_umi": float(pd.to_numeric(sub["total_umi_raw_audit"], errors="coerce").median()) if n else np.nan,
            "median_detected_genes": float(pd.to_numeric(sub["n_genes_by_counts"], errors="coerce").median()) if n else np.nan,
            "median_final_LAM_identity_score": median_score,
            "positive_score_shift_vs_other_positive_median": median_score - other_median if np.isfinite(median_score) and np.isfinite(other_median) else np.nan,
            "median_competing_lineage_penalty": float(pd.to_numeric(sub["competing_lineage_penalty"], errors="coerce").median()) if n else np.nan,
            "low_identity_evidence_only_fraction": float(category_counts.get("low_identity_evidence_only", 0) / n) if n else np.nan,
            "competing_penalty_only_fraction": float(category_counts.get("competing_penalty_only", 0) / n) if n else np.nan,
            "both_fraction": float(category_counts.get("both", 0) / n) if n else np.nan,
            "calibration_threshold_only_fraction": float(category_counts.get("calibration_threshold_only", 0) / n) if n else np.nan,
            "primary_root_cause_category": largest_category if n else "no_positive_reference",
            "lodo_core_or_boundary_recovery": lodo_core_both,
            "lodo_minus_full_recovery_shift": lodo_core_both - full_core_both if np.isfinite(lodo_core_both) and np.isfinite(full_core_both) else np.nan,
            "lodo_core_threshold": float(lodo_row["lodo_core_threshold"].iloc[0]) if not lodo_row.empty and "lodo_core_threshold" in lodo_row else np.nan,
            "full_core_threshold": float(thresholds["LAM_identity_score_core_threshold"]),
            "positive_score_shift_flag": "lower" if np.isfinite(median_score) and np.isfinite(other_median) and median_score < other_median else "not_lower_or_unavailable",
            "marker_dropout_flag": "possible" if n and float((pd.to_numeric(sub["identity_markers_detected"], errors="coerce") < 2).mean()) >= 0.5 else "not_dominant",
            "competing_penalty_flag": "possible" if n and float(pd.to_numeric(sub["competing_lineage_penalty"], errors="coerce").median()) > 0 else "not_dominant",
        })
    root_cause = pd.DataFrame(root_rows)
    root_cause.to_csv(output_dir / "root_cause_by_dataset.csv", index=False)

    gse = root_cause[root_cause["dataset"].eq(FOCUS_DATASET)].iloc[0] if (root_cause["dataset"] == FOCUS_DATASET).any() else None
    report_lines = [
        "# Stage 17 — Cross-dataset identity calibration audit",
        "",
        "This is a read-only audit of Stage 16. It does not change the formal candidate gate, retrain scVI, or recluster cells.",
        "",
        f"Positive-reference cells audited: {len(component):,}",
        f"Focus dataset: {FOCUS_DATASET}",
        f"Formal 777-gene signature at audit time: {formal_manifest.get('status')}",
        "",
        "## GSE190260 root-cause summary",
        "",
    ]
    if gse is not None:
        report_lines.extend(
            [
                f"- Positive references: {int(gse['n_positive_reference']):,}",
                f"- Stage 16 core recovered: {int(gse['core_recovered']):,} ({float(gse['positive_reference_core_recovery']):.3f})",
                f"- Stage 16 core+boundary recovered: {int(gse['core_recovered'] + gse['boundary_recovered']):,} ({float(gse['positive_reference_core_or_boundary_recovery']):.3f})",
                f"- Missed: {int(gse['missed']):,}",
                f"- Primary failure category: {gse['primary_root_cause_category']}",
                f"- Identity-marker dropout (<2 detected): {float(gse['identity_marker_dropout_fraction_lt2']):.3f}",
                f"- Median final score: {float(gse['median_final_LAM_identity_score']):.3f}; shift vs other positive references: {float(gse['positive_score_shift_vs_other_positive_median']):.3f}",
                f"- Median competing-lineage penalty: {float(gse['median_competing_lineage_penalty']):.3f}",
                f"- LODO core+boundary recovery: {float(gse['lodo_core_or_boundary_recovery']):.3f}",
            ]
        )
    report_lines.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "Failure categories are computed from the Stage 16 component thresholds and exclusion flags, not from manual state inspection. The core3_like label remains an upstream reference label; it is not treated as a new formal candidate assignment. A newly available formal signature is supplemental in this audit and is not inserted into the historical Stage 16 score.",
            "",
            "## Outputs",
            "",
            "- positive_reference_failures.csv: cell-level reference decomposition and failure reason.",
            "- component_scores_by_dataset.csv: positive-reference score distributions.",
            "- marker_detection_by_dataset.csv: raw-count dropout and depth audit.",
            "- competing_lineage_penalty_audit.csv: lineage-specific penalty summary.",
            "- counterfactual_calibration.csv: raw, within-dataset z-score and percentile score audit.",
            "- lodo_recovery.csv: explicit held-out positive-reference recovery counts.",
            "- root_cause_by_dataset.csv: dataset-level attribution table.",
        ]
    )
    (output_dir / "identity_calibration_audit.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    audit_manifest = {
        "script": str(Path(__file__).resolve()),
        "read_only_stage16": True,
        "no_scvi_training": True,
        "no_reclustering": True,
        "focus_dataset": FOCUS_DATASET,
        "datasets": DATASETS,
        "stage16_inputs": {
            "cell_identity_evidence": str(PROJECT_ROOT / "results/stage16/cell_identity_evidence.csv"),
            "reference_calibration": str(PROJECT_ROOT / "results/stage16/reference_calibration.csv"),
            "leave_one_dataset_out_validation": str(PROJECT_ROOT / "results/stage16/leave_one_dataset_out_validation.csv"),
        },
        "prepared_h5ad": str(prepared_path),
        "formal_signature_manifest": formal_manifest,
        "raw_score_manifest": raw_manifest,
        "positive_reference_definition": "Stage 16 positive_reference == True; no new reference label created",
        "counterfactual_negative_scope": "Stage 16 cell_identity_evidence strong_competing_reference LAM rows; normal/control rows are not in the fixed cell-level input table",
        "counterfactual_methods": ["raw_stage16_score", "dataset_z_score", "dataset_percentile_score"],
        "outputs": sorted(path.name for path in output_dir.iterdir() if path.is_file()),
    }
    (output_dir / "stage17_manifest.json").write_text(json.dumps(audit_manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    print(f"Audited {len(component):,} Stage 16 positive-reference cells")
    print(f"Focus dataset {FOCUS_DATASET}: {int((component['dataset'].eq(FOCUS_DATASET) & component['missed']).sum()):,} missed positive references")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
