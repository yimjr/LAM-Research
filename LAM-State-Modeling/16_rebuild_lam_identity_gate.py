#!/usr/bin/env python3
"""Rebuild the LAM candidate identity gate without reclustering or scVI training.

The script scores every ``condition == LAM`` cell in the prepared AnnData.
It separates identity anchors (PMEL/MLANA/MITF and CORE2/CORE3-derived
modules) from supportive markers (ACTA2/ESR1/VEGFD/CTSK), and reports
competing lineage evidence as continuous scores.  Thresholds are calibrated
only from independent upstream positive references, the normal reference, and
strong competing-lineage references.  Existing consensus states are joined
only at the final diagnostic step.

No existing AnnData or Stage 1--13 result is written by this script.  New
outputs are written below ``results/stage16``.
"""

from __future__ import annotations

import argparse
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
ALIASES = {"FIGF": "VEGFD"}
MELANOCYTIC_GENES = ["PMEL", "MLANA", "MITF"]
SUPPORT_GENES = ["ACTA2", "ESR1", "VEGFD", "CTSK"]
COMPETING_GENES = {
    "ciliated": [
        "FOXJ1", "PIFO", "TPPP3", "CFAP43", "CFAP44", "CFAP54", "CFAP57",
        "DNAH5", "DNAH9", "DNAH11", "DNAI1", "DNAI2", "HYDIN", "RSPH1",
        "RSPH4A", "RSPH9", "CCDC39", "CCDC40", "C20orf85",
    ],
    "AT2": ["SFTPC", "SFTPA1", "SFTPA2", "NAPSA"],
    "macrophage": [
        "LST1", "TYROBP", "FCER1G", "C1QA", "C1QB", "C1QC", "FABP4",
        "AIF1", "CTSS", "APOE", "TREM2",
    ],
    "endothelial": [
        "PECAM1", "EMCN", "VWF", "CCL21", "FLT1", "KDR", "FLT4", "RAMP2",
        "CA4", "ENG",
    ],
    "fibroblast": [
        "COL1A1", "COL1A2", "DCN", "LUM", "PI16", "COL3A1", "COL6A1",
        "COL6A2", "PDGFRA", "CFD", "C7",
    ],
    "mesothelial": ["ITLN1", "MSLN", "WT1", "PDPN", "MUC16", "UPK3B", "CALB2"],
    # ACTA2/MYH11 are retained here as contextual evidence, never as a sole
    # exclusion rule.  The assignment logic only uses this module for a
    # conditional penalty when identity anchors are weak.
    "pericyte_VSMC": [
        "RGS5", "PDGFRB", "CSPG4", "COX4I2", "MYH11", "ACTA2", "TAGLN",
        "CNN1", "DES", "CALD1", "MYL9", "TPM2",
    ],
}
NON_PERICYTE_LINEAGES = ["ciliated", "AT2", "macrophage", "endothelial", "fibroblast", "mesothelial"]


def as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def bool_array(values: pd.Series) -> np.ndarray:
    return values.map(as_bool).to_numpy(dtype=bool)


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
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


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


def resolve_file(config: dict[str, Any], relative: str, extra: list[Path] | None = None) -> Path | None:
    roots = configured_roots(config) + list(extra or [])
    for root in roots:
        path = root / relative
        if path.exists():
            return path
    return None


def resolve_prepared_path(config: dict[str, Any]) -> Path:
    path = PROJECT_ROOT / str(config["outputs"]["prepared_h5ad"])
    if path.exists():
        return path
    raise FileNotFoundError(f"Prepared AnnData not found: {path}")


def load_lamcore_programs(config: dict[str, Any]) -> tuple[dict[str, list[str]], dict[str, Any]]:
    path = resolve_file(config, "config/known_lam_programs.yaml")
    if path is None:
        return {}, {"status": "not_available", "path": ""}
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    programs: dict[str, list[str]] = {}
    for entry in payload.get("programs", []):
        name = str(entry.get("program_name", ""))
        if name in {"CORE2", "CORE3_identity"}:
            programs[name] = unique_genes([str(gene) for gene in entry.get("genes", [])])
    return programs, {"status": "available", "path": str(path), "programs": programs}


def load_formal_signature(config: dict[str, Any]) -> tuple[list[str], dict[str, Any]]:
    path = resolve_file(config, "data/raw/reference/LAM_core_signature_genes.csv")
    if path is None:
        return [], {"status": "not_available", "path": ""}
    table = pd.read_csv(path)
    gene_column = "Gene" if "Gene" in table.columns else str(table.columns[0])
    genes = unique_genes(table[gene_column].dropna().astype(str).tolist())
    return genes, {"status": "available", "path": str(path), "gene_column": gene_column, "n_genes": len(genes)}


def dense_values(matrix: Any) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def module_scores_from_log_x(
    prepared: ad.AnnData,
    modules: dict[str, list[str]],
    block_size: int = 10000,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Calculate mean log1p-normalized expression in bounded row blocks."""

    var_names = pd.Index(prepared.var_names.astype(str))
    lookup: dict[str, str] = {}
    for actual in var_names:
        lookup.setdefault(canonical_gene(actual), str(actual))
    resolved: dict[str, list[str]] = {}
    union: list[str] = []
    for name, genes in modules.items():
        present = [lookup[canonical_gene(gene)] for gene in unique_genes(genes) if canonical_gene(gene) in lookup]
        resolved[name] = list(dict.fromkeys(present))
        for gene in resolved[name]:
            if gene not in union:
                union.append(gene)
    scores = {name: np.full(prepared.n_obs, np.nan, dtype=np.float32) for name in modules}
    if not union:
        return scores, {"status": "no_module_genes_in_prepared", "resolved_genes": resolved}
    union_index = {gene: i for i, gene in enumerate(union)}
    module_indices = {name: [union_index[gene] for gene in genes] for name, genes in resolved.items()}
    for start in range(0, prepared.n_obs, block_size):
        stop = min(start + block_size, prepared.n_obs)
        block = prepared[start:stop, union].to_memory()
        values = dense_values(block.X)
        for name, indices in module_indices.items():
            if indices:
                scores[name][start:stop] = np.nanmean(values[:, indices], axis=1)
    return scores, {
        "status": "ok",
        "resolved_genes": resolved,
        "n_union_genes": len(union),
        "expression_source": "prepared AnnData X (library-size normalized log1p)",
        "block_size": block_size,
    }


def robust_z_by_dataset(values: np.ndarray, datasets: np.ndarray) -> tuple[np.ndarray, dict[str, dict[str, float]]]:
    result = np.full(values.shape, np.nan, dtype=np.float32)
    parameters: dict[str, dict[str, float]] = {}
    for dataset in sorted(set(datasets.astype(str))):
        mask = datasets.astype(str) == dataset
        finite = values[mask][np.isfinite(values[mask])]
        if finite.size == 0:
            parameters[dataset] = {"median": np.nan, "iqr": np.nan, "scale": np.nan, "n": 0}
            continue
        median = float(np.median(finite))
        q25, q75 = np.quantile(finite, [0.25, 0.75])
        iqr = float(q75 - q25)
        scale = iqr if iqr > 1e-8 else float(np.std(finite))
        scale = scale if scale > 1e-8 else 1.0
        result[mask] = ((values[mask] - median) / scale).astype(np.float32)
        parameters[dataset] = {"median": median, "iqr": iqr, "scale": scale, "n": int(finite.size)}
    return result, parameters


def nanmean_columns(values: list[np.ndarray]) -> np.ndarray:
    matrix = np.column_stack(values)
    valid = np.isfinite(matrix)
    denominator = valid.sum(axis=1)
    result = np.divide(np.nansum(matrix, axis=1), denominator, out=np.full(len(matrix), np.nan), where=denominator > 0)
    return result.astype(np.float32)


def sigmoid(values: np.ndarray) -> np.ndarray:
    clipped = np.clip(values, -30, 30)
    return 1.0 / (1.0 + np.exp(-clipped))


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


def quantile_or_default(values: np.ndarray, quantile: float, default: float = 0.0) -> float:
    finite = values[np.isfinite(values)]
    return float(np.quantile(finite, quantile)) if finite.size else default


def calibrate_thresholds(frame: pd.DataFrame, train_datasets: set[str], label: str) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    dataset_values = frame["dataset"].astype(str).to_numpy()
    in_train_lam = frame["condition"].astype(str).eq("LAM").to_numpy() & np.isin(dataset_values, list(train_datasets))
    normal = frame["normal_reference"].to_numpy(dtype=bool)
    positive = frame["positive_reference"].to_numpy(dtype=bool) & in_train_lam
    negative = frame["negative_reference"].to_numpy(dtype=bool) & (normal | in_train_lam)
    negative &= ~positive
    if positive.sum() == 0 or negative.sum() == 0:
        raise RuntimeError(f"Cannot calibrate {label}: positive={positive.sum()} negative={negative.sum()}")

    score_threshold, score_metrics = best_threshold(
        frame.loc[positive, "LAM_identity_score"].to_numpy(dtype=float),
        frame.loc[negative, "LAM_identity_score"].to_numpy(dtype=float),
    )
    identity_threshold, identity_metrics = best_threshold(
        frame.loc[positive, "identity_anchor_score"].to_numpy(dtype=float),
        frame.loc[negative, "identity_anchor_score"].to_numpy(dtype=float),
    )
    support_threshold, support_metrics = best_threshold(
        frame.loc[positive, "support_score"].to_numpy(dtype=float),
        frame.loc[negative, "support_score"].to_numpy(dtype=float),
    )
    # Boundary thresholds deliberately preserve a lower-evidence band around
    # the independently calibrated positive reference rather than selecting a
    # target number of old consensus states.
    boundary_score = min(score_threshold, quantile_or_default(frame.loc[positive, "LAM_identity_score"].to_numpy(dtype=float), 0.10, score_threshold))
    boundary_identity = min(identity_threshold, quantile_or_default(frame.loc[positive, "identity_anchor_score"].to_numpy(dtype=float), 0.10, identity_threshold))
    boundary_support = quantile_or_default(frame.loc[positive, "support_score"].to_numpy(dtype=float), 0.10, support_threshold)
    competing_cutoffs = {
        lineage: quantile_or_default(frame.loc[negative, f"{lineage}_score_z"].to_numpy(dtype=float), 0.95, 2.0)
        for lineage in NON_PERICYTE_LINEAGES
    }
    pericyte_cutoff = quantile_or_default(frame.loc[negative, "pericyte_VSMC_score_z"].to_numpy(dtype=float), 0.95, 2.0)
    thresholds = {
        "calibration_label": label,
        "train_datasets": sorted(train_datasets),
        "n_positive_reference": int(positive.sum()),
        "n_negative_reference": int(negative.sum()),
        "LAM_identity_score_core_threshold": score_threshold,
        "identity_anchor_core_threshold": identity_threshold,
        "support_core_threshold": support_threshold,
        "LAM_identity_score_boundary_threshold": boundary_score,
        "identity_anchor_boundary_threshold": boundary_identity,
        "support_boundary_threshold": boundary_support,
        "competing_lineage_95pct_cutoffs": competing_cutoffs,
        "pericyte_VSMC_95pct_cutoff": pericyte_cutoff,
        "score_metrics": score_metrics,
        "identity_metrics": identity_metrics,
        "support_metrics": support_metrics,
    }
    audit_rows: list[dict[str, Any]] = []
    for metric_name, threshold, metrics in [
        ("LAM_identity_score", score_threshold, score_metrics),
        ("identity_anchor_score", identity_threshold, identity_metrics),
        ("support_score", support_threshold, support_metrics),
    ]:
        audit_rows.append({"calibration": label, "metric": metric_name, "threshold": threshold, **metrics, "n_positive": int(positive.sum()), "n_negative": int(negative.sum())})
    for lineage, threshold in competing_cutoffs.items():
        audit_rows.append({"calibration": label, "metric": f"{lineage}_score_z_negative_q95", "threshold": threshold, "n_positive": int(positive.sum()), "n_negative": int(negative.sum())})
    audit_rows.append({"calibration": label, "metric": "pericyte_VSMC_score_z_negative_q95", "threshold": pericyte_cutoff, "n_positive": int(positive.sum()), "n_negative": int(negative.sum())})
    return thresholds, audit_rows


def assign_classes(frame: pd.DataFrame, thresholds: dict[str, Any]) -> pd.DataFrame:
    result = frame.copy()
    non_pericyte_dominant = np.zeros(len(result), dtype=bool)
    dominant_lineage = np.full(len(result), "none", dtype=object)
    dominant_value = np.full(len(result), -np.inf, dtype=float)
    for lineage in NON_PERICYTE_LINEAGES:
        values = result[f"{lineage}_score_z"].to_numpy(dtype=float)
        cutoff = float(thresholds["competing_lineage_95pct_cutoffs"][lineage])
        non_pericyte_dominant |= np.isfinite(values) & (values >= cutoff)
        better = np.isfinite(values) & (values > dominant_value)
        dominant_lineage[better] = lineage
        dominant_value[better] = values[better]
    pericyte = result["pericyte_VSMC_score_z"].to_numpy(dtype=float)
    identity = result["identity_anchor_score"].to_numpy(dtype=float)
    pericyte_only_exclusion = (
        np.isfinite(pericyte)
        & (pericyte >= float(thresholds["pericyte_VSMC_95pct_cutoff"]))
        & (identity < float(thresholds["identity_anchor_boundary_threshold"]))
    )
    competing_exclusion = non_pericyte_dominant | pericyte_only_exclusion
    core = (
        (result["LAM_identity_score"].to_numpy(dtype=float) >= float(thresholds["LAM_identity_score_core_threshold"]))
        & (identity >= float(thresholds["identity_anchor_core_threshold"]))
        & (result["support_score"].to_numpy(dtype=float) >= float(thresholds["support_core_threshold"]))
        & ~competing_exclusion
    )
    boundary = (
        ~core
        & (result["LAM_identity_score"].to_numpy(dtype=float) >= float(thresholds["LAM_identity_score_boundary_threshold"]))
        & (identity >= float(thresholds["identity_anchor_boundary_threshold"]))
        & (result["support_score"].to_numpy(dtype=float) >= float(thresholds["support_boundary_threshold"]))
        & ~competing_exclusion
    )
    result["LAM_core_candidate"] = core
    result["LAM_boundary_candidate"] = boundary
    result["non_LAM_like"] = ~(core | boundary)
    result["dominant_competing_lineage"] = dominant_lineage
    result["dominant_competing_lineage_score_z"] = dominant_value
    result["competing_lineage_exclusion"] = competing_exclusion
    result["pericyte_VSMC_conditional_exclusion"] = pericyte_only_exclusion
    result["identity_gate_reason"] = np.select(
        [
            core,
            boundary,
            competing_exclusion,
            identity < float(thresholds["identity_anchor_boundary_threshold"]),
            result["support_score"].to_numpy(dtype=float) < float(thresholds["support_boundary_threshold"]),
        ],
        [
            "LAM_identity_plus_support",
            "LAM_identity_boundary_evidence",
            "competing_lineage_dominant",
            "insufficient_LAM_identity_evidence",
            "insufficient_additional_support",
        ],
        default="low_composite_identity_score",
    )
    return result


def score_summary_rows(frame: pd.DataFrame, cohort_name: str, dataset: str, score_names: list[str]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for score_name in score_names:
        values = frame[score_name].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        if values.size == 0:
            continue
        rows.append(
            {
                "dataset": dataset,
                "cohort": cohort_name,
                "score": score_name,
                "n_cells": int(len(frame)),
                "n_nonnull": int(values.size),
                "mean": float(np.mean(values)),
                "median": float(np.median(values)),
                "q05": float(np.quantile(values, 0.05)),
                "q25": float(np.quantile(values, 0.25)),
                "q75": float(np.quantile(values, 0.75)),
                "q95": float(np.quantile(values, 0.95)),
                "LAM_core_fraction": float(frame["LAM_core_candidate"].mean()) if "LAM_core_candidate" in frame else np.nan,
                "LAM_boundary_fraction": float(frame["LAM_boundary_candidate"].mean()) if "LAM_boundary_candidate" in frame else np.nan,
                "non_LAM_like_fraction": float(frame["non_LAM_like"].mean()) if "non_LAM_like" in frame else np.nan,
            }
        )
    return rows


def compact_counts(values: pd.Series) -> str:
    counts = values.fillna("NA").astype(str).value_counts().to_dict()
    return json.dumps({str(key): int(value) for key, value in counts.items()}, ensure_ascii=False, sort_keys=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/state_modeling.yaml"))
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "results/stage16"))
    parser.add_argument("--block-size", type=int, default=10000)
    args = parser.parse_args()
    config = load_config(Path(args.config).resolve())
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    prepared_path = resolve_prepared_path(config)
    prepared = ad.read_h5ad(prepared_path, backed="r")
    obs_columns = set(prepared.obs.columns)
    metadata_columns = [
        column
        for column in [
            "cell_id", "source_cell_id", "dataset", "condition", "patient_id", "specimen_id",
            "assay", "cell_type", "total_counts", "n_genes_by_counts", "lam_candidate",
        ]
        if column in obs_columns
    ]
    for dataset in DATASETS:
        metadata_columns.extend(
            [
                f"upstream_{dataset}_candidate_source_author_style",
                f"upstream_{dataset}_candidate_source_formal_signature",
                f"upstream_{dataset}_core3_core3_identity_score",
                f"upstream_{dataset}_core3_core3_like",
                f"upstream_{dataset}_adata_lamcore_score_777_consistency",
            ]
        )
    metadata_columns = list(dict.fromkeys(column for column in metadata_columns if column in obs_columns))
    metadata = prepared.obs[metadata_columns].copy().reset_index(names="analysis_cell_id")
    metadata["analysis_cell_id"] = metadata["analysis_cell_id"].astype(str)
    metadata["dataset"] = metadata["dataset"].astype(str)
    metadata["condition"] = metadata["condition"].astype(str)

    lamcore_programs, lamcore_manifest = load_lamcore_programs(config)
    formal_genes, formal_manifest = load_formal_signature(config)
    core2_genes = lamcore_programs.get("CORE2", ["PMEL", "MLANA", "MITF", "TYR", "DCT", "GPNMB", "CTSK", "VEGFD"])
    core3_genes = lamcore_programs.get("CORE3_identity", ["ACTA2", "PMEL", "MLANA", "MITF", "CTSK", "VEGFD", "ESR1"])
    modules: dict[str, list[str]] = {
        "melanocytic": MELANOCYTIC_GENES,
        "support": SUPPORT_GENES,
        "lamcore_CORE2": core2_genes,
        "CORE3_identity": core3_genes,
    }
    if formal_genes:
        modules["LAMCORE_777"] = formal_genes
    modules.update(COMPETING_GENES)
    module_scores, module_manifest = module_scores_from_log_x(prepared, modules, block_size=int(args.block_size))

    n_rows = len(metadata)
    datasets = metadata["dataset"].to_numpy(dtype=str)
    condition = metadata["condition"].to_numpy(dtype=str)
    for name, values in module_scores.items():
        metadata[f"{name}_score_raw"] = values

    # Dataset-specific upstream fields are aligned into one canonical namespace.
    for name in ["source_author_style", "source_formal_signature", "core3_like"]:
        metadata[name] = False
    for name in ["core3_identity_upstream", "lamcore_777_upstream"]:
        metadata[name] = np.nan
    for dataset in DATASETS:
        mask = datasets == dataset
        for canonical, suffix in [
            ("source_author_style", "candidate_source_author_style"),
            ("source_formal_signature", "candidate_source_formal_signature"),
            ("core3_like", "core3_core3_like"),
        ]:
            column = f"upstream_{dataset}_{suffix}"
            if column in metadata:
                values = bool_array(metadata.loc[mask, column])
                metadata.loc[mask, canonical] = values
        score_column = f"upstream_{dataset}_core3_core3_identity_score"
        if score_column in metadata:
            metadata.loc[mask, "core3_identity_upstream"] = pd.to_numeric(metadata.loc[mask, score_column], errors="coerce").to_numpy()
        formal_column = f"upstream_{dataset}_adata_lamcore_score_777_consistency"
        if formal_column in metadata:
            metadata.loc[mask, "lamcore_777_upstream"] = pd.to_numeric(metadata.loc[mask, formal_column], errors="coerce").to_numpy()

    # Continuous components.  Z scores are calculated within dataset across
    # all prepared cells, including normal/control cells, so no old state is
    # involved in calibration.
    z_parameters: dict[str, Any] = {}
    for name in ["melanocytic", "support", "lamcore_CORE2", "CORE3_identity", "LAMCORE_777"] + list(COMPETING_GENES):
        raw_column = f"{name}_score_raw"
        if raw_column in metadata:
            z, params = robust_z_by_dataset(metadata[raw_column].to_numpy(dtype=float), datasets)
            metadata[f"{name}_score_z"] = z
            z_parameters[name] = params
    upstream_core3_z, params = robust_z_by_dataset(metadata["core3_identity_upstream"].to_numpy(dtype=float), datasets)
    metadata["core3_identity_upstream_z"] = upstream_core3_z
    z_parameters["core3_identity_upstream"] = params
    upstream_777_z, params = robust_z_by_dataset(metadata["lamcore_777_upstream"].to_numpy(dtype=float), datasets)
    metadata["lamcore_777_upstream_z"] = upstream_777_z
    z_parameters["lamcore_777_upstream"] = params

    metadata["melanocytic_identity_score"] = metadata["melanocytic_score_raw"]
    metadata["lamcore_identity_score"] = nanmean_columns(
        [
            metadata["lamcore_CORE2_score_z"].to_numpy(dtype=float),
            metadata["CORE3_identity_score_z"].to_numpy(dtype=float),
            metadata["core3_identity_upstream_z"].to_numpy(dtype=float),
            metadata["lamcore_777_upstream_z"].to_numpy(dtype=float),
        ]
    )
    metadata["melanocytic_identity_score_z"] = metadata["melanocytic_score_z"]
    metadata["lamcore_identity_score_z"] = metadata["lamcore_identity_score"]
    metadata["identity_anchor_score"] = nanmean_columns(
        [metadata["melanocytic_identity_score_z"].to_numpy(dtype=float), metadata["lamcore_identity_score_z"].to_numpy(dtype=float)]
    )
    metadata["support_score"] = metadata["support_score_z"]
    metadata["non_pericyte_competing_score"] = np.nanmax(
        metadata[[f"{name}_score_z" for name in NON_PERICYTE_LINEAGES]].to_numpy(dtype=float), axis=1
    )
    pericyte_z = metadata["pericyte_VSMC_score_z"].to_numpy(dtype=float)
    identity_anchor = metadata["identity_anchor_score"].to_numpy(dtype=float)
    pericyte_penalty = np.maximum(pericyte_z, 0.0) * (1.0 - sigmoid(identity_anchor))
    non_pericyte_penalty = np.maximum(metadata["non_pericyte_competing_score"].to_numpy(dtype=float), 0.0)
    metadata["pericyte_VSMC_conditional_penalty"] = pericyte_penalty
    metadata["competing_lineage_penalty"] = non_pericyte_penalty + 0.25 * pericyte_penalty
    metadata["LAM_identity_score"] = (
        0.60 * metadata["identity_anchor_score"].to_numpy(dtype=float)
        + 0.25 * metadata["support_score"].to_numpy(dtype=float)
        - 0.20 * non_pericyte_penalty
        - 0.10 * pericyte_penalty
    ).astype(np.float32)

    metadata["normal_reference"] = condition != "LAM"
    metadata["positive_reference"] = (
        metadata["source_author_style"].to_numpy(dtype=bool)
        | metadata["source_formal_signature"].to_numpy(dtype=bool)
        | metadata["core3_like"].to_numpy(dtype=bool)
    ) & (condition == "LAM")
    metadata["positive_reference_type"] = np.select(
        [
            metadata["source_author_style"].to_numpy(dtype=bool),
            metadata["source_formal_signature"].to_numpy(dtype=bool),
            metadata["core3_like"].to_numpy(dtype=bool),
        ],
        ["source_author_style", "source_formal_signature", "upstream_CORE3_like"],
        default="none",
    )
    # Strong competing references are signature-defined, not state-defined.
    # A positive upstream LAM reference is excluded from the negative set.
    metadata["strong_competing_reference"] = (
        (metadata["non_pericyte_competing_score"].to_numpy(dtype=float) >= 2.0)
        | (pericyte_z >= 2.0)
    ) & ~metadata["positive_reference"].to_numpy(dtype=bool)
    metadata["negative_reference"] = (
        metadata["normal_reference"].to_numpy(dtype=bool)
        | metadata["strong_competing_reference"].to_numpy(dtype=bool)
    ) & ~metadata["positive_reference"].to_numpy(dtype=bool)

    frame = metadata[condition == "LAM"].copy().reset_index(drop=True)
    # Calibration needs normal/competing reference rows too, so preserve a
    # separate all-cell frame with the same score columns.
    all_frame = metadata.copy()
    full_thresholds, calibration_rows = calibrate_thresholds(all_frame, set(DATASETS), "full_four_dataset")
    full_assigned = assign_classes(frame, full_thresholds)

    # LODO calibrations: held-out data never contributes positive/negative
    # reference values to its own threshold selection.
    lodo_rows: list[dict[str, Any]] = []
    lodo_thresholds: dict[str, Any] = {}
    for held_out in DATASETS:
        train = set(DATASETS) - {held_out}
        label = f"leave_out_{held_out}"
        thresholds, rows = calibrate_thresholds(all_frame, train, label)
        lodo_thresholds[held_out] = thresholds
        calibration_rows.extend(rows)
        held_mask = frame["dataset"].astype(str).eq(held_out).to_numpy()
        held = assign_classes(frame.loc[held_mask].copy(), thresholds)
        held_positive = held["positive_reference"].to_numpy(dtype=bool)
        held_negative = held["negative_reference"].to_numpy(dtype=bool)
        lodo_rows.append(
            {
                "held_out_dataset": held_out,
                "train_datasets": json.dumps(sorted(train)),
                "n_held_out_lam_cells": int(len(held)),
                "n_positive_reference_held_out": int(held_positive.sum()),
                "positive_reference_core_recall": float(held.loc[held_positive, "LAM_core_candidate"].mean()) if held_positive.any() else np.nan,
                "positive_reference_boundary_or_core_recall": float((held.loc[held_positive, "LAM_core_candidate"] | held.loc[held_positive, "LAM_boundary_candidate"]).mean()) if held_positive.any() else np.nan,
                "negative_reference_held_out": int(held_negative.sum()),
                "negative_reference_non_LAM_like_fraction": float(held.loc[held_negative, "non_LAM_like"].mean()) if held_negative.any() else np.nan,
                "LAM_core_fraction": float(held["LAM_core_candidate"].mean()),
                "LAM_boundary_fraction": float(held["LAM_boundary_candidate"].mean()),
                "non_LAM_like_fraction": float(held["non_LAM_like"].mean()),
                "mean_LAM_identity_score": float(held["LAM_identity_score"].mean()),
                "median_LAM_identity_score": float(held["LAM_identity_score"].median()),
                "n_competing_lineage_excluded": int(held["competing_lineage_exclusion"].sum()),
            }
        )

    # Diagnostic output is restricted to LAM cells, as requested.
    evidence_columns = [
        "analysis_cell_id", "cell_id", "source_cell_id", "dataset", "condition", "patient_id", "specimen_id",
        "assay", "cell_type", "total_counts", "n_genes_by_counts", "source_author_style", "source_formal_signature",
        "core3_like", "positive_reference", "positive_reference_type", "melanocytic_identity_score",
        "melanocytic_identity_score_z", "lamcore_identity_score", "lamcore_identity_score_z", "core3_identity_upstream",
        "core3_identity_upstream_z", "lamcore_777_upstream", "lamcore_777_upstream_z", "support_score",
        "support_score_raw", "identity_anchor_score", "LAM_identity_score", "competing_lineage_penalty",
        "non_pericyte_competing_score", "pericyte_VSMC_conditional_penalty", "strong_competing_reference",
        "negative_reference", "LAM_core_candidate", "LAM_boundary_candidate", "non_LAM_like",
        "dominant_competing_lineage", "dominant_competing_lineage_score_z", "competing_lineage_exclusion",
        "pericyte_VSMC_conditional_exclusion", "identity_gate_reason",
    ]
    for lineage in list(COMPETING_GENES):
        evidence_columns.extend([f"{lineage}_score_raw", f"{lineage}_score_z"])
    evidence_columns = list(dict.fromkeys(column for column in evidence_columns if column in full_assigned.columns))
    full_assigned[evidence_columns].to_csv(output_dir / "cell_identity_evidence.csv", index=False)

    score_names = [
        "melanocytic_identity_score", "lamcore_identity_score", "identity_anchor_score", "support_score",
        "LAM_identity_score", "non_pericyte_competing_score", "pericyte_VSMC_conditional_penalty",
    ] + [f"{lineage}_score_z" for lineage in COMPETING_GENES]
    summary_rows: list[dict[str, Any]] = []
    cohort_masks = {
        "LAM_all": condition == "LAM",
        "author_style_positive": metadata["source_author_style"].to_numpy(dtype=bool),
        "formal_signature_positive": metadata["source_formal_signature"].to_numpy(dtype=bool),
        "independent_CORE3_positive": metadata["core3_like"].to_numpy(dtype=bool),
        "normal_reference": metadata["normal_reference"].to_numpy(dtype=bool),
        "strong_competing_reference": metadata["strong_competing_reference"].to_numpy(dtype=bool),
    }
    for cohort, mask in cohort_masks.items():
        for dataset in sorted(set(datasets[mask])):
            sub = metadata.loc[mask & (datasets == dataset)].copy()
            summary_rows.extend(score_summary_rows(sub, cohort, dataset, score_names))
    pd.DataFrame(summary_rows).to_csv(output_dir / "identity_score_by_dataset.csv", index=False)

    calibration_df = pd.DataFrame(calibration_rows)
    calibration_df.to_csv(output_dir / "reference_calibration.csv", index=False)
    pd.DataFrame(lodo_rows).to_csv(output_dir / "leave_one_dataset_out_validation.csv", index=False)

    # Freeze the new gate first; only now map it back to the existing 20 states.
    state_path = PROJECT_ROOT / str(config["outputs"]["step7_dir"]) / "state_consensus_with_upstream_annotations.csv"
    if state_path.exists():
        old_states = pd.read_csv(state_path, dtype=str)[["analysis_cell_id", "consensus_state"]].drop_duplicates("analysis_cell_id")
        old_states = old_states.merge(
            full_assigned[["analysis_cell_id", "LAM_core_candidate", "LAM_boundary_candidate", "non_LAM_like", "source_author_style", "LAM_identity_score", "dominant_competing_lineage"]],
            on="analysis_cell_id",
            how="left",
        )
        state_rows: list[dict[str, Any]] = []
        for state, sub in old_states.groupby("consensus_state", observed=True):
            dominant = sub["dominant_competing_lineage"].value_counts().index[0] if len(sub) else "none"
            state_rows.append(
                {
                    "consensus_state": str(state),
                    "cells": int(len(sub)),
                    "LAM_core_fraction": float(sub["LAM_core_candidate"].mean()),
                    "LAM_boundary_fraction": float(sub["LAM_boundary_candidate"].mean()),
                    "non_LAM_like_fraction": float(sub["non_LAM_like"].mean()),
                    "author_supported_fraction": float(sub["source_author_style"].map(as_bool).mean()),
                    "median_LAM_identity_score": float(pd.to_numeric(sub["LAM_identity_score"], errors="coerce").median()),
                    "dominant_competing_lineage": str(dominant),
                    "competing_lineage_counts": compact_counts(sub["dominant_competing_lineage"]),
                }
            )
        state_summary = pd.DataFrame(state_rows)
        if not state_summary.empty:
            state_summary["_sort"] = pd.to_numeric(state_summary["consensus_state"], errors="coerce")
            state_summary = state_summary.sort_values(["_sort", "consensus_state"]).drop(columns="_sort")
        state_summary.to_csv(output_dir / "new_candidate_by_old_state.csv", index=False)
    else:
        pd.DataFrame(columns=["consensus_state", "cells", "LAM_core_fraction", "LAM_boundary_fraction", "non_LAM_like_fraction", "author_supported_fraction", "median_LAM_identity_score", "dominant_competing_lineage", "competing_lineage_counts"]).to_csv(output_dir / "new_candidate_by_old_state.csv", index=False)

    assignment_columns = [
        "analysis_cell_id", "cell_id", "source_cell_id", "dataset", "condition", "patient_id", "assay",
        "source_author_style", "source_formal_signature", "core3_like", "positive_reference", "positive_reference_type",
        "melanocytic_identity_score", "lamcore_identity_score", "support_score", "identity_anchor_score",
        "LAM_identity_score", "LAM_core_candidate", "LAM_boundary_candidate", "non_LAM_like",
        "identity_gate_reason", "dominant_competing_lineage", "competing_lineage_exclusion",
        "pericyte_VSMC_conditional_exclusion",
    ]
    full_assigned[assignment_columns].to_csv(output_dir / "new_candidate_assignment.csv", index=False)

    manifest = {
        "script": str(Path(__file__).resolve()),
        "read_only_existing_artifacts": True,
        "prepared_h5ad": str(prepared_path),
        "n_prepared_cells": int(len(metadata)),
        "n_lam_cells": int((condition == "LAM").sum()),
        "n_normal_or_control_cells": int((condition != "LAM").sum()),
        "n_datasets_lam": int(len(set(datasets[condition == "LAM"]))),
        "n_old_states_used_only_posthoc": int(old_states["consensus_state"].nunique()) if state_path.exists() else 0,
        "no_scvi_training": True,
        "no_reclustering": True,
        "gene_aliases": ALIASES,
        "identity_anchor_genes": MELANOCYTIC_GENES,
        "support_genes": SUPPORT_GENES,
        "competing_gene_sets": COMPETING_GENES,
        "lamcore_program_manifest": lamcore_manifest,
        "formal_signature_manifest": formal_manifest,
        "module_manifest": module_manifest,
        "robust_z_parameters": z_parameters,
        "positive_reference_definition": "source_author_style OR source_formal_signature OR upstream CORE3_like, condition=LAM",
        "negative_reference_definition": "normal/control OR strong signature-defined competing lineage, excluding positive reference",
        "assignment_formula": "0.60*identity_anchor + 0.25*support - 0.20*non_pericyte_penalty - 0.10*conditional_pericyte_penalty",
        "full_thresholds": full_thresholds,
        "n_positive_reference": int(metadata["positive_reference"].sum()),
        "n_negative_reference": int(metadata["negative_reference"].sum()),
        "full_assignment_counts": {column: int(full_assigned[column].sum()) for column in ["LAM_core_candidate", "LAM_boundary_candidate", "non_LAM_like", "competing_lineage_exclusion"]},
        "output_files": sorted(path.name for path in output_dir.iterdir() if path.is_file()),
    }
    (output_dir / "stage16_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    report = [
        "# Stage 16 — LAM candidate identity gate audit",
        "",
        f"LAM cells scored: {int((condition == 'LAM').sum()):,}",
        f"Normal/control cells used for calibration: {int((condition != 'LAM').sum()):,}",
        f"Independent positive references: {int(metadata['positive_reference'].sum()):,}",
        f"Negative references: {int(metadata['negative_reference'].sum()):,}",
        "",
        "## Frozen evidence model",
        "",
        "Identity anchors are PMEL/MLANA/MITF plus CORE2/CORE3 module evidence and available upstream CORE3 continuous scores. ACTA2/ESR1/VEGFD/CTSK are supportive evidence only. Competing lineage modules are continuous penalties; pericyte/VSMC is conditional on weak LAM identity and is not a standalone exclusion rule.",
        "",
        "The formal 777-gene LAMCORE CSV was not available in the resolved input roots, so it was not fabricated or used. The report records this as an unavailable optional reference.",
        "",
        "## Assignment counts",
        "",
        f"- LAM_core_candidate: {int(full_assigned['LAM_core_candidate'].sum()):,}",
        f"- LAM_boundary_candidate: {int(full_assigned['LAM_boundary_candidate'].sum()):,}",
        f"- non_LAM_like: {int(full_assigned['non_LAM_like'].sum()):,}",
        f"- competing-lineage exclusion: {int(full_assigned['competing_lineage_exclusion'].sum()):,}",
        "",
        "Thresholds were calibrated without consensus_state and the existing states were joined only for the final diagnostic table. No Step 7–13 script and no scVI training was called.",
    ]
    (output_dir / "identity_gate_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    prepared.file.close()
    print(f"Scored {int((condition == 'LAM').sum()):,} LAM cells")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()

