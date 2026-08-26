"""Shared read-only loaders and expression-state utilities."""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


PROJECT_ROOT = Path(__file__).resolve().parents[1]


@dataclass
class ExpressionData:
    obs: pd.DataFrame
    target_genes: list[str]
    available_genes: list[str]
    raw_counts: np.ndarray
    normalized: np.ndarray
    gene_to_col: dict[str, int]
    raw_counts_available: bool
    source_path: str


def load_yaml(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text())


def load_project_config() -> dict[str, Any]:
    return load_yaml(PROJECT_ROOT / "config" / "project.yaml")


def load_source_manifest() -> dict[str, Any]:
    return load_yaml(PROJECT_ROOT / "config" / "source_manifest.yaml")


def load_signatures() -> dict[str, dict[str, Any]]:
    return load_yaml(PROJECT_ROOT / "config" / "signatures.yaml")


def resolve_source(path: str | Path, source_root: str | Path | None = None) -> Path:
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    root = Path(source_root or load_source_manifest()["source_root"])
    if not root.is_absolute():
        root = PROJECT_ROOT / root
    return root / candidate


def project_relative(path: str | Path) -> str:
    """Return a path relative to this project's root for manifests and reports."""

    candidate = Path(path)
    resolved = candidate if candidate.is_absolute() else PROJECT_ROOT / candidate
    return Path(os.path.relpath(resolved.resolve(), PROJECT_ROOT.resolve())).as_posix()


def ensure_output_path(path: Path) -> Path:
    resolved = path.resolve()
    if PROJECT_ROOT.resolve() not in resolved.parents and resolved != PROJECT_ROOT.resolve():
        raise ValueError(f"Refusing to write outside project root: {resolved}")
    resolved.parent.mkdir(parents=True, exist_ok=True)
    return resolved


def materialize_matrix(matrix: Any) -> np.ndarray:
    if hasattr(matrix, "to_memory"):
        matrix = matrix.to_memory()
    if sp.issparse(matrix):
        return matrix.toarray().astype(float, copy=False)
    return np.asarray(matrix, dtype=float)


def _gene_symbols(adata: ad.AnnData) -> np.ndarray:
    for column in ("gene_symbol_upper", "gene_symbol"):
        if column in adata.var.columns:
            return adata.var[column].astype(str).str.upper().to_numpy()
    return adata.var_names.astype(str).str.upper().to_numpy()


def load_expression(path: Path, target_genes: list[str]) -> ExpressionData:
    """Read only selected genes from an h5ad and close the source file."""

    target = [str(g).upper() for g in dict.fromkeys(target_genes)]
    source = ad.read_h5ad(path, backed="r")
    try:
        symbols = _gene_symbols(source)
        gene_to_source_col: dict[str, int] = {}
        for index, symbol in enumerate(symbols):
            gene_to_source_col.setdefault(symbol, index)
        available = [gene for gene in target if gene in gene_to_source_col]
        indices = [gene_to_source_col[gene] for gene in available]
        if indices:
            view = source[:, indices]
            raw_available = "counts" in source.layers
            matrix = view.layers["counts"] if raw_available else view.X
            values = materialize_matrix(matrix)
        else:
            raw_available = False
            values = np.zeros((source.n_obs, 0), dtype=float)
        obs = source.obs.copy()
        if "total_counts" in obs.columns:
            totals = pd.to_numeric(obs["total_counts"], errors="coerce").fillna(0).to_numpy(float)
        else:
            totals = values.sum(axis=1)
        totals = np.maximum(totals, 1.0)
        if raw_available:
            normalized = np.log1p(values / totals[:, None] * 1e4)
        else:
            # Existing converted objects are already log-normalized when a
            # counts layer is absent. The manifest records this limitation.
            normalized = values
        return ExpressionData(
            obs=obs,
            target_genes=target,
            available_genes=available,
            raw_counts=values,
            normalized=normalized,
            gene_to_col={gene: i for i, gene in enumerate(available)},
            raw_counts_available=raw_available,
            source_path=str(path),
        )
    finally:
        if getattr(source, "file", None) is not None:
            source.file.close()


def load_pool_labels(path: Path | None, obs_index: pd.Index) -> pd.DataFrame:
    result = pd.DataFrame(index=obs_index)
    for column in ("pool_high_confidence", "pool_broad_lam_like", "pool_unrestricted_lam"):
        result[column] = False
    if path is None or not path.exists():
        result["identity_pool"] = "other"
        return result
    table = pd.read_csv(path)
    cell_column = "cell_id" if "cell_id" in table.columns else table.columns[0]
    table = table.set_index(cell_column)
    table.index = table.index.astype(str)
    aligned = table.reindex(obs_index.astype(str))
    for column in result.columns:
        if column in aligned:
            result[column] = aligned[column].fillna(False).astype(bool).to_numpy()
    result["identity_pool"] = np.select(
        [result["pool_high_confidence"], result["pool_broad_lam_like"], result["pool_unrestricted_lam"]],
        ["high_confidence", "broad_lam_like", "unrestricted_lam"],
        default="other",
    )
    return result


def add_metadata(expr: ExpressionData, dataset: str, spec: dict[str, Any], pool_path: Path | None) -> pd.DataFrame:
    obs = expr.obs.copy()
    obs.index = obs.index.astype(str)
    patient_column = spec.get("patient_column", "patient_id")
    if patient_column in obs:
        obs["patient_id"] = obs[patient_column].astype(str)
    elif "donor_id" in obs:
        obs["patient_id"] = obs["donor_id"].astype(str)
    else:
        obs["patient_id"] = "unknown_patient"
    obs["dataset"] = dataset
    obs["assay_label"] = obs[spec["assay_column"]].astype(str) if spec.get("assay_column") in obs else "unknown"
    obs["source_cell_id"] = obs.index
    pools = load_pool_labels(pool_path, obs.index)
    return pd.concat([obs, pools], axis=1)


def expression_intervals(
    expr: ExpressionData,
    positive_quantile: float,
    min_positive_observations: int,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for gene in expr.available_genes:
        col = expr.gene_to_col[gene]
        positive = expr.normalized[:, col][expr.raw_counts[:, col] > 0]
        if len(positive) >= min_positive_observations:
            low_cut = float(np.quantile(positive, positive_quantile))
            high_cut = float(np.quantile(positive, 1.0 - positive_quantile))
            stable = bool(np.isfinite(low_cut) and np.isfinite(high_cut) and high_cut > low_cut)
        else:
            low_cut = np.nan
            high_cut = np.nan
            stable = False
        rows.append({
            "gene": gene,
            "gene_in_matrix": True,
            "n_positive_observations": int(len(positive)),
            "low_expression_upper_bound": low_cut,
            "high_expression_lower_bound": high_cut,
            "expression_interval_stable": stable,
            "interval_basis": "nonzero normalized expression within dataset x assay",
        })
    for gene in expr.target_genes:
        if gene not in expr.gene_to_col:
            rows.append({
                "gene": gene,
                "gene_in_matrix": False,
                "n_positive_observations": 0,
                "low_expression_upper_bound": np.nan,
                "high_expression_lower_bound": np.nan,
                "expression_interval_stable": False,
                "interval_basis": "not_assayed",
            })
    return pd.DataFrame(rows).drop_duplicates("gene")


def gene_state_matrix(expr: ExpressionData, intervals: pd.DataFrame) -> np.ndarray:
    lookup = intervals.set_index("gene").to_dict("index")
    states = np.full((expr.obs.shape[0], len(expr.target_genes)), "not_assayed", dtype=object)
    for target_index, gene in enumerate(expr.target_genes):
        if gene not in expr.gene_to_col:
            continue
        col = expr.gene_to_col[gene]
        counts = expr.raw_counts[:, col]
        values = expr.normalized[:, col]
        state = np.full(len(counts), "not_detected", dtype=object)
        positive = counts > 0
        row = lookup[gene]
        if row["expression_interval_stable"]:
            state[positive & (values < row["low_expression_upper_bound"])] = "detected_low"
            state[positive & (values >= row["high_expression_lower_bound"])] = "detected_high"
            state[positive & (values >= row["low_expression_upper_bound"]) & (values < row["high_expression_lower_bound"])] = "detected_unclassified"
        else:
            state[positive] = "detected_unclassified"
        states[:, target_index] = state
    return states


def score_modules(
    expr: ExpressionData,
    signatures: dict[str, dict[str, Any]],
    min_available_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    scores = pd.DataFrame(index=expr.obs.index)
    availability_rows: list[dict[str, Any]] = []
    for module, spec in signatures.items():
        genes = [str(g).upper() for g in spec.get("genes", [])]
        available = [g for g in genes if g in expr.gene_to_col]
        cols = [expr.gene_to_col[g] for g in available]
        scores[f"module_{module}"] = expr.normalized[:, cols].mean(axis=1) if cols else np.nan
        detected = (expr.raw_counts[:, cols] > 0).sum(axis=1) if cols else np.zeros(expr.obs.shape[0], dtype=int)
        scores[f"module_{module}__n_available"] = len(available)
        scores[f"module_{module}__n_detected"] = detected
        if not available:
            status = "not_assayed"
        elif len(available) < min_available_genes:
            status = "insufficient_coverage"
        elif len(available) < len(genes):
            status = "partial_coverage"
        else:
            status = "analyzable"
        scores[f"module_{module}__status"] = status
        availability_rows.append({
            "module": module,
            "n_genes_requested": len(genes),
            "n_genes_available": len(available),
            "available_genes": ";".join(available),
            "missing_genes": ";".join(g for g in genes if g not in available),
            "module_status": status,
        })
    return scores, pd.DataFrame(availability_rows)


def quality_status(obs: pd.DataFrame, config: dict[str, Any]) -> pd.Series:
    min_genes = config["expression"]["low_quality_min_genes"]
    min_counts = config["expression"]["low_quality_min_counts"]
    genes = pd.to_numeric(obs.get("n_genes_by_counts", pd.Series(np.nan, index=obs.index)), errors="coerce")
    counts = pd.to_numeric(obs.get("total_counts", pd.Series(np.nan, index=obs.index)), errors="coerce")
    return pd.Series(np.where((genes < min_genes) | (counts < min_counts), "low_quality", "usable"), index=obs.index)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    ensure_output_path(path).write_text(json.dumps(payload, ensure_ascii=False, indent=2, default=str))
