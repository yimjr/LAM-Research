#!/usr/bin/env python3
"""Audit the inherited high-confidence LAM candidate identity.

This script is intentionally read-only with respect to existing analysis
artifacts.  It reads the Step 7 consensus annotation, the four original
candidate_pool_labels.csv files, and the prepared AnnData.  New audit tables
are written only below ``results/stage15``.

Important provenance note
-------------------------
    The State Modeling import has already canonicalized FIGF -> VEGFD in the
    AnnData counts layer.  The upstream ``marker_expr_*`` columns are retained
    from the upstream marker extraction and are therefore used to audit the
    original candidate rule and the original FIGF/VEGFD detections.  The
current counts layer is also read for all canonical genes that remain in the
prepared AnnData; its VEGFD value is the alias-merged count and is not used as
an independent original FIGF count.
"""

from __future__ import annotations

import argparse
import json
import os
from collections import Counter
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
import yaml


PROJECT_ROOT = Path(__file__).resolve().parent
DATASETS = ["GSE135851", "GSE190260", "GSE217108", "GSE302356"]
MARKER_GENES = ["PMEL", "MLANA", "MITF", "ACTA2", "ESR1", "FIGF", "VEGFD", "CTSK"]
BOOL_TRUE = {"1", "true", "t", "yes", "y"}
BOOL_FALSE = {"0", "false", "f", "no", "n", "", "nan", "none", "na", "<na>"}


def as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    text = str(value).strip().lower()
    if text in BOOL_TRUE:
        return True
    if text in BOOL_FALSE:
        return False
    return False


def bool_series(values: pd.Series) -> pd.Series:
    return values.map(as_bool).astype(bool)


def numeric_series(values: pd.Series, default: float = np.nan) -> pd.Series:
    result = pd.to_numeric(values, errors="coerce")
    return result.fillna(default)


def compact_counts(values: pd.Series) -> str:
    counts = values.fillna("NA").astype(str).value_counts(dropna=False).to_dict()
    return json.dumps({str(k): int(v) for k, v in counts.items()}, ensure_ascii=False, sort_keys=True)


def json_value(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, default=str)


def load_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def configured_roots(config: dict[str, Any]) -> list[Path]:
    roots = [PROJECT_ROOT / str(root) for root in config.get("input_roots", [])]
    # The current Windows-mounted workspace uses this sibling spelling.  It
    # is a read-only fallback and does not change the project configuration.
    roots.extend(
        [
            Path("/mnt/e/LAM-Research/LAM-Cell-Research"),
            Path("/mnt/e/LAM-Research/data-temp"),
            PROJECT_ROOT / "data/upstream",
        ]
    )
    unique: list[Path] = []
    seen: set[str] = set()
    for root in roots:
        key = str(root.resolve())
        if key not in seen:
            unique.append(root)
            seen.add(key)
    return unique


def resolve_candidate_path(config: dict[str, Any], dataset: str) -> Path | None:
    annotation_dir = str(config["datasets"][dataset]["annotation_dir"])
    for root in configured_roots(config):
        path = root / annotation_dir / "candidate_pool_labels.csv"
        if path.exists():
            return path
    return None


def resolve_prepared_path(config: dict[str, Any]) -> Path:
    configured = PROJECT_ROOT / str(config["outputs"]["prepared_h5ad"])
    if configured.exists():
        return configured
    candidates = [
        PROJECT_ROOT / "data/processed/state_model_prepared.h5ad",
        PROJECT_ROOT / "data/upstream/state_model_prepared.h5ad",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("state_model_prepared.h5ad was not found")


def resolve_original_dataset_path(config: dict[str, Any], dataset: str) -> Path | None:
    """Resolve the pre-alias source AnnData used for raw-count auditing."""

    spec = config["datasets"][dataset]
    candidates = [str(value) for value in spec.get("h5ad_candidates", [])]
    if dataset == "GSE135851":
        candidates.extend(
            [
                "GSE135851_core_reproduction.h5ad",
                "GSE135851_core_baseline_qc.h5ad",
            ]
        )
    else:
        candidates.append(f"{dataset}.h5ad")
    roots = configured_roots(config) + [Path("/mnt/e/LAM-Research/data-temp")]
    for root in roots:
        for candidate in candidates:
            path = root / candidate
            if path.exists():
                return path
        # data-temp stores the files by basename rather than under the source
        # project's data/processed path.
        for basename in [
            "GSE135851_core_reproduction.h5ad" if dataset == "GSE135851" else f"{dataset}.h5ad",
            "GSE135851_core_baseline_qc.h5ad" if dataset == "GSE135851" else "",
        ]:
            if basename:
                path = root / basename
                if path.exists():
                    return path
    return None


def state_sort(frame: pd.DataFrame, column: str = "consensus_state") -> pd.DataFrame:
    result = frame.copy()
    result["_numeric_state"] = pd.to_numeric(result[column], errors="coerce")
    return result.sort_values(["_numeric_state", column], na_position="last").drop(columns="_numeric_state")


def dense_values(matrix: Any) -> np.ndarray:
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def read_selected_counts(
    prepared: ad.AnnData,
    selected_ids: pd.Index,
    genes: list[str],
) -> tuple[np.ndarray, list[str], dict[str, Any]]:
    """Read only selected rows and marker columns from a backed AnnData."""

    if "counts" not in prepared.layers:
        return np.full((len(selected_ids), len(genes)), np.nan, dtype=np.float32), [], {
            "status": "missing_counts_layer",
            "available_genes": [],
        }

    var_names = pd.Index(prepared.var_names.astype(str))
    available = [gene for gene in genes if gene in var_names]
    result = np.full((len(selected_ids), len(genes)), np.nan, dtype=np.float32)
    if not available:
        return result, [], {"status": "no_marker_genes_in_counts", "available_genes": []}

    row_positions = pd.Index(prepared.obs_names.astype(str)).get_indexer(selected_ids)
    if (row_positions < 0).any():
        raise RuntimeError("Selected consensus cells are not all present in prepared AnnData")
    order = np.argsort(row_positions)
    sorted_rows = row_positions[order]
    try:
        subset = prepared[sorted_rows, available].to_memory()
        values_sorted = dense_values(subset.layers["counts"])
        values = np.empty_like(values_sorted)
        values[order] = values_sorted
        for j, gene in enumerate(available):
            result[:, genes.index(gene)] = values[:, j]
    except Exception as exc:  # pragma: no cover - only reached for malformed H5AD indexing
        return result, available, {
            "status": "counts_read_error",
            "available_genes": available,
            "error": f"{type(exc).__name__}: {exc}",
        }
    return result, available, {"status": "ok", "available_genes": available}


def read_original_raw_counts(
    config: dict[str, Any],
    current: pd.DataFrame,
    marker_column_presence: dict[str, list[str]],
) -> tuple[np.ndarray, dict[str, Any]]:
    """Read original, pre-alias raw counts one dataset at a time.

    The input is deliberately processed sequentially.  Only selected
    consensus cells and the eight marker genes are materialized in memory.
    Missing genes are represented as zero when the original marker column was
    present (this matches ``column_values`` in the upstream code); genes that
    were not part of a dataset's marker panel remain NaN.
    """

    result = np.full((len(current), len(MARKER_GENES)), np.nan, dtype=np.float32)
    audits: dict[str, Any] = {}
    for dataset in DATASETS:
        row_mask = current["dataset"].eq(dataset).to_numpy()
        panel = set(marker_column_presence.get(dataset, []))
        path = resolve_original_dataset_path(config, dataset)
        if path is None:
            audits[dataset] = {"status": "not_available", "path": "", "marker_genes": []}
            continue
        try:
            source = ad.read_h5ad(path, backed="r")
            source_index = pd.Index(source.obs_names.astype(str))
            source_ids = pd.Index(current.loc[row_mask, "source_cell_id"].astype(str))
            positions = source_index.get_indexer(source_ids)
            if (positions < 0).any():
                audits[dataset] = {
                    "status": "source_cell_id_missing",
                    "path": str(path),
                    "n_missing": int((positions < 0).sum()),
                }
                source.file.close()
                continue
            if "counts" not in source.layers:
                audits[dataset] = {
                    "status": "missing_counts_layer",
                    "path": str(path),
                    "marker_genes": [],
                }
                source.file.close()
                continue

            source_var_names = pd.Index(source.var_names.astype(str))
            var_lookup: dict[str, str] = {}
            for actual in source_var_names:
                var_lookup.setdefault(str(actual).upper(), str(actual))
            for column in ["gene_symbol", "gene_symbol_upper"]:
                if column in source.var:
                    for actual, symbol in zip(source_var_names, source.var[column].astype(str)):
                        var_lookup.setdefault(str(symbol).upper(), str(actual))
            actual_genes = [var_lookup.get(gene.upper()) for gene in MARKER_GENES]
            available_actual = [actual for actual in actual_genes if actual is not None]
            raw_values = np.zeros((int(row_mask.sum()), len(MARKER_GENES)), dtype=np.float32)
            if available_actual:
                order = np.argsort(positions)
                subset = source[positions[order], available_actual].to_memory()
                values_sorted = dense_values(subset.layers["counts"])
                values = np.empty_like(values_sorted)
                values[order] = values_sorted
                for j, actual in enumerate(available_actual):
                    gene = next(gene for gene, candidate in zip(MARKER_GENES, actual_genes) if candidate == actual)
                    raw_values[:, MARKER_GENES.index(gene)] = values[:, j]
            for j, gene in enumerate(MARKER_GENES):
                if gene not in panel:
                    raw_values[:, j] = np.nan
            result[row_mask, :] = raw_values
            audits[dataset] = {
                "status": "ok",
                "path": str(path),
                "n_cells": int(row_mask.sum()),
                "marker_panel": sorted(panel),
                "genes_present_in_source_counts": [gene for gene, actual in zip(MARKER_GENES, actual_genes) if actual is not None],
                "genes_absent_in_source_counts_treated_as_zero": [gene for gene, actual in zip(MARKER_GENES, actual_genes) if actual is None and gene in panel],
            }
            source.file.close()
        except Exception as exc:
            audits[dataset] = {
                "status": "read_error",
                "path": str(path),
                "error": f"{type(exc).__name__}: {exc}",
            }
    return result, audits


def marker_bins(values: np.ndarray) -> np.ndarray:
    bins = np.full(values.shape, "not_available", dtype=object)
    finite = np.isfinite(values)
    bins[finite & (values == 0)] = "0"
    bins[finite & (values == 1)] = "1"
    bins[finite & (values == 2)] = "2"
    bins[finite & (values >= 3)] = ">=3"
    bins[finite & (values < 0)] = "negative"
    integer_mask = finite & (values >= 0) & (np.floor(values) != values)
    bins[integer_mask] = "non_integer"
    return bins


def make_marker_pattern(values: np.ndarray) -> str:
    detected = [gene for gene, value in zip(MARKER_GENES, values) if np.isfinite(value) and value > 0]
    return " + ".join(detected) if detected else "none_detected"


def make_count_json(values: np.ndarray, only_detected: bool = False) -> str:
    payload: dict[str, float] = {}
    for gene, value in zip(MARKER_GENES, values):
        if not np.isfinite(value):
            continue
        if only_detected and value <= 0:
            continue
        payload[gene] = int(value) if float(value).is_integer() else float(value)
    return json_value(payload)


def load_source_tables(config: dict[str, Any], current: pd.DataFrame) -> tuple[pd.DataFrame, list[dict[str, Any]]]:
    records: list[dict[str, Any]] = []
    file_audit: list[dict[str, Any]] = []
    tables: dict[str, pd.DataFrame] = {}
    for dataset in DATASETS:
        path = resolve_candidate_path(config, dataset)
        if path is None:
            file_audit.append({"dataset": dataset, "path": "", "status": "missing"})
            tables[dataset] = pd.DataFrame()
            continue
        table = pd.read_csv(path, dtype=str)
        duplicate_ids = int(table["cell_id"].astype(str).duplicated(keep=False).sum()) if "cell_id" in table else -1
        if "cell_id" not in table:
            raise ValueError(f"{path} is missing cell_id")
        table["cell_id"] = table["cell_id"].astype(str)
        tables[dataset] = table.drop_duplicates("cell_id", keep="first").set_index("cell_id")
        subset = current[current["dataset"].astype(str).eq(dataset)]
        missing = int((~subset["source_cell_id"].astype(str).isin(tables[dataset].index)).sum())
        file_audit.append(
            {
                "dataset": dataset,
                "path": str(path),
                "status": "available",
                "n_rows": int(len(table)),
                "duplicate_id_rows": duplicate_ids,
                "n_consensus_cells": int(len(subset)),
                "n_consensus_ids_missing": missing,
            }
        )

    source_columns = [
        "condition",
        "donor_id",
        "assay",
        "classical_marker_count",
        "known_marker_genes_detected",
        "source_author_style",
        "source_formal_signature",
        "pool_high_confidence",
        "pool_broad_lam_like",
        "pool_unrestricted_lam",
        "candidate_reason",
    ]
    for row in current.itertuples(index=False):
        table = tables.get(str(row.dataset), pd.DataFrame())
        if not table.empty and str(row.source_cell_id) in table.index:
            source_row = table.loc[str(row.source_cell_id)]
            if isinstance(source_row, pd.DataFrame):
                source_row = source_row.iloc[0]
            record = {column: source_row.get(column, pd.NA) for column in source_columns}
        else:
            record = {column: pd.NA for column in source_columns}
        records.append(record)
    return pd.DataFrame(records), file_audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/state_modeling.yaml"))
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "results/stage15"))
    args = parser.parse_args()
    config = load_config(Path(args.config).resolve())
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    assignment_path = PROJECT_ROOT / str(config["outputs"]["step7_dir"]) / "state_consensus_with_upstream_annotations.csv"
    if not assignment_path.exists():
        raise FileNotFoundError(f"Missing Step 7 merged annotation: {assignment_path}")
    current = pd.read_csv(assignment_path, dtype=str)
    required_current = {
        "consensus_state",
        "candidate_reason",
        "source_author_style",
        "source_formal_signature",
        "known_marker_combo_ge2",
        "dataset",
        "patient_id",
        "source_cell_id",
        "analysis_cell_id",
    }
    missing_current = sorted(required_current - set(current.columns))
    if missing_current:
        raise ValueError(f"Step 7 merged annotation is missing columns: {missing_current}")
    current = current.reset_index(drop=True)
    current["dataset"] = current["dataset"].astype(str)
    current["source_cell_id"] = current["source_cell_id"].astype(str)
    current["analysis_cell_id"] = current["analysis_cell_id"].astype(str)
    current["_row"] = np.arange(len(current), dtype=int)

    prepared_path = resolve_prepared_path(config)
    prepared = ad.read_h5ad(prepared_path, backed="r")
    prepared_index = pd.Index(prepared.obs_names.astype(str))
    selected_ids = pd.Index(current["analysis_cell_id"].astype(str))
    missing_prepared = sorted(set(selected_ids) - set(prepared_index))
    if missing_prepared:
        raise RuntimeError(f"{len(missing_prepared)} consensus IDs are absent from prepared AnnData")

    source, file_audit = load_source_tables(config, current)
    source = source.add_prefix("source_")

    all_obs_columns = set(prepared.obs.columns)
    metadata_columns = [
        "source_cell_id",
        "dataset",
        "condition",
        "patient_id",
        "assay",
        "upstream_pool_high_confidence",
        "upstream_pool_broad_lam_like",
        "upstream_pool_unrestricted_lam",
        "lam_candidate",
        "boundary",
    ]
    for dataset in DATASETS:
        metadata_columns.extend(
            [
                f"upstream_{dataset}_adata_marker_expr_{gene}"
                for gene in MARKER_GENES
                if f"upstream_{dataset}_adata_marker_expr_{gene}" in all_obs_columns
            ]
        )
        metadata_columns.extend(
            [
                f"upstream_{dataset}_adata_known_marker_genes_detected",
                f"upstream_{dataset}_adata_known_marker_combo_ge2",
                f"upstream_{dataset}_candidate_classical_marker_count",
                f"upstream_{dataset}_candidate_known_marker_genes_detected",
                f"upstream_{dataset}_candidate_source_author_style",
                f"upstream_{dataset}_candidate_source_formal_signature",
                f"upstream_{dataset}_candidate_candidate_reason",
            ]
        )
    metadata_columns = list(dict.fromkeys(column for column in metadata_columns if column in all_obs_columns))
    metadata = prepared.obs.reindex(selected_ids)[metadata_columns].reset_index(drop=True)
    metadata_missing = int(metadata["source_cell_id"].isna().sum()) if "source_cell_id" in metadata else len(metadata)

    # Marker values are the pre-normalization upstream marker columns.  Core
    # GSE135851 has seven columns; external datasets have all eight.
    marker_values = np.full((len(current), len(MARKER_GENES)), np.nan, dtype=np.float32)
    marker_column_presence: dict[str, list[str]] = {}
    for dataset in DATASETS:
        row_mask = current["dataset"].eq(dataset).to_numpy()
        present: list[str] = []
        for j, gene in enumerate(MARKER_GENES):
            column = f"upstream_{dataset}_adata_marker_expr_{gene}"
            if column not in metadata:
                continue
            present.append(gene)
            values = pd.to_numeric(metadata[column], errors="coerce").to_numpy(dtype=float)
            marker_values[row_mask, j] = np.nan_to_num(values[row_mask], nan=0.0)
        marker_column_presence[dataset] = present

    detected = np.isfinite(marker_values) & (marker_values > 0)
    recomputed_marker_count = detected.sum(axis=1).astype(int)
    recomputed_combo = recomputed_marker_count >= 2
    marker_pattern = np.array([make_marker_pattern(row) for row in marker_values], dtype=object)
    marker_counts_json = np.array([make_count_json(row) for row in marker_values], dtype=object)

    # The current, alias-canonicalized counts are read separately.  This is a
    # small 5,378 x marker subset, not the full 131k x 20k matrix.
    canonical_counts, available_count_genes, counts_audit = read_selected_counts(
        prepared, selected_ids, MARKER_GENES
    )
    # For UMI-level diagnostics, use the original per-dataset AnnData files
    # from data-temp/source data.  These retain FIGF and VEGFD before the
    # State Modeling FIGF -> VEGFD merge.
    original_raw_counts, original_counts_audit = read_original_raw_counts(
        config, current, marker_column_presence
    )
    raw_detected = np.isfinite(original_raw_counts) & (original_raw_counts > 0)
    upstream_detected = detected.copy()
    raw_detection_mismatch = np.zeros(len(current), dtype=bool)
    for dataset in DATASETS:
        row_mask = current["dataset"].eq(dataset).to_numpy()
        panel = marker_column_presence.get(dataset, [])
        for gene in panel:
            j = MARKER_GENES.index(gene)
            raw_detection_mismatch[row_mask] |= raw_detected[row_mask, j] != upstream_detected[row_mask, j]
    raw_support_complete = np.all(np.isfinite(original_raw_counts) | ~detected, axis=1)
    observed_raw = original_raw_counts[np.isfinite(original_raw_counts)]
    raw_count_validity = {
        "n_marker_values": int(observed_raw.size),
        "finite": bool(np.isfinite(observed_raw).all()),
        "nonnegative": bool((observed_raw >= 0).all()) if observed_raw.size else True,
        "integer_valued": bool(np.isclose(observed_raw, np.round(observed_raw), rtol=0, atol=1e-6).all()) if observed_raw.size else True,
        "n_non_integer": int(np.count_nonzero(~np.isclose(observed_raw, np.round(observed_raw), rtol=0, atol=1e-6))),
        "n_negative": int(np.count_nonzero(observed_raw < 0)),
    }

    audit = pd.concat(
        [current, source, metadata.add_prefix("prepared_")], axis=1
    ).copy()

    audit["source_condition"] = audit["source_condition"].fillna(audit["prepared_condition"])
    audit["condition_for_recalculation"] = audit["source_condition"].fillna(current.get("condition", ""))
    audit["source_author_style_bool"] = bool_series(audit["source_source_author_style"])
    audit["source_formal_signature_bool"] = bool_series(audit["source_source_formal_signature"])
    audit["source_pool_high_confidence_bool"] = bool_series(audit["source_pool_high_confidence"])
    audit["prepared_pool_high_confidence_bool"] = bool_series(audit["prepared_upstream_pool_high_confidence"])
    # Use the dataset-specific metadata field rather than the first dataset's
    # field; the first field is only populated for GSE135851.
    prepared_combo = np.zeros(len(audit), dtype=bool)
    prepared_known_detected = np.full(len(audit), np.nan, dtype=float)
    prepared_candidate_count = np.full(len(audit), np.nan, dtype=float)
    prepared_candidate_known = np.full(len(audit), np.nan, dtype=float)
    prepared_author = np.zeros(len(audit), dtype=bool)
    prepared_formal = np.zeros(len(audit), dtype=bool)
    prepared_reason = np.full(len(audit), "", dtype=object)
    for dataset in DATASETS:
        row_mask = audit["dataset"].eq(dataset).to_numpy()
        combo_col = f"prepared_upstream_{dataset}_adata_known_marker_combo_ge2"
        known_col = f"prepared_upstream_{dataset}_adata_known_marker_genes_detected"
        count_col = f"prepared_upstream_{dataset}_candidate_classical_marker_count"
        candidate_known_col = f"prepared_upstream_{dataset}_candidate_known_marker_genes_detected"
        author_col = f"prepared_upstream_{dataset}_candidate_source_author_style"
        formal_col = f"prepared_upstream_{dataset}_candidate_source_formal_signature"
        reason_col = f"prepared_upstream_{dataset}_candidate_candidate_reason"
        if combo_col in audit:
            prepared_combo[row_mask] = bool_series(audit.loc[row_mask, combo_col]).to_numpy()
        if known_col in audit:
            prepared_known_detected[row_mask] = numeric_series(audit.loc[row_mask, known_col]).to_numpy()
        if count_col in audit:
            prepared_candidate_count[row_mask] = numeric_series(audit.loc[row_mask, count_col]).to_numpy()
        if candidate_known_col in audit:
            prepared_candidate_known[row_mask] = numeric_series(audit.loc[row_mask, candidate_known_col]).to_numpy()
        if author_col in audit:
            prepared_author[row_mask] = bool_series(audit.loc[row_mask, author_col]).to_numpy()
        if formal_col in audit:
            prepared_formal[row_mask] = bool_series(audit.loc[row_mask, formal_col]).to_numpy()
        if reason_col in audit:
            prepared_reason[row_mask] = (
                audit.loc[row_mask, reason_col].astype(object).fillna("").astype(str).to_numpy()
            )
    audit["prepared_adata_known_marker_combo_ge2_bool"] = prepared_combo
    audit["prepared_adata_known_marker_genes_detected"] = prepared_known_detected
    audit["prepared_candidate_classical_marker_count"] = prepared_candidate_count
    audit["prepared_candidate_known_marker_genes_detected"] = prepared_candidate_known
    audit["prepared_candidate_source_author_style_bool"] = prepared_author
    audit["prepared_candidate_source_formal_signature_bool"] = prepared_formal
    audit["prepared_candidate_reason"] = prepared_reason

    source_known = numeric_series(audit["source_known_marker_genes_detected"])
    source_marker_count = numeric_series(audit["source_classical_marker_count"])
    condition_lam = audit["condition_for_recalculation"].astype(str).str.upper().eq("LAM").to_numpy()
    source_reason = audit["source_candidate_reason"].fillna("").astype(str).to_numpy()
    source_known_combo = source_known.to_numpy() >= 2
    author = audit["source_author_style_bool"].to_numpy()
    formal = audit["source_formal_signature_bool"].to_numpy()
    recomputed_high = condition_lam & (author | formal | recomputed_combo)
    recomputed_reason = np.select(
        [
            author | formal,
            recomputed_combo,
            recomputed_marker_count >= 1,
            source_known.to_numpy() >= 1,
            condition_lam,
        ],
        [
            "author_or_formal_candidate",
            "known_marker_combo_and_marker_support",
            "weak_classical_marker",
            "known_marker_detection",
            "LAM_condition_guardrail",
        ],
        default="none",
    )
    audit["marker_count_recomputed"] = recomputed_marker_count
    audit["known_marker_combo_ge2_recomputed"] = recomputed_combo
    audit["known_marker_genes_detected_recomputed"] = recomputed_marker_count
    audit["pool_high_confidence_recomputed"] = recomputed_high
    audit["candidate_reason_recomputed"] = recomputed_reason
    audit["marker_pattern"] = marker_pattern
    audit["marker_pair"] = np.array(
        [pattern if int(count) == 2 else "" for pattern, count in zip(marker_pattern, recomputed_marker_count)],
        dtype=object,
    )
    # Replace the provisional expression-value labels with explicit names:
    # marker_expression_values_json is the upstream value used for >0
    # detection, while marker_counts_json is the original raw-count value.
    audit["marker_expression_values_json"] = marker_counts_json
    audit["marker_counts_json"] = np.array(
        [make_count_json(row) for row in original_raw_counts], dtype=object
    )
    audit["supporting_marker_counts_json"] = np.array(
        [make_count_json(row, only_detected=True) for row in original_raw_counts], dtype=object
    )
    only_one_umi = np.array(
        [
            bool(
                source_reason[i] == "known_marker_combo_and_marker_support"
            )
            and bool(recomputed_combo[i])
            and bool(raw_support_complete[i])
            and bool(raw_detected[i].any())
            and bool(np.all(original_raw_counts[i, detected[i]] == 1))
            for i in range(len(audit))
        ],
        dtype=bool,
    )
    audit["only_1_umi_detections"] = only_one_umi
    audit["raw_support_counts_complete"] = raw_support_complete
    audit["raw_detection_vs_upstream_marker_expr_mismatch"] = raw_detection_mismatch
    audit["source_combo_from_known_detected"] = source_known_combo
    audit["prepared_combo_matches_recomputed"] = prepared_combo == recomputed_combo

    # Alias audit.  The original count is the number of positive upstream
    # marker columns.  Collapsing FIGF and VEGFD removes one count only when
    # both are positive in the same cell.
    figf_index = MARKER_GENES.index("FIGF")
    vegfd_index = MARKER_GENES.index("VEGFD")
    figf_positive = detected[:, figf_index]
    vegfd_positive = detected[:, vegfd_index]
    duplicate_pass = figf_positive & vegfd_positive
    alias_corrected_count = recomputed_marker_count - duplicate_pass.astype(int)
    alias_corrected_combo = alias_corrected_count >= 2
    alias_drop = recomputed_combo & ~alias_corrected_combo
    audit["FIGF_positive"] = figf_positive
    audit["VEGFD_positive"] = vegfd_positive
    audit["FIGF_VEGFD_duplicate_pass"] = duplicate_pass
    audit["marker_count_alias_corrected"] = alias_corrected_count
    audit["known_marker_combo_ge2_alias_corrected"] = alias_corrected_combo
    audit["original_ge2_but_alias_corrected_lt2"] = alias_drop

    # Field-level identity/merge checks.  The candidate source has no explicit
    # combo column, so combo is checked against both the inherited AnnData
    # value and an independent recomputation from its original marker columns.
    merged_author = bool_series(current["source_author_style"])
    merged_formal = bool_series(current["source_formal_signature"])
    merged_combo = bool_series(current["known_marker_combo_ge2"])
    merged_reason = current["candidate_reason"].fillna("").astype(str)
    audit["mismatch_pool_high_confidence_source_vs_prepared"] = (
        audit["source_pool_high_confidence_bool"].to_numpy() != audit["prepared_pool_high_confidence_bool"].to_numpy()
    )
    audit["mismatch_candidate_reason_source_vs_merged"] = source_reason != merged_reason.to_numpy()
    audit["mismatch_author_source_vs_merged"] = author != merged_author.to_numpy()
    audit["mismatch_formal_source_vs_merged"] = formal != merged_formal.to_numpy()
    audit["mismatch_combo_merged_vs_prepared"] = merged_combo.to_numpy() != prepared_combo
    audit["mismatch_combo_merged_vs_recomputed"] = merged_combo.to_numpy() != recomputed_combo
    audit["mismatch_marker_count_source_vs_recomputed"] = (
        source_marker_count.to_numpy() != recomputed_marker_count
    )
    audit["mismatch_known_detected_source_vs_recomputed"] = source_known.to_numpy() != recomputed_marker_count
    audit["mismatch_source_vs_recomputed_reason"] = source_reason != recomputed_reason
    mismatch_columns = [column for column in audit.columns if column.startswith("mismatch_")]
    audit["field_inconsistency"] = audit[mismatch_columns].any(axis=1)

    # Canonical count JSON is deliberately separate from the upstream raw
    # marker JSON: VEGFD here is the post-alias count and FIGF is unavailable.
    audit["canonical_counts_json"] = [make_count_json(row) for row in canonical_counts]
    audit["canonical_count_genes_available"] = json_value(available_count_genes)

    cell_columns = [
        "analysis_cell_id",
        "source_cell_id",
        "consensus_state",
        "dataset",
        "patient_id",
        "source_condition",
        "candidate_reason",
        "source_candidate_reason",
        "source_author_style",
        "source_source_author_style",
        "source_formal_signature",
        "source_source_formal_signature",
        "source_pool_high_confidence",
        "source_classical_marker_count",
        "source_known_marker_genes_detected",
        "known_marker_combo_ge2",
        "marker_count_recomputed",
        "known_marker_combo_ge2_recomputed",
        "pool_high_confidence_recomputed",
        "candidate_reason_recomputed",
        "marker_pattern",
        "marker_pair",
        "marker_expression_values_json",
        "marker_counts_json",
        "supporting_marker_counts_json",
        "only_1_umi_detections",
        "raw_support_counts_complete",
        "raw_detection_vs_upstream_marker_expr_mismatch",
        "FIGF_positive",
        "VEGFD_positive",
        "FIGF_VEGFD_duplicate_pass",
        "marker_count_alias_corrected",
        "original_ge2_but_alias_corrected_lt2",
        "canonical_counts_json",
        "field_inconsistency",
    ]
    audit.loc[state_sort(audit).index, cell_columns].to_csv(
        output_dir / "candidate_identity_cell_audit.csv", index=False
    )

    # 1. Annotation/ID merge audit summary.
    merge_rows = [
        {"check": "consensus_cell_count", "value": int(len(current)), "status": "ok"},
        {"check": "consensus_duplicate_analysis_cell_id", "value": int(current["analysis_cell_id"].duplicated().sum()), "status": "ok" if current["analysis_cell_id"].duplicated().sum() == 0 else "fail"},
        {"check": "consensus_duplicate_source_cell_id_within_dataset", "value": int(current.duplicated(["dataset", "source_cell_id"]).sum()), "status": "ok" if current.duplicated(["dataset", "source_cell_id"]).sum() == 0 else "fail"},
        {"check": "missing_prepared_annotation_rows", "value": metadata_missing, "status": "ok" if metadata_missing == 0 else "fail"},
        {"check": "missing_source_candidate_ids", "value": int(source["source_candidate_reason"].isna().sum()), "status": "ok" if source["source_candidate_reason"].isna().sum() == 0 else "fail"},
        {"check": "field_inconsistency_cells", "value": int(audit["field_inconsistency"].sum()), "status": "ok" if audit["field_inconsistency"].sum() == 0 else "fail"},
    ]
    for column in mismatch_columns:
        merge_rows.append(
            {
                "check": column,
                "value": int(audit[column].sum()),
                "status": "ok" if int(audit[column].sum()) == 0 else "fail",
            }
        )
    pd.DataFrame(merge_rows).to_csv(output_dir / "annotation_merge_audit.csv", index=False)

    # 2. Recalculation summary.
    recalculation_rows = []
    recalculation_checks = {
        "source_classical_marker_count_vs_recomputed": "mismatch_marker_count_source_vs_recomputed",
        "source_known_marker_genes_detected_vs_recomputed": "mismatch_known_detected_source_vs_recomputed",
        "source_pool_high_confidence_vs_recomputed": "source_pool_high_confidence_bool",
        "source_candidate_reason_vs_recomputed": "mismatch_source_vs_recomputed_reason",
        "prepared_adata_combo_vs_recomputed": "mismatch_combo_merged_vs_recomputed",
    }
    for check, mismatch_column in recalculation_checks.items():
        if check == "source_pool_high_confidence_vs_recomputed":
            mismatch = audit["source_pool_high_confidence_bool"].to_numpy() != recomputed_high
        else:
            mismatch = audit[mismatch_column].to_numpy()
        recalculation_rows.append(
            {
                "check": check,
                "n_cells": int(len(audit)),
                "n_mismatch": int(mismatch.sum()),
                "fraction_mismatch": float(mismatch.mean()) if len(mismatch) else 0.0,
                "status": "ok" if not mismatch.any() else "fail",
            }
        )
    pd.DataFrame(recalculation_rows).to_csv(output_dir / "rule_recalculation_audit.csv", index=False)

    # 3. Alias audit by dataset and state.
    audit["marker_combo_candidate"] = source_reason == "known_marker_combo_and_marker_support"

    def alias_summary(group_keys: list[str]) -> pd.DataFrame:
        rows: list[dict[str, Any]] = []
        for key, frame in audit.groupby(group_keys, observed=True, dropna=False):
            if not isinstance(key, tuple):
                key = (key,)
            row = dict(zip(group_keys, [str(value) for value in key]))
            combo_n = int(frame["marker_combo_candidate"].sum())
            row.update(
                {
                    "cells": int(len(frame)),
                    "FIGF_positive": int(frame["FIGF_positive"].sum()),
                    "VEGFD_positive": int(frame["VEGFD_positive"].sum()),
                    "FIGF_and_VEGFD_positive": int(frame["FIGF_VEGFD_duplicate_pass"].sum()),
                    "duplicate_pass_fraction": float(frame["FIGF_VEGFD_duplicate_pass"].mean()),
                    "original_marker_count_ge2": int(frame["known_marker_combo_ge2_recomputed"].sum()),
                    "alias_corrected_marker_count_ge2": int(frame["known_marker_combo_ge2_alias_corrected"].sum()),
                    "original_ge2_but_alias_corrected_lt2": int(frame["original_ge2_but_alias_corrected_lt2"].sum()),
                    "alias_drop_fraction_among_original_ge2": float(
                        frame["original_ge2_but_alias_corrected_lt2"].sum()
                        / frame["known_marker_combo_ge2_recomputed"].sum()
                    )
                    if frame["known_marker_combo_ge2_recomputed"].sum()
                    else 0.0,
                    "duplicate_pass_fraction_among_marker_combo_candidates": float(
                        frame.loc[frame["marker_combo_candidate"], "FIGF_VEGFD_duplicate_pass"].mean()
                    )
                    if combo_n
                    else 0.0,
                    "marker_combo_candidates": combo_n,
                }
            )
            rows.append(row)
        return state_sort(pd.DataFrame(rows), group_keys[-1]) if "consensus_state" in group_keys else pd.DataFrame(rows)

    alias_summary(["dataset"]).to_csv(output_dir / "alias_audit_by_dataset.csv", index=False)
    alias_summary(["consensus_state"]).to_csv(output_dir / "alias_audit_by_state.csv", index=False)

    # 4. Marker combinations/patterns.
    pattern_rows: list[dict[str, Any]] = []
    for pattern, frame in audit.groupby("marker_pattern", observed=True, dropna=False):
        pattern_rows.append(
            {
                "marker_pattern": str(pattern),
                "marker_pair": str(pattern) if int(frame["marker_count_recomputed"].iloc[0]) == 2 else "",
                "detected_marker_count": int(frame["marker_count_recomputed"].iloc[0]),
                "cells": int(len(frame)),
                "marker_combo_candidate_cells": int(frame["marker_combo_candidate"].sum()),
                "states": int(frame["consensus_state"].nunique()),
                "datasets": int(frame["dataset"].nunique()),
                "patients": int(frame["patient_id"].nunique()),
                "state_counts": compact_counts(frame["consensus_state"]),
                "dataset_counts": compact_counts(frame["dataset"]),
                "patient_counts": compact_counts(frame["patient_id"]),
            }
        )
    pattern_table = pd.DataFrame(pattern_rows).sort_values(["marker_combo_candidate_cells", "cells"], ascending=False)
    pattern_table.to_csv(output_dir / "marker_patterns.csv", index=False)

    # 5. Raw marker UMI distributions, globally and by state.
    umi_rows: list[dict[str, Any]] = []
    state_umi_rows: list[dict[str, Any]] = []
    for j, gene in enumerate(MARKER_GENES):
        bins = marker_bins(original_raw_counts[:, j])
        for bin_name, n in Counter(bins.tolist()).items():
            available_n = int(np.isfinite(original_raw_counts[:, j]).sum())
            umi_rows.append(
                {
                    "marker": gene,
                    "count_bin": str(bin_name),
                    "cells": int(n),
                    "available_cells": available_n,
                    "fraction_of_available_cells": float(n / available_n) if available_n else 0.0,
                    "value_source": "original_dataset_layers_counts",
                }
            )
        for state, frame in audit.groupby("consensus_state", observed=True):
            positions = frame.index.to_numpy(dtype=int)
            state_bins = marker_bins(original_raw_counts[positions, j])
            available_n = int(np.isfinite(original_raw_counts[positions, j]).sum())
            for bin_name, n in Counter(state_bins.tolist()).items():
                state_umi_rows.append(
                    {
                        "consensus_state": str(state),
                        "marker": gene,
                        "count_bin": str(bin_name),
                        "cells": int(n),
                        "available_cells": available_n,
                        "fraction_of_available_cells": float(n / available_n) if available_n else 0.0,
                        "value_source": "original_dataset_layers_counts",
                    }
                )
    pd.DataFrame(umi_rows).to_csv(output_dir / "marker_umi_distribution.csv", index=False)
    state_umi = pd.DataFrame(state_umi_rows)
    if not state_umi.empty:
        state_umi = state_sort(state_umi, "consensus_state")
    state_umi.to_csv(output_dir / "state_marker_umi_distribution.csv", index=False)

    # 6. State-level diagnostic container.
    state_rows: list[dict[str, Any]] = []
    for state, frame in audit.groupby("consensus_state", observed=True, dropna=False):
        combo_frame = frame[frame["marker_combo_candidate"]]
        pattern_counts = frame["marker_pattern"].value_counts().head(5).to_dict()
        duplicate_n = int(frame["FIGF_VEGFD_duplicate_pass"].sum())
        combo_n = int(len(combo_frame))
        low_umi_combo_n = int(combo_frame["only_1_umi_detections"].sum())
        state_rows.append(
            {
                "consensus_state": str(state),
                "cells": int(len(frame)),
                "author_style_cells": int(frame["source_author_style_bool"].sum()),
                "formal_signature_cells": int(frame["source_formal_signature_bool"].sum()),
                "author_or_formal_fraction": float(
                    (frame["source_author_style_bool"] | frame["source_formal_signature_bool"]).mean()
                ),
                "marker_combo_candidate_cells": combo_n,
                "marker_combo_fraction": float(combo_n / len(frame)) if len(frame) else 0.0,
                "dominant_marker_combinations": json_value(pattern_counts),
                "FIGF_VEGFD_duplicate_pass_cells": duplicate_n,
                "FIGF_VEGFD_duplicate_pass_fraction": float(frame["FIGF_VEGFD_duplicate_pass"].mean()),
                "FIGF_VEGFD_duplicate_pass_fraction_among_marker_combo": float(
                    combo_frame["FIGF_VEGFD_duplicate_pass"].mean()
                )
                if combo_n
                else 0.0,
                "median_detected_marker_count": float(frame["marker_count_recomputed"].median()),
                "fraction_supported_only_by_1_UMI_detections": float(low_umi_combo_n / combo_n) if combo_n else 0.0,
                "only_1_UMI_marker_combo_cells": low_umi_combo_n,
                "original_ge2_but_alias_corrected_lt2_cells": int(frame["original_ge2_but_alias_corrected_lt2"].sum()),
                "datasets": int(frame["dataset"].nunique()),
                "patients": int(frame["patient_id"].nunique()),
            }
        )
    state_summary = state_sort(pd.DataFrame(state_rows))
    state_summary.to_csv(output_dir / "state_identity_summary.csv", index=False)

    # Root-cause evidence is deliberately quantitative and does not change
    # any candidate rule.
    n_combo = int(audit["marker_combo_candidate"].sum())
    n_author_formal = int((author | formal).sum())
    n_alias_drop = int(alias_drop.sum())
    n_only_one = int(audit.loc[audit["marker_combo_candidate"], "only_1_umi_detections"].sum())
    n_combo_with_nonspecific = int(
        audit.loc[audit["marker_combo_candidate"], "marker_pattern"].map(
            lambda pattern: sum(gene in str(pattern).split(" + ") for gene in ["ACTA2", "CTSK", "VEGFD"]) >= 2
        ).sum()
    )
    annotation_pass = int(audit["field_inconsistency"].sum()) == 0 and not any(
        int(item["value"]) > 0
        for item in merge_rows
        if item["check"] in {"missing_prepared_annotation_rows", "missing_source_candidate_ids"}
    )
    root_cause_rows = [
        {
            "root_cause": "A_annotation_or_ID_merge_error",
            "evidence_metric": "field_inconsistency_cells",
            "value": int(audit["field_inconsistency"].sum()),
            "fraction_of_5378": float(audit["field_inconsistency"].mean()),
            "interpretation": "not_supported_by_this_audit" if annotation_pass else "requires_investigation",
        },
        {
            "root_cause": "B_FIGF_VEGFD_alias_duplicate_counting",
            "evidence_metric": "original_ge2_but_alias_corrected_lt2",
            "value": n_alias_drop,
            "fraction_of_marker_combo_candidates": float(n_alias_drop / n_combo) if n_combo else 0.0,
            "interpretation": "quantified_contribution" if n_alias_drop else "no_cells_observed",
        },
        {
            "root_cause": "C_marker_combo_gate_specificity",
            "evidence_metric": "marker_combo_candidates_with_only_1_UMI_support",
            "value": n_only_one,
            "fraction_of_marker_combo_candidates": float(n_only_one / n_combo) if n_combo else 0.0,
            "interpretation": "rule_permissiveness_evidence_to_review" if n_only_one else "not_observed",
        },
        {
            "root_cause": "C_marker_combo_gate_specificity",
            "evidence_metric": "marker_combo_candidates_with_at_least_two_of_ACTA2_CTSK_VEGFD",
            "value": n_combo_with_nonspecific,
            "fraction_of_marker_combo_candidates": float(n_combo_with_nonspecific / n_combo) if n_combo else 0.0,
            "interpretation": "pattern_reported_for_biological_review",
        },
        {
            "root_cause": "context",
            "evidence_metric": "author_or_formal_supported_cells",
            "value": n_author_formal,
            "fraction_of_5378": float(n_author_formal / len(audit)) if len(audit) else 0.0,
            "interpretation": "upstream_source_flag_support",
        },
    ]
    pd.DataFrame(root_cause_rows).to_csv(output_dir / "root_cause_evidence.csv", index=False)

    manifest = {
        "script": str(Path(__file__).resolve()),
        "read_only_existing_artifacts": True,
        "n_cells": int(len(audit)),
        "n_states": int(audit["consensus_state"].nunique()),
        "state_consensus_input": str(assignment_path),
        "prepared_h5ad": str(prepared_path),
        "candidate_files": file_audit,
        "marker_genes": MARKER_GENES,
        "marker_columns_present_by_dataset": marker_column_presence,
        "counts_audit": counts_audit,
        "original_counts_audit": original_counts_audit,
        "original_raw_marker_count_validity": raw_count_validity,
        "counts_layer_provenance": "current canonicalized prepared AnnData; FIGF is absent and VEGFD is alias-merged",
        "marker_expr_provenance": "upstream marker_expr columns retained in prepared AnnData; used only for exact >0 rule reproduction",
        "raw_marker_count_provenance": "original per-dataset AnnData layers[counts] from resolved source/data-temp files",
        "annotation_merge_pass": bool(annotation_pass),
        "n_marker_combo_candidates": n_combo,
        "n_author_or_formal_supported": n_author_formal,
        "n_alias_duplicate_pass": int(duplicate_pass.sum()),
        "n_original_ge2_but_alias_corrected_lt2": n_alias_drop,
        "n_marker_combo_only_1_umi": n_only_one,
        "outputs": sorted(path.name for path in output_dir.iterdir() if path.is_file()),
    }
    with (output_dir / "candidate_identity_audit.json").open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, ensure_ascii=False, indent=2, default=str)

    report_lines = [
        "# Candidate identity audit",
        "",
        f"- Audited cells: {len(audit):,}",
        f"- Consensus states: {audit['consensus_state'].nunique()}",
        f"- Marker-combo candidate reason: {n_combo:,} cells",
        f"- Author/formal support: {n_author_formal:,} cells",
        f"- Formal signature support alone: {int(formal.sum()):,} cells",
        f"- Annotation/ID merge field inconsistencies: {int(audit['field_inconsistency'].sum()):,}",
        f"- FIGF+VEGFD duplicate passes: {int(duplicate_pass.sum()):,}",
        f"- Original >=2 but alias-corrected <2: {n_alias_drop:,}",
        f"- Marker-combo candidates supported only by 1-UMI detections: {n_only_one:,}",
        "",
        "The audit does not modify candidate rules or existing Stage 1–13 artifacts.",
        "FIGF/VEGFD detections for rule reproduction come from retained upstream marker_expr columns; UMI bins come from original per-dataset layers[counts].",
        "",
        "Root-cause evidence is in root_cause_evidence.csv; state-level diagnostics are in state_identity_summary.csv.",
    ]
    (output_dir / "candidate_identity_audit.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    prepared.file.close()
    print(f"Audited {len(audit):,} cells across {audit['consensus_state'].nunique()} states")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
