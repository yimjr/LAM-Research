"""Shared data and effect-size utilities for the LAM research project."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]

# Shared project paths.  Keeping candidate-generation outputs separate from
# downstream candidate-analysis data prevents later reruns from recreating a
# mixed results/candidates directory.
CANDIDATE_RESULTS = ROOT / "results" / "candidates"
CANDIDATE_ANALYSIS_DATA = ROOT / "data" / "processed" / "candidate_analysis"
CANDIDATE_PROGRAMS = CANDIDATE_ANALYSIS_DATA / "programs"
CANDIDATE_DECONVOLUTION = CANDIDATE_ANALYSIS_DATA / "deconvolution"
CANDIDATE_DRUG_TARGETS = CANDIDATE_ANALYSIS_DATA / "drug_targets"
CANDIDATE_VALIDATION = CANDIDATE_ANALYSIS_DATA / "validation"
CANDIDATE_AUDIT = CANDIDATE_ANALYSIS_DATA / "audit"
LINCS_CANDIDATE_INTERMEDIATE = ROOT / "data" / "processed" / "LINCS" / "candidate_intermediate"
CANDIDATE_ANALYSIS_REPORTS = ROOT / "reports" / "06_candidate_analysis"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n")


def normalize_gene_index(index: Iterable[object]) -> pd.Index:
    values = pd.Index([str(x).strip() for x in index])
    return values.str.replace(r"\\.\\d+$", "", regex=True)


def collapse_duplicate_genes(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()
    table.index = normalize_gene_index(table.index)
    table = table.apply(pd.to_numeric, errors="coerce")
    table = table.groupby(level=0, sort=False).mean()
    return table


def signed_ratio(d1: pd.Series, d0: pd.Series, epsilon: float) -> pd.Series:
    denominator = d0.where(d0.abs() >= epsilon)
    return d1.divide(denominator)


def classify_ratio(ratio: pd.Series) -> pd.Series:
    result = pd.Series("unclassified", index=ratio.index, dtype="object")
    result.loc[ratio.lt(0)] = "direction_reversal"
    result.loc[ratio.ge(0) & ratio.lt(0.2)] = "near_complete_rescue"
    result.loc[ratio.ge(0.2) & ratio.lt(0.8)] = "partial_rescue_residual"
    result.loc[ratio.ge(0.8) & ratio.le(1.2)] = "persistent_residual"
    result.loc[ratio.gt(1.2)] = "worsened_residual"
    result.loc[ratio.isna()] = "unclassified"
    return result


def mean_by_group(table: pd.DataFrame, samples: pd.DataFrame, group_cols: list[str]) -> pd.DataFrame:
    aligned = samples.set_index("sample_id").loc[table.columns]
    grouped = table.T.groupby([aligned[col] for col in group_cols], sort=False).mean().T
    grouped.columns = ["|".join(map(str, col if isinstance(col, tuple) else (col,))) for col in grouped.columns]
    return grouped
