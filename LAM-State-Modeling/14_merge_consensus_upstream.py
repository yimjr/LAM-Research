#!/usr/bin/env python3
"""Merge Step 7 states with original upstream cell annotations and summarize."""

from __future__ import annotations

import argparse
import json
import os

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd

from state_modeling_utils import PROJECT_ROOT, annotation_directory, as_bool, load_config, write_json


REQUIRED_COLUMNS = [
    "consensus_state",
    "cell_type",
    "candidate_reason",
    "source_author_style",
    "source_formal_signature",
    "known_marker_combo_ge2",
    "doublet_score",
    "doublet_predicted",
    "dataset",
    "patient_id",
]


def compact_counts(values: pd.Series) -> str:
    counts = values.fillna("NA").astype(str).value_counts(dropna=False).to_dict()
    return json.dumps({str(key): int(value) for key, value in counts.items()}, ensure_ascii=False, sort_keys=True)


def candidate_path(config: dict, dataset: str):
    directory = annotation_directory(config, dataset)
    if directory is None:
        return None
    path = directory / "candidate_pool_labels.csv"
    return path if path.exists() else None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    stage7_dir = PROJECT_ROOT / config["outputs"]["step7_dir"]
    assignment_path = stage7_dir / "state_consensus_assignments.csv"
    scvi_path = PROJECT_ROOT / config["outputs"]["scvi_h5ad"]
    assignments = pd.read_csv(assignment_path, dtype=str)
    scvi = ad.read_h5ad(scvi_path, backed="r")
    scvi_index = pd.Index(scvi.obs.index.astype(str))
    selected = scvi.obs.reindex(assignments["scvi_obs_index"].astype(str)).copy()
    if selected.isna().all(axis=1).any():
        raise RuntimeError("Could not map consensus assignments to scVI AnnData metadata")
    selected.index = assignments.index
    merged = assignments.copy()
    merged["source_cell_id"] = selected["source_cell_id"].astype(str).to_numpy()
    for column in ["cell_type", "doublet_score", "doublet_predicted", "dataset", "patient_id"]:
        if column in selected:
            merged[column] = selected[column].to_numpy()
    merged["cell_type"] = merged["cell_type"].fillna("unknown").astype(str)
    merged["dataset"] = merged["dataset"].fillna(assignments["dataset"]).astype(str)
    merged["patient_id"] = merged["patient_id"].fillna(assignments["patient_id"]).astype(str)
    merged["doublet_score"] = pd.to_numeric(merged["doublet_score"], errors="coerce")
    merged["doublet_predicted"] = merged["doublet_predicted"].map(as_bool).astype(bool)

    audit_rows = []
    for dataset in sorted(merged["dataset"].unique()):
        path = candidate_path(config, dataset)
        subset = merged["dataset"].eq(dataset)
        if path is None:
            merged.loc[subset, "candidate_reason"] = ""
            merged.loc[subset, "source_author_style"] = ""
            merged.loc[subset, "source_formal_signature"] = ""
            merged.loc[subset, "known_marker_combo_ge2"] = False
            audit_rows.append({"dataset": dataset, "path": "", "n_cells": int(subset.sum()), "n_candidate_matched": 0, "status": "not_available"})
            continue
        candidate = pd.read_csv(path, dtype={"cell_id": str}).drop_duplicates("cell_id", keep="last").set_index("cell_id")
        source_ids = merged.loc[subset, "source_cell_id"].astype(str)
        aligned = candidate.reindex(source_ids.to_numpy())
        aligned.index = merged.index[subset]
        for column in ["candidate_reason", "source_author_style", "source_formal_signature"]:
            merged.loc[subset, column] = aligned[column].fillna("").astype(str).to_numpy() if column in aligned else ""
        combo_column = f"upstream_{dataset}_adata_known_marker_combo_ge2"
        if combo_column in selected:
            merged.loc[subset, "known_marker_combo_ge2"] = selected.loc[subset, combo_column].map(as_bool).to_numpy()
        elif "known_marker_combo_ge2" in aligned:
            merged.loc[subset, "known_marker_combo_ge2"] = aligned["known_marker_combo_ge2"].map(as_bool).fillna(False).to_numpy()
        else:
            merged.loc[subset, "known_marker_combo_ge2"] = False
        audit_rows.append({
            "dataset": dataset,
            "path": str(path),
            "n_cells": int(subset.sum()),
            "n_candidate_rows": int(len(candidate)),
            "n_candidate_matched": int(source_ids.isin(candidate.index).sum()),
            "status": "available",
        })

    for column in REQUIRED_COLUMNS:
        if column not in merged:
            merged[column] = "" if column not in {"known_marker_combo_ge2", "doublet_predicted"} else False
    merged["known_marker_combo_ge2"] = merged["known_marker_combo_ge2"].map(as_bool).astype(bool)
    merged["candidate_annotation_matched"] = (
        merged["candidate_reason"].astype(str).ne("")
        | merged["source_author_style"].astype(str).ne("")
        | merged["source_formal_signature"].astype(str).ne("")
    )
    ordered = REQUIRED_COLUMNS + [
        "source_cell_id",
        "scvi_obs_index",
        "analysis_cell_id",
        "assay",
        "candidate_annotation_matched",
        "upstream_pool_high_confidence",
        "upstream_pool_broad_lam_like",
        "upstream_pool_unrestricted_lam",
        "within_coassignment",
        "margin",
        "edge_score",
    ]
    ordered = [column for column in ordered if column in merged.columns]
    merged[ordered].to_csv(PROJECT_ROOT / config["outputs"]["consensus_upstream_annotations"], index=False)

    state_rows = []
    for state, frame in merged.groupby("consensus_state", observed=True):
        row = {
            "consensus_state": str(state),
            "cells": int(len(frame)),
            "patients": int(frame["patient_id"].nunique()),
            "datasets": int(frame["dataset"].nunique()),
            "cell_type_counts": compact_counts(frame["cell_type"]),
            "candidate_reason_counts": compact_counts(frame["candidate_reason"]),
            "source_author_style_counts": compact_counts(frame["source_author_style"]),
            "source_formal_signature_counts": compact_counts(frame["source_formal_signature"]),
            "known_marker_combo_ge2_cells": int(frame["known_marker_combo_ge2"].sum()),
            "doublet_predicted_cells": int(frame["doublet_predicted"].sum()),
            "mean_doublet_score": float(pd.to_numeric(frame["doublet_score"], errors="coerce").mean()),
            "candidate_annotation_matched_cells": int(frame["candidate_annotation_matched"].sum()),
        }
        state_rows.append(row)
    summary = pd.DataFrame(state_rows)
    summary["_state_sort"] = pd.to_numeric(summary["consensus_state"], errors="coerce")
    summary = summary.sort_values(["_state_sort", "consensus_state"]).drop(columns="_state_sort").reset_index(drop=True)
    summary.to_csv(PROJECT_ROOT / config["outputs"]["consensus_state_summary"], index=False)

    for column in ["cell_type", "candidate_reason", "source_author_style", "source_formal_signature"]:
        counts = merged.groupby(["consensus_state", column], observed=True).size().reset_index(name="cells")
        counts.to_csv(stage7_dir / f"state_by_{column}.csv", index=False)
    write_json(stage7_dir / "consensus_upstream_merge_manifest.json", {
        "n_cells": int(len(merged)),
        "n_states": int(merged["consensus_state"].nunique()),
        "required_columns": REQUIRED_COLUMNS,
        "dataset_audit": audit_rows,
        "merge_key": "scvi_obs_index -> source_cell_id -> original candidate_pool_labels.cell_id",
        "known_marker_combo_source": "dataset-specific upstream AnnData metadata when available",
    })
    print(f"Merged {len(merged)} cells across {merged['consensus_state'].nunique()} states")


if __name__ == "__main__":
    main()
