#!/usr/bin/env python3
"""Step 8: leave-one-patient/dataset robustness on the fixed scVI latent."""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd

from state_modeling_utils import (
    PROJECT_ROOT,
    latent_leiden_labels,
    load_config,
    partition_cluster_matches,
    write_json,
)


def grid_match(left: np.ndarray, right: np.ndarray, left_name: str, right_name: str, **context: object) -> pd.DataFrame:
    frame = partition_cluster_matches(left, right)
    frame.insert(0, "left_partition", left_name)
    frame.insert(1, "right_partition", right_name)
    for key, value in reversed(list(context.items())):
        frame.insert(0, key, value)
    return frame


def best_by_consensus(matches: pd.DataFrame, state_ids: list[str], prefix: str) -> pd.DataFrame:
    rows = []
    for state in state_ids:
        subset = matches[matches["cluster_b"].astype(str).eq(str(state))]
        if subset.empty:
            rows.append({f"{prefix}_best_jaccard": np.nan, f"{prefix}_best_intersection": 0, f"{prefix}_matched_cluster": ""})
        else:
            best = subset.sort_values(["jaccard", "intersection"], ascending=False).iloc[0]
            rows.append({
                f"{prefix}_best_jaccard": float(best["jaccard"]),
                f"{prefix}_best_intersection": int(best["intersection"]),
                f"{prefix}_matched_cluster": str(best["cluster_a"]),
            })
    result = pd.DataFrame(rows, index=pd.Index(state_ids, name="consensus_state")).reset_index()
    return result


def patient_support(frame: pd.DataFrame, cluster_column: str) -> pd.DataFrame:
    counts = frame.groupby([cluster_column, "patient_id"], observed=True).size().reset_index(name="cells")
    return counts[counts["cells"] >= 5].groupby(cluster_column, observed=True)["patient_id"].nunique().rename("qualified_patients")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    out = PROJECT_ROOT / config["outputs"]["step8_dir"]
    out.mkdir(parents=True, exist_ok=True)
    consensus_path = PROJECT_ROOT / config["outputs"]["consensus_h5ad"]
    obj = ad.read_h5ad(consensus_path)
    if "X_scVI" not in obj.obsm:
        raise RuntimeError("BLOCKED_INPUT: consensus_state.h5ad has no X_scVI")
    if "consensus_state" not in obj.obs:
        raise RuntimeError("BLOCKED_INPUT: consensus_state.h5ad has no consensus_state")
    latent = np.asarray(obj.obsm["X_scVI"], dtype=np.float32)
    index = obj.obs.index
    full_consensus = obj.obs["consensus_state"].astype(str).to_numpy()
    full_reference_series = latent_leiden_labels(
        latent,
        index,
        int(config["step8"]["reference_n_neighbors"]),
        float(config["step8"]["reference_resolution"]),
        int(config["random_seed"]),
        key_added="leiden_reference",
    )
    full_reference = full_reference_series.to_numpy(dtype=str)
    state_ids = sorted(np.unique(full_consensus), key=lambda value: (len(value), value))
    all_obs = obj.obs.copy()
    all_obs["full_consensus_state"] = full_consensus
    all_obs["full_reference_cluster"] = full_reference

    full_matches = grid_match(
        full_reference,
        full_consensus,
        "full_reference",
        "full_consensus",
        omitted_type="none",
        omitted_id="",
        metric_scope="all_full_high_confidence_cells",
        retained_cells=int(len(obj)),
    )
    full_matches.to_csv(out / "full_reference_consensus_matches.csv", index=False)

    match_frames = []
    summary_frames = []
    assignment_frames = []
    loo_rows = []
    # Patient and dataset LOO are intentionally sequential to keep memory
    # bounded.  No scVI setup or training happens in this script.
    for omitted_type, field in [("patient", "patient_id"), ("dataset", "dataset")]:
        values = sorted(all_obs[field].astype(str).dropna().unique())
        for omitted in values:
            retained = ~all_obs[field].astype(str).eq(str(omitted)).to_numpy()
            if retained.sum() < 2:
                continue
            loo_labels = latent_leiden_labels(
                latent[retained],
                index[retained],
                int(config["step8"]["reference_n_neighbors"]),
                float(config["step8"]["reference_resolution"]),
                int(config["random_seed"]),
                key_added="leiden_loo",
            ).to_numpy(dtype=str)
            restricted_consensus = full_consensus[retained]
            restricted_reference = full_reference[retained]
            context = {
                "omitted_type": omitted_type,
                "omitted_id": str(omitted),
                "metric_scope": "retained_cells_only",
                "retained_cells": int(retained.sum()),
            }
            baseline = grid_match(restricted_reference, restricted_consensus, "restricted_full_reference", "restricted_full_consensus", **context)
            recovery = grid_match(loo_labels, restricted_consensus, "loo_reference", "restricted_full_consensus", **context)
            loo_reference = grid_match(loo_labels, restricted_reference, "loo_reference", "restricted_full_reference", **context)
            match_frames.extend([baseline, recovery, loo_reference])

            retained_frame = all_obs.loc[retained].copy()
            retained_frame["loo_cluster"] = loo_labels
            baseline_best = best_by_consensus(baseline, state_ids, "baseline")
            recovery_best = best_by_consensus(recovery, state_ids, "loo")
            ref_best = best_by_consensus(loo_reference, state_ids, "loo_to_reference")
            summary = baseline_best.merge(recovery_best, on="consensus_state").merge(ref_best, on="consensus_state")
            summary["omitted_type"] = omitted_type
            summary["omitted_id"] = str(omitted)
            summary["retained_cells"] = int(retained.sum())
            summary["baseline_loo_jaccard"] = summary["baseline_best_jaccard"]
            summary["loo_recovery_jaccard"] = summary["loo_best_jaccard"]
            summary["loo_additional_loss"] = summary["baseline_loo_jaccard"] - summary["loo_recovery_jaccard"]
            summary["patient_support"] = summary["loo_matched_cluster"].map(
                patient_support(retained_frame, "loo_cluster").to_dict()
            ).fillna(0).astype(int)
            summary["retained_consensus_patients"] = summary["consensus_state"].map(
                retained_frame.groupby("full_consensus_state", observed=True)["patient_id"].nunique().to_dict()
            ).fillna(0).astype(int)
            summary["scvi_training_called"] = False
            summary_frames.append(summary)

            assignment = retained_frame[["analysis_cell_id", "patient_id", "dataset", "assay"]].copy() if "analysis_cell_id" in retained_frame else retained_frame[["patient_id", "dataset", "assay"]].copy()
            assignment["scvi_obs_index"] = retained_frame.index.astype(str)
            assignment["full_consensus_state"] = restricted_consensus
            assignment["full_reference_cluster"] = restricted_reference
            assignment["loo_cluster"] = loo_labels
            assignment["omitted_type"] = omitted_type
            assignment["omitted_id"] = str(omitted)
            assignment_frames.append(assignment.reset_index(drop=True))
            loo_rows.append({
                "omitted_type": omitted_type,
                "omitted_id": str(omitted),
                "retained_cells": int(retained.sum()),
                "loo_clusters": int(pd.Series(loo_labels).nunique()),
            })
            del loo_labels, retained_frame, baseline, recovery, loo_reference
            gc.collect()

    if match_frames:
        pd.concat(match_frames, ignore_index=True).to_csv(out / "loo_cluster_matches.csv", index=False)
    else:
        pd.DataFrame().to_csv(out / "loo_cluster_matches.csv", index=False)
    if summary_frames:
        pd.concat(summary_frames, ignore_index=True).to_csv(out / "loo_state_summary.csv", index=False)
    else:
        pd.DataFrame().to_csv(out / "loo_state_summary.csv", index=False)
    if assignment_frames:
        pd.concat(assignment_frames, ignore_index=True).to_csv(out / "loo_assignments.csv", index=False)
    else:
        pd.DataFrame().to_csv(out / "loo_assignments.csv", index=False)
    pd.DataFrame(loo_rows).to_csv(out / "loo_runs.csv", index=False)
    write_json(out / "step8_manifest.json", {
        "reference_n_neighbors": int(config["step8"]["reference_n_neighbors"]),
        "reference_resolution": float(config["step8"]["reference_resolution"]),
        "patient_loo_count": int(sum(row["omitted_type"] == "patient" for row in loo_rows)),
        "dataset_loo_count": int(sum(row["omitted_type"] == "dataset" for row in loo_rows)),
        "overlap_scope": "retained_cells_only",
        "split_merge_matching": True,
        "scvi_training_called": False,
    })
    print(f"Step 8 complete: {len(loo_rows)} LOO runs; all overlaps use retained cells")


if __name__ == "__main__":
    main()
