"""Run post hoc threshold, QC, assay, and donor robustness checks."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def discover_from_threshold(obs: pd.DataFrame, cluster_col: str, score_col: str, detected_col: str, quantile: float, min_detected: int) -> tuple[pd.Series, str]:
    summary = obs.groupby(cluster_col, observed=True).agg(
        mean_score=(score_col, "mean"),
        mean_detected=(detected_col, "mean"),
    )
    cutoff = float(summary["mean_score"].quantile(quantile))
    eligible = summary.index[
        (summary["mean_score"] >= cutoff) & (summary["mean_detected"] >= min_detected)
    ].astype(str)
    candidate = obs[cluster_col].astype(str).isin(set(eligible))
    if int(candidate.sum()) < 10:
        cell_cutoff = float(obs[score_col].quantile(quantile))
        candidate = (obs[score_col] >= cell_cutoff) & (obs[detected_col] >= min_detected)
        return candidate.astype(bool), "cell_score_fallback"
    return candidate.astype(bool), "cluster_mean_score"


def donor_program_summary(obs: pd.DataFrame, candidate: pd.Series, program_columns: list[str], min_cells: int) -> tuple[int, list[str]]:
    lam = obs[obs["condition"].astype(str) == "LAM"].copy()
    lam["candidate_tmp"] = candidate.loc[lam.index].to_numpy()
    tested = 0
    passed: list[str] = []
    for program in program_columns:
        higher = 0
        for _, donor in lam.groupby("donor_id", observed=True):
            cand = donor[donor["candidate_tmp"]]
            other = donor[~donor["candidate_tmp"]]
            if len(cand) < min_cells or len(other) < min_cells:
                continue
            if float(cand[program].mean()) > float(other[program].mean()):
                higher += 1
        valid = sum(
            len(donor[donor["candidate_tmp"]]) >= min_cells and len(donor[~donor["candidate_tmp"]]) >= min_cells
            for _, donor in lam.groupby("donor_id", observed=True)
        )
        tested = max(tested, valid)
        if valid == lam["donor_id"].nunique() and higher >= 3:
            passed.append(program)
    return tested, passed


def summarize_filter(obs: pd.DataFrame, retained: pd.Series, label: str, program_columns: list[str], min_cells: int) -> dict:
    lam = obs[(obs["condition"].astype(str) == "LAM") & retained]
    rows = []
    for donor_id, donor in lam.groupby("donor_id", observed=True):
        formal = donor[donor["lamcore_candidate_formal"].astype(bool)]
        rows.append(
            {
                "filter": label,
                "donor_id": donor_id,
                "assay": donor["assay"].iloc[0],
                "retained_cells": len(donor),
                "formal_candidates_retained": len(formal),
                "formal_candidate_fraction": len(formal) / len(donor) if len(donor) else np.nan,
                "mean_formal_score": float(donor["lamcore_score_formal"].mean()),
            }
        )
    return rows


def main() -> None:
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text())
    accession = config["accession"]
    input_path = ROOT / config["paths"]["processed"] / f"{accession}_lam_states.h5ad"
    table_dir = ROOT / config["paths"]["results"] / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_path)
    obs = adata.obs.copy()
    cluster_col = "leiden_r0_8"
    program_columns = [column for column in obs.columns if column.startswith("program_")]
    min_cells = int(config["analysis"]["min_cells_per_donor_for_summary"])

    threshold_rows = []
    for quantile in (0.80, 0.90, 0.95):
        for min_detected in (5, 10, 20):
            candidate, method = discover_from_threshold(
                obs, cluster_col, "lamcore_score_formal", "lamcore_genes_detected_formal", quantile, min_detected
            )
            tested, passed = donor_program_summary(obs, candidate, program_columns, min_cells)
            lam = obs[obs["condition"].astype(str) == "LAM"]
            threshold_rows.append(
                {
                    "cluster_quantile": quantile,
                    "min_formal_genes_detected": min_detected,
                    "method": method,
                    "all_cells_candidates": int(candidate.sum()),
                    "lam_cells_candidates": int(candidate.loc[lam.index].sum()),
                    "lam_candidate_fraction": float(candidate.loc[lam.index].mean()),
                    "lam_donors_tested": tested,
                    "programs_passing_3_of_4": ";".join(passed),
                }
            )
    pd.DataFrame(threshold_rows).to_csv(table_dir / "robustness_threshold_sensitivity.csv", index=False)

    base_mt = obs["mt_pct_limit"].astype(float)
    not_doublet = ~obs["doublet_predicted"].astype(bool)
    primary = obs["analysis_pass"].astype(bool)
    loose = (
        (obs["n_genes_by_counts"] >= 100)
        & (obs["total_counts"] >= 300)
        & (obs["pct_counts_mt"] <= base_mt + 10.0)
        & not_doublet
    )
    strict_limit = np.maximum(base_mt - 5.0, 5.0)
    strict = (
        (obs["n_genes_by_counts"] >= 500)
        & (obs["total_counts"] >= 1000)
        & (obs["pct_counts_mt"] <= strict_limit)
        & not_doublet
    )
    qc_rows = []
    for label, retained in (("loose", loose), ("primary_analysis_pass", primary), ("strict", strict)):
        qc_rows.extend(summarize_filter(obs, retained, label, program_columns, min_cells))
    pd.DataFrame(qc_rows).to_csv(table_dir / "robustness_qc_sensitivity.csv", index=False)

    assay_rows = []
    for (condition, donor_id, assay), group in obs.groupby(["condition", "donor_id", "assay"], observed=True):
        formal = group[group["lamcore_candidate_formal"].astype(bool)]
        row = {
            "condition": condition,
            "donor_id": donor_id,
            "assay": assay,
            "cells": len(group),
            "formal_candidates": len(formal),
            "formal_candidate_fraction": len(formal) / len(group) if len(group) else np.nan,
            "formal_score_mean": float(group["lamcore_score_formal"].mean()),
            "formal_score_median": float(group["lamcore_score_formal"].median()),
        }
        for program in program_columns:
            row[f"{program}_mean"] = float(group[program].mean())
            row[f"{program}_candidate_minus_other"] = (
                float(formal[program].mean() - group.loc[~group["lamcore_candidate_formal"].astype(bool), program].mean())
                if len(formal) >= min_cells and len(group) - len(formal) >= min_cells
                else np.nan
            )
        assay_rows.append(row)
    pd.DataFrame(assay_rows).to_csv(table_dir / "robustness_assay_summary.csv", index=False)

    summary = {
        "threshold_rows": len(threshold_rows),
        "primary_formal_candidates": int(obs["lamcore_candidate_formal"].sum()),
        "qc_retained_cells": {"loose": int(loose.sum()), "primary_analysis_pass": int(primary.sum()), "strict": int(strict.sum())},
        "formal_candidate_counts_by_filter": {
            label: int((obs["lamcore_candidate_formal"] & retained).sum())
            for label, retained in (("loose", loose), ("primary_analysis_pass", primary), ("strict", strict))
        },
        "assay_groups": int(len(assay_rows)),
    }
    (table_dir / "robustness_summary.json").write_text(json.dumps(summary, indent=2))
    manifest_path = ROOT / config["paths"]["manifests"] / "run_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest["steps"] = [step for step in manifest.get("steps", []) if step.get("name") != "robustness_tests"]
    manifest["steps"].append(
        {
            "name": "robustness_tests",
            "completed_at": pd.Timestamp.now(tz="UTC").isoformat(),
            "input": str(input_path.relative_to(ROOT)),
            "outputs": [
                "results/tables/robustness_threshold_sensitivity.csv",
                "results/tables/robustness_qc_sensitivity.csv",
                "results/tables/robustness_assay_summary.csv",
                "results/tables/robustness_summary.json",
            ],
        }
    )
    manifest["status"] = "analysis_completed_formal_signature_robustness_checked"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
