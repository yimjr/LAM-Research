"""Estimate broad immune context states and patient-level candidate associations."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import PROJECT_ROOT, ensure_output_path, write_json


def assign_context(table: pd.DataFrame) -> pd.Series:
    marker_modules = {
        "T_cell_candidate": "module_t_cell_markers__n_detected",
        "NK_candidate": "module_nk_markers__n_detected",
        "macrophage_candidate": "module_macrophage_markers__n_detected",
    }
    counts = pd.DataFrame({name: pd.to_numeric(table.get(column, 0), errors="coerce").fillna(0) for name, column in marker_modules.items()}, index=table.index)
    best = counts.idxmax(axis=1)
    best_value = counts.max(axis=1)
    return pd.Series(np.where(best_value >= 2, best, "unassigned"), index=table.index)


def summarize_dataset(table: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = table.copy()
    table["context_label"] = assign_context(table)
    counts = table.groupby(["dataset", "patient_id", "context_label"], observed=True).size().rename("n_cells").reset_index()
    totals = table.groupby(["dataset", "patient_id"], observed=True).size().rename("n_total_cells").reset_index()
    counts = counts.merge(totals, on=["dataset", "patient_id"], how="left")
    counts["fraction_of_cells"] = counts["n_cells"] / counts["n_total_cells"]

    rows = []
    for (dataset, patient), group in table.groupby(["dataset", "patient_id"], observed=True):
        lam = group[group["pool_high_confidence"].astype(bool)]
        row = {"dataset": dataset, "patient_id": patient, "n_cells": len(group), "n_high_confidence_lamcore": len(lam)}
        for module in ["antigen_associated", "presentation_machinery", "ifn_response", "immune_evasion", "nk_ligands"]:
            column = f"module_{module}"
            row[f"lamcore_mean_{module}"] = float(pd.to_numeric(lam[column], errors="coerce").mean()) if len(lam) else np.nan
        for context in ["T_cell_candidate", "NK_candidate", "macrophage_candidate"]:
            immune = group[group["context_label"].eq(context)]
            for module in ["cd8_cytotoxicity", "t_cell_exhaustion", "nk_state", "macrophage_suppressive"]:
                column = f"module_{module}"
                row[f"{context}_mean_{module}"] = float(pd.to_numeric(immune[column], errors="coerce").mean()) if len(immune) else np.nan
        rows.append(row)
    patient = pd.DataFrame(rows)
    return counts, patient


def associations(patient: pd.DataFrame) -> pd.DataFrame:
    rows = []
    target_pairs = [
        ("lamcore_mean_immune_evasion", "T_cell_candidate_mean_t_cell_exhaustion", "candidate evasion-to-T-cell exhaustion axis"),
        ("lamcore_mean_nk_ligands", "NK_candidate_mean_nk_state", "candidate NK-ligand-to-NK-state axis"),
        ("lamcore_mean_presentation_machinery", "T_cell_candidate_mean_cd8_cytotoxicity", "candidate presentation-to-CD8 cytotoxicity axis"),
        ("lamcore_mean_antigen_associated", "macrophage_candidate_mean_macrophage_suppressive", "candidate antigen-to-macrophage context axis"),
    ]
    for dataset, group in patient.groupby("dataset", observed=True):
        for left, right, label in target_pairs:
            pair = group[[left, right]].dropna()
            if len(pair) >= 3 and pair[left].nunique() > 1 and pair[right].nunique() > 1:
                rho, p_value = spearmanr(pair[left], pair[right])
            else:
                rho, p_value = np.nan, np.nan
            rows.append({
                "dataset": dataset,
                "axis": label,
                "left_variable": left,
                "right_variable": right,
                "n_patients": len(pair),
                "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                "p_value_descriptive": float(p_value) if np.isfinite(p_value) else np.nan,
                "interpretation": "candidate association only; not direct communication or causality",
            })
    return pd.DataFrame(rows)


def main() -> None:
    tables = []
    for path in sorted((PROJECT_ROOT / "results" / "cell_scores").glob("*_cell_visibility_scores.csv")):
        tables.append(pd.read_csv(path))
    if not tables:
        raise FileNotFoundError("Run score_visibility_modules.py first.")
    summaries = []
    patient_tables = []
    for table in tables:
        counts, patient = summarize_dataset(table)
        summaries.append(counts)
        patient_tables.append(patient)
    context = pd.concat(summaries, ignore_index=True)
    patient = pd.concat(patient_tables, ignore_index=True)
    assoc = associations(patient)
    output_dir = PROJECT_ROOT / "results" / "immune_context"
    context.to_csv(ensure_output_path(output_dir / "patient_immune_context_counts.csv"), index=False)
    patient.to_csv(ensure_output_path(output_dir / "patient_lamcore_immune_summary.csv"), index=False)
    assoc.to_csv(ensure_output_path(output_dir / "patient_lamcore_immune_associations.csv"), index=False)
    write_json(PROJECT_ROOT / "manifests" / "immune_context_manifest.json", {
        "cell_type_rule": "highest marker detected-count with at least 2 detected marker genes",
        "unassigned_rule": "fewer than 2 detected markers or unresolved tie",
        "interpretation": "All associations are candidate associations, not direct communication claims.",
    })
    print(f"context rows: {len(context)}; patient rows: {len(patient)}; associations: {len(assoc)}")


if __name__ == "__main__":
    main()

