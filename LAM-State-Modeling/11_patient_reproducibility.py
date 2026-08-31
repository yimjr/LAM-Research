#!/usr/bin/env python3
"""Step 11: convert cell-level states into patient-level evidence."""

from __future__ import annotations

import argparse
import os

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd

from state_modeling_utils import PROJECT_ROOT, load_config, write_json


def state_markers(de: pd.DataFrame, state: str, available_genes: set[str]) -> list[str]:
    if de.empty:
        return []
    frame = de[de["state_id"].astype(str).eq(str(state))].copy()
    if frame.empty:
        return []
    frame["abs_lfc"] = pd.to_numeric(frame["log2FoldChange"], errors="coerce").abs()
    frame["padj_sort"] = pd.to_numeric(frame.get("padj"), errors="coerce").fillna(1.0)
    frame = frame.sort_values(["padj_sort", "abs_lfc"], ascending=[True, False])
    genes = [gene for gene in frame["gene"].astype(str) if gene in available_genes]
    return list(dict.fromkeys(genes))[:30]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    out = PROJECT_ROOT / config["outputs"]["step11_dir"]
    out.mkdir(parents=True, exist_ok=True)
    consensus = ad.read_h5ad(PROJECT_ROOT / config["outputs"]["consensus_h5ad"])
    obs = consensus.obs.copy()
    obs["consensus_state"] = obs["consensus_state"].astype(str)
    obs["patient_id"] = obs["patient_id"].astype(str)
    obs["dataset"] = obs["dataset"].astype(str)
    states = sorted(obs["consensus_state"].unique(), key=lambda value: (len(value), value))
    patients = sorted(obs["patient_id"].unique())
    de_path = PROJECT_ROOT / config["outputs"]["step10_dir"] / "state_de_results.csv"
    de = pd.read_csv(de_path) if de_path.exists() and de_path.stat().st_size else pd.DataFrame()

    x = consensus.X
    gene_names = consensus.var_names.astype(str).tolist()
    gene_positions = {gene: i for i, gene in enumerate(gene_names)}
    patient_rows = []
    min_cells = int(config["step11"].get("min_patient_cells", 5))
    for state in states:
        markers = state_markers(de, state, set(gene_names))
        marker_pos = [gene_positions[gene] for gene in markers]
        for patient in patients:
            patient_mask = obs["patient_id"].eq(patient).to_numpy()
            state_mask = patient_mask & obs["consensus_state"].eq(state).to_numpy()
            state_cells = int(state_mask.sum())
            patient_cells = int(patient_mask.sum())
            if marker_pos and state_cells:
                values = x[state_mask][:, marker_pos]
                if hasattr(values, "toarray"):
                    values = values.toarray()
                signature = float(np.asarray(values).mean())
            else:
                signature = np.nan
            patient_rows.append({
                "state_id": state,
                "patient_id": patient,
                "dataset_count": int(obs.loc[state_mask, "dataset"].nunique()) if state_cells else 0,
                "cells": state_cells,
                "patient_total_lam_cells": patient_cells,
                "fraction": state_cells / patient_cells if patient_cells else 0.0,
                "signature_score": signature,
                "signature_genes": ";".join(markers),
                "support": bool(state_cells >= min_cells),
            })
    patient_table = pd.DataFrame(patient_rows)
    patient_table.to_csv(out / "patient_state_matrix.csv", index=False)

    step7_path = PROJECT_ROOT / config["outputs"]["step7_dir"] / "state_stability_summary.csv"
    structural = pd.read_csv(step7_path) if step7_path.exists() else pd.DataFrame()
    step8_path = PROJECT_ROOT / config["outputs"]["step8_dir"] / "loo_state_summary.csv"
    loo = pd.read_csv(step8_path) if step8_path.exists() else pd.DataFrame()
    rows = []
    for state in states:
        state_patients = patient_table[patient_table["state_id"].eq(state)]
        supported = state_patients[state_patients["support"]]
        state_de = de[de["state_id"].astype(str).eq(state)] if not de.empty else pd.DataFrame()
        if not state_de.empty:
            concordance = pd.to_numeric(state_de["patient_direction_concordance"], errors="coerce")
            biological_direction = float(concordance.mean()) if concordance.notna().any() else np.nan
            significant = pd.to_numeric(state_de["padj"], errors="coerce").le(float(config["step10"]["fdr_alpha"]))
            biological_de_support = float(significant.mean()) if len(significant) else 0.0
            lfc_abs = pd.to_numeric(state_de["log2FoldChange"], errors="coerce").abs()
            biological_effect_support = float((lfc_abs >= float(config["step10"]["min_abs_log2fc"])).mean()) if len(lfc_abs) else 0.0
        else:
            biological_direction = np.nan
            biological_de_support = 0.0
            biological_effect_support = 0.0
        if not loo.empty:
            state_loo = loo[loo["consensus_state"].astype(str).eq(state)]
            patient_loo = state_loo[state_loo["omitted_type"].eq("patient")]
            dataset_loo = state_loo[state_loo["omitted_type"].eq("dataset")]
            baseline = pd.to_numeric(patient_loo["baseline_loo_jaccard"], errors="coerce")
            recovery = pd.to_numeric(patient_loo["loo_recovery_jaccard"], errors="coerce")
            additional = pd.to_numeric(patient_loo["loo_additional_loss"], errors="coerce")
            structural_baseline = float(baseline.mean()) if baseline.notna().any() else np.nan
            structural_recovery = float(recovery.mean()) if recovery.notna().any() else np.nan
            structural_loss = float(additional.mean()) if additional.notna().any() else np.nan
            dataset_recovery = pd.to_numeric(dataset_loo["loo_recovery_jaccard"], errors="coerce")
            dataset_support = float(dataset_recovery.mean()) if dataset_recovery.notna().any() else np.nan
        else:
            structural_baseline = np.nan
            structural_recovery = np.nan
            structural_loss = np.nan
            dataset_support = np.nan
        direct_stability = np.nan
        if not structural.empty and structural["consensus_state"].astype(str).eq(state).any():
            direct_stability = float(structural.loc[structural["consensus_state"].astype(str).eq(state), "mean_within_coassignment"].iloc[0])
        components = [value for value in [direct_stability, structural_recovery] if np.isfinite(value)]
        structural_score = float(np.mean(components)) if components else np.nan
        bio_components = [value for value in [biological_direction, biological_de_support, biological_effect_support] if np.isfinite(value)]
        biological_score = float(np.mean(bio_components)) if bio_components else np.nan
        supported_patient_ids = set(supported["patient_id"].astype(str))
        supported_dataset_count = int(obs.loc[
            obs["patient_id"].isin(supported_patient_ids) & obs["consensus_state"].eq(state),
            "dataset",
        ].nunique())
        rows.append({
            "state_id": state,
            "supported_patients": int(len(supported)),
            "patient_coverage": float(len(supported) / len(patients)) if patients else np.nan,
            "dataset_coverage": supported_dataset_count,
            "structural_stability": structural_score,
            "direct_within_coassignment": direct_stability,
            "mean_loo_recovery": structural_recovery,
            "mean_full_reference_baseline_jaccard": structural_baseline,
            "mean_loo_additional_loss": structural_loss,
            "mean_dataset_loo_recovery": dataset_support,
            "biological_reproducibility": biological_score,
            "patient_direction_concordance": biological_direction,
            "de_support_fraction": biological_de_support,
            "effect_size_support_fraction": biological_effect_support,
            "de_genes": int(len(state_de)),
        })
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "state_reproducibility_summary.csv", index=False)
    write_json(out / "step11_manifest.json", {
        "n_states": len(states),
        "n_patients": len(patients),
        "structural_stability_is_continuous": True,
        "biological_reproducibility_is_continuous": True,
        "scvi_training_called": False,
    })
    (PROJECT_ROOT / "reports/stage11_patient_reproducibility.md").write_text(
        "# Stage 11 patient-level reproducibility\n\n"
        "Structural stability and biological reproducibility are reported as continuous evidence dimensions.\n"
        "Patient-level support is based on patient × state counts and state-specific DE direction/effect evidence.\n",
        encoding="utf-8",
    )
    print(f"Step 11 complete: {len(states)} states × {len(patients)} patients")


if __name__ == "__main__":
    main()
