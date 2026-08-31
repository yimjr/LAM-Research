#!/usr/bin/env python3
"""Step 13: assemble a continuous-evidence LAM state atlas."""

from __future__ import annotations

import argparse
import json
import os

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd

from state_modeling_utils import PROJECT_ROOT, load_config, write_json


def read_csv(path):
    return pd.read_csv(path) if path.exists() and path.stat().st_size else pd.DataFrame()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    out = PROJECT_ROOT / config["outputs"]["step13_dir"]
    out.mkdir(parents=True, exist_ok=True)
    stage7 = read_csv(PROJECT_ROOT / config["outputs"]["step7_dir"] / "state_stability_summary.csv")
    stage8 = read_csv(PROJECT_ROOT / config["outputs"]["step8_dir"] / "loo_state_summary.csv")
    stage10 = read_csv(PROJECT_ROOT / config["outputs"]["step10_dir"] / "state_de_results.csv")
    programs = read_csv(PROJECT_ROOT / config["outputs"]["step10_dir"] / "state_program_scores.csv")
    stage11 = read_csv(PROJECT_ROOT / config["outputs"]["step11_dir"] / "state_reproducibility_summary.csv")
    stage12 = read_csv(PROJECT_ROOT / config["outputs"]["step12_dir"] / "state_auxiliary_summary.csv")
    hierarchy = ad.read_h5ad(PROJECT_ROOT / config["outputs"]["hierarchy_h5ad"])
    states = sorted(hierarchy.obs["consensus_state"].astype(str).unique(), key=lambda value: (len(value), value))

    rows = []
    for state in states:
        state_loo = stage8[stage8["consensus_state"].astype(str).eq(state)] if not stage8.empty else pd.DataFrame()
        patient_loo = state_loo[state_loo["omitted_type"].eq("patient")] if not state_loo.empty else pd.DataFrame()
        dataset_loo = state_loo[state_loo["omitted_type"].eq("dataset")] if not state_loo.empty else pd.DataFrame()
        state_programs = programs[programs["state_id"].astype(str).eq(state)] if not programs.empty else pd.DataFrame()
        if not state_programs.empty:
            state_programs = state_programs.sort_values("delta_state_minus_rest", key=lambda series: series.abs(), ascending=False)
            correspondence = ";".join(state_programs["program_name"].astype(str).head(5))
            program_signal = pd.to_numeric(state_programs["delta_state_minus_rest"], errors="coerce").abs().max()
            genes_present = pd.to_numeric(state_programs["genes_present"], errors="coerce").max()
        else:
            correspondence = ""
            program_signal = np.nan
            genes_present = 0
        state_de = stage10[stage10["state_id"].astype(str).eq(state)] if not stage10.empty else pd.DataFrame()
        has_upstream = bool(state_programs.shape[0] and genes_present > 0)
        hierarchy_state = hierarchy.obs[hierarchy.obs["consensus_state"].astype(str).eq(state)]
        parent = str(hierarchy_state["parent_state"].mode().iloc[0]) if not hierarchy_state.empty and "parent_state" in hierarchy_state and not hierarchy_state["parent_state"].mode().empty else ""
        role = str(hierarchy_state["substate_role"].mode().iloc[0]) if not hierarchy_state.empty and "substate_role" in hierarchy_state and not hierarchy_state["substate_role"].mode().empty else ""
        row = {
            "state_id": state,
            "parent_state": parent,
            "substate_role": role,
            "loo_patient_support": float(pd.to_numeric(patient_loo.get("patient_support"), errors="coerce").mean()) if not patient_loo.empty else np.nan,
            "loo_dataset_support": float(pd.to_numeric(dataset_loo.get("loo_recovery_jaccard"), errors="coerce").mean()) if not dataset_loo.empty else np.nan,
            "full_reference_baseline_jaccard": float(pd.to_numeric(state_loo.get("baseline_loo_jaccard"), errors="coerce").mean()) if not state_loo.empty else np.nan,
            "loo_additional_loss": float(pd.to_numeric(state_loo.get("loo_additional_loss"), errors="coerce").mean()) if not state_loo.empty else np.nan,
            "mean_matched_jaccard": float(pd.to_numeric(patient_loo.get("loo_recovery_jaccard"), errors="coerce").mean()) if not patient_loo.empty else np.nan,
            "dataset_coverage": int(hierarchy_state["dataset"].nunique()) if not hierarchy_state.empty else 0,
            "normal_distance": np.nan,
            "boundary_connectivity": np.nan,
            "upstream_program_correspondence": correspondence,
            "novel_or_unexplained": not has_upstream,
            "de_genes": int(len(state_de)),
            "program_signal": float(program_signal) if np.isfinite(program_signal) else np.nan,
        }
        if not stage11.empty and stage11["state_id"].astype(str).eq(state).any():
            source = stage11.loc[stage11["state_id"].astype(str).eq(state)].iloc[0]
            for column in ["structural_stability", "biological_reproducibility", "patient_direction_concordance"]:
                row[column] = source.get(column, np.nan)
            row["patient_coverage"] = source.get("patient_coverage", np.nan)
        else:
            row.update({"structural_stability": np.nan, "biological_reproducibility": np.nan, "patient_direction_concordance": np.nan, "patient_coverage": np.nan})
        if not stage12.empty and stage12["state_id"].astype(str).eq(state).any():
            source = stage12.loc[stage12["state_id"].astype(str).eq(state)].iloc[0]
            row["normal_distance"] = source.get("normal_mean_distance", np.nan)
            row["boundary_connectivity"] = source.get("mean_neighbor_states", np.nan)
        rows.append(row)
    atlas = pd.DataFrame(rows)
    numeric = ["structural_stability", "biological_reproducibility", "mean_matched_jaccard", "patient_direction_concordance"]
    atlas["continuous_evidence_score"] = atlas[numeric].mean(axis=1, skipna=True)
    atlas.to_csv(out / "state_atlas.csv", index=False)
    hypotheses = atlas.sort_values("continuous_evidence_score", ascending=False).head(int(config["step13"].get("hypothesis_top_n", 10))).copy()
    hypotheses["hypothesis_basis"] = "continuous structural/biological evidence; not a confidence threshold"
    hypotheses.to_csv(out / "state_hypothesis_candidates.csv", index=False)
    payload = {
        "scope": "LAM-State-Modeling Step 7-13",
        "n_states": int(len(atlas)),
        "no_hard_confidence_threshold": True,
        "scvi_training_called": False,
        "states": atlas.replace({np.nan: None}).to_dict(orient="records"),
    }
    (out / "state_atlas.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    atlas_h5ad = hierarchy.copy()
    atlas_by_state = atlas.set_index("state_id")
    for column in ["structural_stability", "biological_reproducibility", "loo_patient_support", "loo_dataset_support", "normal_distance", "boundary_connectivity", "continuous_evidence_score"]:
        atlas_h5ad.obs[f"atlas_{column}"] = atlas_h5ad.obs["consensus_state"].astype(str).map(atlas_by_state[column]).to_numpy()
    atlas_h5ad.uns["step13_atlas"] = {
        "n_states": int(len(atlas)),
        "no_hard_confidence_threshold": True,
        "hypothesis_candidate_count": int(len(hypotheses)),
        "scvi_training_called": False,
    }
    atlas_h5ad.write(PROJECT_ROOT / config["outputs"]["atlas_h5ad"])
    write_json(out / "step13_manifest.json", atlas_h5ad.uns["step13_atlas"])
    (PROJECT_ROOT / "reports/state_atlas.md").write_text(
        "# LAM State Atlas\n\n"
        f"- Consensus states retained: {len(atlas)}\n"
        "- Structural stability and biological reproducibility remain continuous evidence dimensions.\n"
        "- No hard high/medium/low confidence threshold was applied.\n"
        "- Hypothesis candidates are an evidence-based reading list, not a new state definition.\n",
        encoding="utf-8",
    )
    print(f"Step 13 complete: atlas assembled for {len(atlas)} states")


if __name__ == "__main__":
    main()
