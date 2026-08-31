from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "18_validate_state15_anchor.py"
OUTPUT_DIR = PROJECT_DIR / "results/stage18"


def test_stage18_freezes_state15_without_retraining_or_reclustering():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "no_reclustering" in text
    assert "no_scvi_training" in text
    assert "State 15 expected 200 cells" in text
    assert "scvi.model.SCVI" not in text
    assert "sc.tl.leiden" not in text


def test_stage18_formal_reference_and_frozen_object_are_recorded():
    summary = json.loads((OUTPUT_DIR / "state15_anchor_summary.json").read_text(encoding="utf-8"))
    assert summary["target_state"] == "15"
    assert summary["target_cell_count"] == 200
    assert summary["formal_signature_manifest"]["status"] == "available"
    assert summary["formal_signature_manifest"]["n_genes"] == 777
    assert summary["no_reclustering"] is True
    assert summary["no_scvi_training"] is True


def test_stage18_has_requested_cohorts_and_comparators():
    profile = pd.read_csv(OUTPUT_DIR / "state15_marker_profile.csv")
    assert {"State_15", "State_18", "State_20", "State_12", "State_7", "State_5", "boundary", "normal"}.issubset(
        set(profile["cohort"])
    )
    lamcore = pd.read_csv(OUTPUT_DIR / "state15_lamcore_summary.csv")
    assert {"State_15", "boundary", "normal"}.issubset(set(lamcore["cohort"]))
    assert int(lamcore.loc[lamcore["cohort"].eq("State_15"), "n_cells"].iloc[0]) == 200


def test_stage18_author_enrichment_is_overall_dataset_and_patient_stratified():
    enrichment = pd.read_csv(OUTPUT_DIR / "state15_author_enrichment.csv")
    assert {"overall", "dataset", "patient"}.issubset(set(enrichment["stratum_level"]))
    overall = enrichment.query("stratum_level == 'overall' and stratum == 'all'").iloc[0]
    assert float(overall["enrichment_fold"]) > 1
    assert float(overall["fisher_pvalue_greater"]) < 1e-10


def test_stage18_has_pseudobulk_and_latent_neighborhood_outputs():
    required = [
        "state15_patient_pseudobulk.csv",
        "state15_patient_consistency.csv",
        "state15_latent_neighbors.csv",
        "state15_latent_neighbor_edges.csv",
        "state15_latent_distance_by_cell.csv",
        "state15_latent_distance_gradient.csv",
        "state15_anchor_report.md",
    ]
    assert all((OUTPUT_DIR / name).exists() for name in required)
    neighbors = pd.read_csv(OUTPUT_DIR / "state15_latent_neighbors.csv")
    assert "State_15" in set(neighbors["neighbor_cohort"])
