from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "21_validate_state15_manifold.py"
OUTPUT_DIR = PROJECT_DIR / "results/stage21"


def test_stage21_is_anchor_excluded_and_does_not_train_or_recluster():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "anchor_excluded_from_primary_gradient" in text
    assert "no_scvi_training" in text
    assert "no_reclustering" in text
    assert "no_candidate_gate_change" in text
    assert "scvi.model.SCVI" not in text
    assert "sc.tl.leiden" not in text


def test_stage21_manifest_freezes_anchor_and_validation_object():
    manifest = json.loads((OUTPUT_DIR / "stage21_manifest.json").read_text(encoding="utf-8"))
    assert manifest["anchor_cell_count"] == 200
    assert manifest["validation_cell_count"] == 22061
    assert manifest["candidate_null_pool_count"] == 5178
    assert manifest["latent_key"] == "X_scVI"
    assert manifest["no_scvi_training"] is True
    assert manifest["no_reclustering"] is True
    assert manifest["no_candidate_gate_change"] is True
    assert manifest["distance_contract"]["column"] == "distance_to_state15"


def test_stage21_lamcore_independence_audit_has_expected_components():
    audit = pd.read_csv(OUTPUT_DIR / "lamcore_independence_gene_audit.csv")
    assert len(audit) == 777
    required = {
        "in_scvi_4000_hvg",
        "in_old_candidate_gate_marker",
        "in_neither_scvi_hvg_nor_gate",
        "in_expression_matrix",
    }
    assert required.issubset(audit.columns)
    assert int(audit["in_scvi_4000_hvg"].sum()) == 220
    assert int(audit["in_old_candidate_gate_marker"].sum()) == 7
    assert int(audit["in_neither_scvi_hvg_nor_gate"].sum()) == 554
    assert int(audit["in_expression_matrix"].sum()) == 729


def test_stage21_scores_and_anchor_excluded_gradient_cover_all_validation_cells():
    scores = pd.read_csv(OUTPUT_DIR / "independent_lamcore_scores.csv")
    assert len(scores) == 22061
    assert scores["analysis_cell_id"].is_unique
    assert not scores["current_state"].eq("State_15").any()
    for column in ["LAMCORE_full", "LAMCORE_no_gate", "LAMCORE_outside_scVI", "LAMCORE_independent"]:
        assert column in scores.columns
    gradient = pd.read_csv(OUTPUT_DIR / "non_state15_distance_gradient.csv")
    assert len(gradient) == 6
    assert gradient["n_cells"].sum() == 22061
    smooth = pd.read_csv(OUTPUT_DIR / "distance_score_smooth.csv")
    assert len(smooth) == 20


def test_stage21_models_have_all_scopes_and_dataset_patient_outputs():
    models = pd.read_csv(OUTPUT_DIR / "gradient_models.csv")
    assert set(models["scope"]) == {"all_non_state15", "non_state15_candidates"}
    assert set(models["score_name"]) == {"LAMCORE_full", "LAMCORE_no_gate", "LAMCORE_outside_scVI", "LAMCORE_independent"}
    assert len(pd.read_csv(OUTPUT_DIR / "dataset_independent_gradient.csv")) == 16
    patients = pd.read_csv(OUTPUT_DIR / "patient_independent_gradient.csv")
    assert len(patients) == 48
    assert patients["patient"].nunique() == 12
    assert "n_State15" in patients.columns
    assert "gradient_class" in patients.columns


def test_stage21_matched_null_has_500_replicates_and_empirical_result():
    null = pd.read_csv(OUTPUT_DIR / "matched_anchor_null.csv")
    assert len(null) == 500
    assert null["replicate"].is_unique
    assert (null["sampling_with_replacement_cells"] == 0).all()
    manifest = json.loads((OUTPUT_DIR / "stage21_manifest.json").read_text(encoding="utf-8"))
    assert manifest["matched_anchor_null"]["repetitions"] == 500
    assert manifest["checkpoint_details"]["matched_anchor_empirical_pvalue"] <= 0.05


def test_stage21_state16_boundary_connectivity_and_checkpoint_are_diagnostic_only():
    state16 = pd.read_csv(OUTPUT_DIR / "state16_distance_gradient.csv")
    assert set(state16.loc[state16["scope"].eq("pooled"), "distance_segment"]) == {"near", "mid", "far"}
    assert state16["patient"].nunique() >= 1
    boundary = pd.read_csv(OUTPUT_DIR / "boundary_independent_gradient.csv")
    assert len(boundary) == 6
    connectivity = pd.read_csv(OUTPUT_DIR / "manifold_connectivity.csv")
    assert (connectivity["source_state"] == "State_15").any()
    manifest = json.loads((OUTPUT_DIR / "stage21_manifest.json").read_text(encoding="utf-8"))
    assert manifest["checkpoint"] == "state15_lam_rich_gradient_but_not_robust_manifold"
