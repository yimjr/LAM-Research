from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "19_state15_cross_patient_audit.py"
OUTPUT_DIR = PROJECT_DIR / "results/stage19"


def test_stage19_reuses_state_labels_without_reclustering_or_scvi_training():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "no_reclustering" in text
    assert "no_scvi_training" in text
    assert "scvi.model.SCVI" not in text
    assert "sc.tl.leiden" not in text


def test_stage19_patient_composition_quantifies_lam1163_enrichment():
    composition = pd.read_csv(OUTPUT_DIR / "state15_patient_composition.csv")
    row = composition.loc[composition["patient_id"].eq("LAM1163")].iloc[0]
    assert int(row["candidate_pool_total_cells"]) == 500
    assert int(row["state15_cells"]) == 127
    assert abs(float(row["candidate_pool_fraction"]) - 500 / 5378) < 1e-8
    assert abs(float(row["state15_fraction"]) - 127 / 200) < 1e-8
    assert abs(float(row["enrichment"]) - 6.830060) < 1e-5


def test_stage19_author_availability_does_not_call_unassayed_datasets_negative():
    availability = pd.read_csv(OUTPUT_DIR / "author_annotation_availability.csv")
    status = dict(zip(availability["dataset"], availability["author_style_annotation_status"]))
    assert status["GSE135851"] == "available"
    assert status["GSE190260"] == "not_assayed"
    assert status["GSE217108"] == "not_assayed"
    assert status["GSE302356"] == "not_assayed"
    assayed = pd.read_csv(OUTPUT_DIR / "state15_author_enrichment_assayed.csv")
    assert set(assayed["dataset"]) == {"GSE135851"}


def test_stage19_patient_matched_and_lam1163_sensitivity_support_profile_persistence():
    matched = pd.read_csv(OUTPUT_DIR / "state15_patient_matched_comparison.csv")
    lamcore = matched.loc[matched["feature"].eq("LAMCORE_777")]
    assert len(lamcore) == 7
    assert (lamcore["delta_mean"] > 0).all()

    sensitivity = pd.read_csv(OUTPUT_DIR / "state15_without_LAM1163.csv")
    lamcore_sensitivity = sensitivity.loc[sensitivity["feature"].eq("LAMCORE_777")]
    assert set(lamcore_sensitivity["comparator_group"]) == {"State_18", "State_20", "State_12", "State_7", "State_5"}
    assert (lamcore_sensitivity["target_median"] > lamcore_sensitivity["comparator_median"]).all()

    manifest = json.loads((OUTPUT_DIR / "stage19_manifest.json").read_text(encoding="utf-8"))
    assert manifest["sensitivity"]["state15_cells_after_removal"] == 73
    assert manifest["stage19_conclusion"].startswith("B_")


def test_stage19_lopo_has_seven_held_out_patients_and_five_comparators():
    lopo = pd.read_csv(OUTPUT_DIR / "state15_lopo_validation.csv")
    assert lopo["held_out_patient"].nunique() == 7
    assert lopo["comparison_group"].nunique() == 5
    assert len(lopo) == 35
    assert float(lopo["reference_latent_closer"].mean()) > 0.9
