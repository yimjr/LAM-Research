from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "20_state15_centered_manifold.py"
OUTPUT_DIR = PROJECT_DIR / "results/stage20"


def test_stage20_reuses_frozen_scvi_geometry_without_training_or_reclustering():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "no_scvi_training" in text
    assert "no_reclustering" in text
    assert "no_candidate_gate_change" in text
    assert "scvi.model.SCVI" not in text
    assert "sc.tl.leiden" not in text
    assert "X_scVI" in text


def test_stage20_manifest_freezes_state15_and_main_cohort():
    manifest = json.loads((OUTPUT_DIR / "stage20_manifest.json").read_text(encoding="utf-8"))
    assert manifest["frozen_input"]["state15_cell_count"] == 200
    assert manifest["cohort"]["main_cell_count"] == 22261
    assert manifest["cohort"]["candidate_cell_count"] == 5378
    assert manifest["cohort"]["boundary_cell_count"] == 16883
    assert manifest["frozen_input"]["latent_key"] == "X_scVI"
    assert manifest["no_scvi_training"] is True
    assert manifest["no_reclustering"] is True
    assert manifest["no_candidate_gate_change"] is True


def test_stage20_distance_outputs_have_expected_scope_and_features():
    distances = pd.read_csv(OUTPUT_DIR / "state15_cell_distances.csv")
    manifold = pd.read_csv(OUTPUT_DIR / "state15_centered_manifold.csv")
    assert len(distances) == 22261
    assert len(manifold) == 22261
    assert distances["analysis_cell_id"].is_unique
    assert int(distances["current_state"].eq("State_15").sum()) == 200
    assert set(distances["analysis_role"]) == {"primary_candidate", "boundary"}
    assert set(distances.loc[distances["analysis_role"].eq("boundary"), "current_state"].fillna("")) == {""}
    for column in [
        "nearest_state15_distance",
        "mean_5_nearest_state15_distance",
        "mean_15_nearest_state15_distance",
        "state15_centroid_distance",
        "state15_neighbor_fraction",
        "state15_neighbor_count_k30",
    ]:
        assert column in distances.columns
    assert set(manifold["analysis_role"]) == {"primary_candidate", "boundary"}


def test_stage20_distance_bins_and_gradients_cover_the_main_cohort():
    bins = pd.read_csv(OUTPUT_DIR / "state15_distance_bins.csv")
    expected_bins = {"0-10%", "10-20%", "20-40%", "40-60%", "60-80%", "80-100%"}
    assert set(bins["distance_bin"]) == expected_bins
    assert bins.groupby("distance_bin")["cells"].sum().to_dict() == {
        "0-10%": 2226,
        "10-20%": 2226,
        "20-40%": 4452,
        "40-60%": 4452,
        "60-80%": 4452,
        "80-100%": 4453,
    }
    identity = pd.read_csv(OUTPUT_DIR / "state15_identity_gradient.csv")
    lineage = pd.read_csv(OUTPUT_DIR / "state15_lineage_gradient.csv")
    assert set(identity["distance_bin"]) == expected_bins
    assert set(lineage["distance_bin"]) == expected_bins
    assert "LAMCORE_777_median" in identity.columns
    assert "macrophage_median" in lineage.columns


def test_stage20_state16_and_auxiliary_outputs_are_present():
    coexpression = pd.read_csv(OUTPUT_DIR / "state16_lam_immune_coexpression.csv")
    categories = set(coexpression.loc[coexpression["row_type"].eq("category"), "category"])
    assert categories == {
        "LAM-high / immune-low",
        "LAM-high / immune-high",
        "LAM-low / immune-high",
        "LAM-low / immune-low",
    }
    assert int(coexpression.loc[coexpression["row_type"].eq("category"), "cells"].sum()) == 396
    assert len(pd.read_csv(OUTPUT_DIR / "state16_cell_audit.csv")) == 396
    technical = pd.read_csv(OUTPUT_DIR / "state16_doublet_audit.csv")
    technical_counts = dict(zip(technical["technical_cohort"], technical["n_cells"]))
    assert technical_counts["State_15"] == 200
    assert technical_counts["State_16"] == 396
    assert len(pd.read_csv(OUTPUT_DIR / "boundary_state15_projection.csv")) == 16883
    assert len(pd.read_csv(OUTPUT_DIR / "patient_gradient_consistency.csv")) == 12
    assert len(pd.read_csv(OUTPUT_DIR / "dataset_gradient_consistency.csv")) == 4
    assert len(pd.read_csv(OUTPUT_DIR / "normal_remote_summary.csv")) == 1


def test_stage20_checkpoint_is_geometry_only_and_does_not_add_normal_to_main_table():
    manifest = json.loads((OUTPUT_DIR / "stage20_manifest.json").read_text(encoding="utf-8"))
    assert manifest["distance_geometry"]["reference"] == "frozen State 15 cells only"
    assert manifest["distance_geometry"]["candidate_boundary_graph_k"] == 30
    assert "normal_reference_scope" in manifest["frozen_input"]
    manifold = pd.read_csv(OUTPUT_DIR / "state15_centered_manifold.csv")
    assert not manifold["dataset"].astype(str).str.contains("normal", case=False).any()
