from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "22_state15_local_branch_analysis.py"
OUTPUT_DIR = PROJECT_DIR / "results/stage22"


def test_stage22_is_geometry_only():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "no_scvi_training" in text
    assert "no_reclustering" in text
    assert "no_candidate_gate_change" in text
    assert "no_state15_redefinition" in text
    assert "scvi.model.SCVI" not in text
    assert "sc.tl.leiden" not in text


def test_stage22_manifest_freezes_expected_scope():
    manifest = json.loads((OUTPUT_DIR / "stage22_manifest.json").read_text(encoding="utf-8"))
    assert manifest["main_cell_count"] == 22261
    assert manifest["candidate_cell_count"] == 5378
    assert manifest["boundary_cell_count"] == 16883
    assert manifest["anchor_cell_count"] == 200
    assert manifest["graph_k"] == 30
    assert manifest["no_scvi_training"] is True
    assert manifest["no_reclustering"] is True
    assert manifest["no_candidate_gate_change"] is True


def test_stage22_local_graph_has_frozen_state15_and_hop_layers():
    graph = pd.read_csv(OUTPUT_DIR / "state15_local_graph_cells.csv")
    assert len(graph) == 22261
    assert graph["analysis_cell_id"].is_unique
    assert int(graph["state"].eq("State_15").sum()) == 200
    assert set(graph["min_graph_hop_to_State15"]) == {-1, 0, 1, 2, 3}
    assert int(graph["min_graph_hop_to_State15"].eq(0).sum()) == 200
    assert graph["number_of_State15_neighbors"].ge(0).all()
    assert graph["fraction_of_State15_neighbors"].between(0, 1).all()


def test_stage22_branch_selection_is_automatic_and_excludes_boundary():
    connectivity = pd.read_csv(OUTPUT_DIR / "state15_state_connectivity.csv")
    branches = pd.read_csv(OUTPUT_DIR / "branch_candidates.csv")
    assert "boundary" in set(connectivity["state"])
    assert set(branches["source_state"]) == {"State_16", "State_12", "State_20", "State_7"}
    assert (branches["1hop_cells"] >= 10).all()
    assert (branches["patient_count"] >= 2).all()
    assert not branches["source_state"].eq("boundary").any()


def test_stage22_state16_geometry_and_patient_replication_are_separate_from_scores():
    position = pd.read_csv(OUTPUT_DIR / "state16_branch_position.csv")
    assert len(position) == 396
    assert set(position["distance_segment"]) == {"near", "mid", "far"}
    gradient = pd.read_csv(OUTPUT_DIR / "state16_branch_gradient.csv")
    segment = gradient[gradient["row_type"].eq("segment")]
    assert set(segment["distance_segment"]) == {"near", "mid", "far"}
    patients = pd.read_csv(OUTPUT_DIR / "state16_patient_branch_consistency.csv")
    assert len(patients) == 7
    assert (patients["n_branch_cells"] >= 10).all()
    assert "LAMCORE_direction" in patients.columns


def test_stage22_branch_outputs_cover_all_selected_branches():
    branches = pd.read_csv(OUTPUT_DIR / "branch_candidates.csv")
    ids = set(branches["branch_id"])
    gradients = pd.read_csv(OUTPUT_DIR / "all_branch_gradients.csv")
    assert set(gradients["branch_id"]) == ids
    patient = pd.read_csv(OUTPUT_DIR / "branch_patient_consistency.csv")
    assert set(patient["branch_id"]).issubset(ids)
    extension = pd.read_csv(OUTPUT_DIR / "boundary_branch_extension.csv")
    assert ids.issubset(set(extension["branch_id"]))
    evidence = pd.read_csv(OUTPUT_DIR / "branch_evidence_summary.csv")
    assert set(evidence["branch_id"]) == ids
    assert set(evidence["evidence_label"]).issubset(
        {"LAM_like_branch_candidate", "LAM_to_lineage_transition_candidate", "ordinary_lineage_adjacency", "ambiguous_branch"}
    )


def test_stage22_boundary_assignment_does_not_create_new_lam_labels():
    boundary = pd.read_csv(OUTPUT_DIR / "boundary_local_branch_assignment.csv")
    assert len(boundary) == 16883
    assert set(boundary["state"].dropna()) == {"boundary"}
    assert "unresolved" in set(boundary["branch_assignment"])
    assert (boundary["branch_assignment"] != "LAM").all()


def test_stage22_branch_matched_null_has_500_replicates_per_branch():
    null = pd.read_csv(OUTPUT_DIR / "branch_matched_null.csv")
    assert len(null) == 4 * 500
    assert null.groupby("branch_id")["replicate"].nunique().to_dict() == {
        "branch_01": 500,
        "branch_02": 500,
        "branch_03": 500,
        "branch_04": 500,
    }
    assert set(null["distance_metric"]) == {"Stage20 distance_to_state15"}


def test_stage22_checkpoint_is_local_branch_diagnostic():
    manifest = json.loads((OUTPUT_DIR / "stage22_manifest.json").read_text(encoding="utf-8"))
    assert manifest["checkpoint"] == "local_branched_lam_manifold_candidate"
    assert manifest["branch_count"] == 4
