from __future__ import annotations

import json
from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "23_visualize_state15_latent_space.py"
OUTPUT_DIR = PROJECT_DIR / "results/stage23_visualization"


def test_stage23_is_visualization_only():
    text = SCRIPT.read_text(encoding="utf-8")
    assert "no_scvi_training" in text
    assert "no_reclustering" in text
    assert "no_candidate_gate_change" in text
    assert "scvi.model.SCVI" not in text
    assert "sc.tl.leiden" not in text
    assert "consensus_state" not in text


def test_stage23_manifest_freezes_scope_and_sources():
    manifest = json.loads((OUTPUT_DIR / "visualization_manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage"] == 23
    assert manifest["main_cells"] == 22261
    assert manifest["candidate_cells"] == 5378
    assert manifest["boundary_cells"] == 16883
    assert manifest["state15_cells"] == 200
    assert manifest["latent"]["latent_key"] == "X_scVI"
    assert manifest["local_graph"]["k"] == 30
    assert manifest["no_scvi_training"] is True
    assert manifest["no_reclustering"] is True
    assert manifest["no_candidate_gate_change"] is True
    assert manifest["no_upstream_artifact_modified"] is True
    assert set(manifest["highlight_states"]) == {"State_15", "State_16", "State_12", "State_20", "State_7"}


def test_stage23_static_and_interactive_outputs_exist():
    expected = {
        "01_global_latent_umap.png", "01_global_latent_umap.pdf",
        "02_state15_local_knn_states.png", "02_state15_local_knn_states.pdf",
        "02_state15_local_knn_lamcore.png", "02_state15_local_knn_lamcore.pdf",
        "03_distance_lamcore_gradient.png", "03_distance_lamcore_gradient.pdf",
        "04_state16_program_heatmap.png", "04_state16_program_heatmap.pdf",
        "05_branch_matched_null.png", "05_branch_matched_null.pdf",
        "3d_global_umap.html", "3d_state15_local_graph.html", "3d_global_pca.html",
    }
    for name in expected:
        path = OUTPUT_DIR / name
        assert path.exists(), name
        assert path.stat().st_size > 1000, name


def test_stage23_html_has_requested_color_modes_and_hover_fields():
    text = (OUTPUT_DIR / "3d_global_umap.html").read_text(encoding="utf-8").replace(r"\u002f", "/")
    for label in ["State", "LAMCORE_independent", "Patient", "Dataset", "Candidate/Boundary", "Branch"]:
        assert label in text
    for field in ["cell_id", "patient", "dataset", "state", "distance_to_State15", "CORE1", "CORE3", "T_NK", "VSMC/pericyte"]:
        assert field in text


def test_stage23_manifest_records_four_selected_branch_directions():
    manifest = json.loads((OUTPUT_DIR / "visualization_manifest.json").read_text(encoding="utf-8"))
    candidates = manifest["latent"]["branch_candidates"]
    assert {row["source_state"] for row in candidates} == {"State_16", "State_12", "State_20", "State_7"}
    assert manifest["local_graph"]["edges_plotted"] <= 100000
    assert manifest["local_graph"]["edges_plotted"] > 0
