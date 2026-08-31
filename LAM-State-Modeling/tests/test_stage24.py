from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
OUTPUT_DIR = PROJECT_DIR / "results/stage24_final"
SCRIPT = PROJECT_DIR / "24_finalize_project.py"


def test_stage24_freezes_existing_core_artifacts():
    manifest = json.loads((OUTPUT_DIR / "stage24_manifest.json").read_text(encoding="utf-8"))
    assert manifest["stage"] == 24
    assert manifest["stage_count"] == 23
    assert manifest["stages_indexed"] == list(range(1, 24))
    assert manifest["core_artifacts_frozen"] is True
    assert manifest["no_scvi_training"] is True
    assert manifest["no_reclustering"] is True
    assert manifest["no_candidate_gate_change"] is True
    assert manifest["stage22_prior_state16_transition_label_withdrawn"] is True
    assert manifest["stage22_current_checkpoint"] == "ordinary_lineage_adjacency_dominates"
    assert manifest["expected_core_artifacts_missing"] == []


def test_stage24_indexes_all_states_and_key_artifacts():
    states = pd.read_csv(OUTPUT_DIR / "state_human_cell_analogue.csv")
    assert states["state_id"].tolist() == list(range(1, 21))
    assert states.loc[states["state_id"] == 15, "consensus_cells"].iloc[0] == 200
    assert states.loc[states["state_id"] == 16, "consensus_cells"].iloc[0] == 396

    artifacts = pd.read_csv(OUTPUT_DIR / "artifact_index.csv")
    indexed = set(artifacts["relative_path"].astype(str))
    for path in [
        "24_finalize_project.py",
        "results/stage7/state_consensus_assignments.csv",
        "results/stage13/state_atlas.csv",
        "results/stage18/state15_anchor_summary.json",
        "results/stage20/state15_centered_manifold.csv",
        "results/stage21/gradient_models.csv",
        "results/stage22/branch_evidence_summary.csv",
        "results/stage22/branch_patient_lopo.csv",
        "results/stage23_visualization/visualization_manifest.json",
    ]:
        assert path in indexed


def test_stage24_preserves_known_history_audits():
    audit = pd.read_csv(OUTPUT_DIR / "narrative_audit.csv")
    text = " ".join(audit["issue"].astype(str))
    assert "Stage20 checkpoint" in text
    assert "Stage6 12 clusters" in text
    assert "Stage16 says formal 777 signature unavailable" in text
    assert "Stage21 independent score" in text

    materials = (OUTPUT_DIR / "final_project_source_materials.md").read_text(encoding="utf-8")
    for phrase in ["State15", "State16", "GSE190260", "LAM1163", "ordinary lineage adjacency"]:
        assert phrase in materials


def test_stage24_is_read_only_synthesis_code():
    source = SCRIPT.read_text(encoding="utf-8")
    assert "scvi.model.SCVI" not in source
    assert "sc.tl.leiden" not in source
    assert "no_scvi_training" in source
    assert "no_reclustering" in source
