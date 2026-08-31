from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT_DIR / "17_identity_calibration_audit.py"


def load_module():
    spec = importlib.util.spec_from_file_location("stage17_audit", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_stage17_is_read_only_and_does_not_train_or_cluster():
    text = SCRIPT.read_text(encoding="utf-8").lower()
    assert "no_scvi_training" in text
    assert "no_reclustering" in text
    assert "stage 16 gate" in text


def test_stage17_failure_categories_follow_decomposition():
    module = load_module()
    thresholds = {
        "identity_anchor_boundary_threshold": 1.0,
        "support_boundary_threshold": 1.0,
        "LAM_identity_score_boundary_threshold": 1.0,
        "LAM_identity_score_core_threshold": 3.0,
    }
    base = {
        "LAM_core_candidate": False,
        "LAM_boundary_candidate": False,
        "identity_anchor_score": 0.0,
        "support_score": 0.0,
        "LAM_identity_score": 0.0,
        "competing_lineage_exclusion": False,
    }
    detailed, root, identity_failed, competing_failed = module.classify_failure(pd.Series(base), thresholds)
    assert (detailed, root, identity_failed, competing_failed) == (
        "insufficient_positive_evidence",
        "low_identity_evidence_only",
        True,
        False,
    )
    base["competing_lineage_exclusion"] = True
    detailed, root, identity_failed, competing_failed = module.classify_failure(pd.Series(base), thresholds)
    assert (detailed, root, identity_failed, competing_failed) == ("both", "both", True, True)


def test_stage17_outputs_cover_four_datasets_and_formal_reference():
    output_dir = PROJECT_DIR / "results/stage17"
    required = [
        "positive_reference_failures.csv",
        "component_scores_by_dataset.csv",
        "marker_detection_by_dataset.csv",
        "competing_lineage_penalty_audit.csv",
        "counterfactual_calibration.csv",
        "lodo_recovery.csv",
        "root_cause_by_dataset.csv",
        "identity_calibration_audit.md",
    ]
    assert all((output_dir / name).exists() for name in required)
    manifest = json.loads((output_dir / "stage17_manifest.json").read_text(encoding="utf-8"))
    assert manifest["formal_signature_manifest"]["status"] == "available"
    assert manifest["formal_signature_manifest"]["n_genes"] == 777
    assert set(pd.read_csv(output_dir / "root_cause_by_dataset.csv")["dataset"]) == {
        "GSE135851",
        "GSE190260",
        "GSE217108",
        "GSE302356",
    }
    failures = pd.read_csv(output_dir / "positive_reference_failures.csv")
    assert len(failures) == 2257
    assert int((failures["dataset"] == "GSE190260").sum()) == 2117
