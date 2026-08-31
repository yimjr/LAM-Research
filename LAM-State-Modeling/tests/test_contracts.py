from __future__ import annotations

from pathlib import Path
import sys

PROJECT_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_DIR))

import numpy as np
import pandas as pd
import anndata as ad
import yaml

from state_modeling_utils import (
    attach_candidate_annotation,
    canonicalize_gene_aliases,
    ensure_counts_layer,
    model_mask,
    validate_integer_counts,
)

def test_float_dtype_integer_valued_counts_are_valid():
    obj = ad.AnnData(X=np.asarray([[0.0, 1.0], [2.0, 0.0]], dtype=np.float32))
    audit, copied = ensure_counts_layer(obj)
    assert audit["valid"] is True
    assert copied is True
    assert audit["dtype"] == "float32"
    assert np.array_equal(np.asarray(obj.layers["counts"]), obj.X)


def test_non_integer_counts_are_invalid():
    obj = ad.AnnData(X=np.asarray([[0.5, 1.0]], dtype=np.float32))
    audit = validate_integer_counts(obj)
    assert audit["valid"] is False
    assert audit["integer_valued"] is False


def test_figf_is_merged_into_vegfd_counts():
    obj = ad.AnnData(X=np.asarray([[1.0, 2.0, 3.0], [0.0, 4.0, 1.0]], dtype=np.float32))
    obj.layers["counts"] = obj.X.copy()
    obj.var_names = ["FIGF", "VEGFD", "ACTA2"]
    obj, audit = canonicalize_gene_aliases(obj, {"FIGF": "VEGFD"})
    assert list(obj.var_names) == ["VEGFD", "ACTA2"]
    assert np.array_equal(np.asarray(obj.layers["counts"]), [[3.0, 3.0], [4.0, 1.0]])
    assert audit["merged_groups"]["VEGFD"] == ["FIGF", "VEGFD"]


def test_only_high_confidence_is_primary_candidate(tmp_path):
    obj = ad.AnnData(X=np.zeros((3, 1), dtype=np.float32))
    obj.obs_names = ["LAM1:c1", "LAM1:c2", "LAM1:c3"]
    obj.obs["source_cell_id"] = obj.obs_names
    table = pd.DataFrame({
        "cell_id": obj.obs_names,
        "pool_high_confidence": [True, False, False],
        "pool_broad_lam_like": [True, True, False],
        "pool_unrestricted_lam": [True, True, True],
        "candidate_reason": ["high", "broad", "guardrail"],
    })
    path = tmp_path / "candidate_pool_labels.csv"
    table.to_csv(path, index=False)
    attach_candidate_annotation(obj, path, "GSE135851")
    assert obj.obs["lam_candidate"].tolist() == [True, False, False]
    assert obj.obs["boundary"].tolist() == [False, True, False]
    assert obj.obs["analysis_role"].tolist() == ["primary_candidate", "boundary", "unrestricted_audit_only"]
    assert model_mask(obj).tolist() == [True, True, False]


def test_selection_contract_has_no_pool_fallback():
    config = yaml.safe_load((PROJECT_DIR / "config/state_modeling.yaml").read_text())
    assert config["selection"]["lam_candidate_source"] == "pool_high_confidence"
    assert config["selection"]["boundary_definition"] == "pool_broad_lam_like AND NOT pool_high_confidence"
    assert config["selection"]["no_pool_fallback"] is True
    assert config["selection"]["unrestricted_is_audit_only"] is True


def test_scvi_contract_uses_dataset_only_and_keeps_assay_metadata():
    config = yaml.safe_load((PROJECT_DIR / "config/state_modeling.yaml").read_text())
    assert config["scvi"]["batch_key"] == "dataset"
    assert config["scvi"]["categorical_covariate_keys"] == []
    assert "assay" not in config["scvi"]["categorical_covariate_keys"]
    assert config["normal_reference"]["optional"] is True


def test_qc_contract_distinguishes_core_inheritance_and_external_recalculation():
    text = (PROJECT_DIR / "03_qc_and_harmonize.py").read_text()
    assert "inherited_LAM-Cell-Research_qc_pass" in text
    assert "external_qc(obj, config)" in text
    config = yaml.safe_load((PROJECT_DIR / "config/state_modeling.yaml").read_text())
    assert config["qc"]["min_genes"] == 200
    assert config["qc"]["min_counts"] == 500
    assert config["qc"]["mt_pct_by_assay"] == {"scrna": 20.0, "snrna": 10.0, "default": 20.0}


def test_nmf_and_scvi_matrix_contracts_are_separate():
    nmf_text = (PROJECT_DIR / "04_baseline_pca_nmf.py").read_text()
    scvi_text = (PROJECT_DIR / "05_train_scvi.py").read_text()
    assert "recreate_log_normalized_x(obj" in nmf_text
    assert 'setup_anndata(obj, layer="counts", batch_key="dataset")' in scvi_text
    assert 'accelerator=accelerator' in scvi_text
    assert 'batch_key="assay"' not in scvi_text


def test_old_state_correspondence_is_posthoc_not_no_go_gate():
    text = (PROJECT_DIR / "06_stage6_checkpoint.py").read_text()
    assert "posthoc_interpretation" in text
    assert "upstream state/program correspondence is descriptive only" in text
    assert 'status = "NO_GO"' in text
    assert "latent structure is multi-cluster and cross-patient" in text


def test_stage6_reclusters_high_confidence_cells_only():
    text = (PROJECT_DIR / "06_stage6_checkpoint.py").read_text()
    assert 'high_obj = obj[high].copy()' in text
    assert 'use_rep="X_scVI"' in text
    assert 'key_added="neighbors_scvi_lam_only"' in text
    assert 'lam_cluster_col = "leiden_scvi_lam_only"' in text
    assert 'stage6["lam_n_neighbors"]' in text
    assert 'stage6["lam_leiden_resolution"]' in text
    assert 'config["preprocess"]["n_neighbors"]' not in text
    assert 'config["preprocess"]["leiden_resolution"]' not in text
    assert "full_cohort_cluster_count_retained_for_audit" in text
    assert "stage6_cluster_patient_counts.csv" in text
    assert "stage6_parameter_grid.csv" in text
    assert "stage6_grid_pairwise_ari.csv" in text
    assert "boundary_connectivity(obj, high_obj, config)" in text
    assert "normal_neighbors(obj, high_obj, config)" in text


def test_preparation_preserves_dataset_specific_upstream_obs_columns():
    text = (PROJECT_DIR / "03_qc_and_harmonize.py").read_text()
    assert 'ad.concat(prepared, join="outer"' in text


def test_step7_uses_nine_equal_weight_configurations_and_seed_audit():
    text = (PROJECT_DIR / "07_consensus_stability.py").read_text()
    assert "config_coassignment / float(n_configs)" in text
    assert "len(seed_labels)" in text
    assert "configuration_equal_weight" in text
    assert "full_distance_used_for_average_linkage" in text
    manifest = PROJECT_DIR / "results/stage7/step7_manifest.json"
    if manifest.exists():
        payload = yaml.safe_load(manifest.read_text())
        assert payload["n_grid_configurations"] == 9
        assert payload["n_partitions"] == 21


def test_step7_uses_complete_final_distance_for_average_linkage():
    text = (PROJECT_DIR / "07_consensus_stability.py").read_text()
    assert "squareform(1.0 - c_final" in text
    assert "linkage(distance, method=\"average\")" in text
    assert "np.float32" in text


def test_step8_restricts_all_loo_comparisons_to_retained_cells():
    text = (PROJECT_DIR / "08_loo_robustness.py").read_text()
    assert "retained = ~all_obs[field]" in text
    assert "restricted_consensus = full_consensus[retained]" in text
    assert "restricted_reference = full_reference[retained]" in text
    assert '"metric_scope": "retained_cells_only"' in text
    assert "loo_additional_loss" in text
    assert "full_reference_consensus_matches.csv" in text


def test_step10_has_one_independent_patient_aware_model_per_state():
    text = (PROJECT_DIR / "10_biology_annotation.py").read_text()
    assert "for state in states" in text
    assert 'design_factors=["patient_id", "group"]' in text
    assert 'contrast=["group", state_group, "Rest_of_LAM"]' in text
    assert '"~ patient_id + group"' in text
    assert "if n_patients < int(config[\"step10\"][\"min_de_patients\"])" in text
    assert "state_vs_rest" not in text


def test_steps_10_to_13_never_train_scvi():
    for name in ["10_biology_annotation.py", "11_patient_reproducibility.py", "12_boundary_normal_validation.py", "13_state_atlas.py"]:
        text = (PROJECT_DIR / name).read_text()
        assert "scvi_training_called" in text
        assert "train" not in text.lower() or "training_called" in text.lower()


def test_boundary_and_normal_are_auxiliary_only():
    text = (PROJECT_DIR / "12_boundary_normal_validation.py").read_text()
    assert '"states_defined_by": "step7 consensus only"' in text
    assert '"boundary_and_normal_participate_in_state_count": False' in text


def test_consensus_upstream_merge_contract_and_outputs():
    text = (PROJECT_DIR / "14_merge_consensus_upstream.py").read_text()
    for column in ["consensus_state", "cell_type", "candidate_reason", "source_author_style", "source_formal_signature", "known_marker_combo_ge2", "doublet_score", "doublet_predicted", "dataset", "patient_id"]:
        assert f'"{column}"' in text
    output = PROJECT_DIR / "results/stage7/state_consensus_with_upstream_annotations.csv"
    summary = PROJECT_DIR / "results/stage7/state_consensus_state_summary.csv"
    if output.exists() and summary.exists():
        cells = pd.read_csv(output)
        states = pd.read_csv(summary)
        assert len(cells) == 5378
        assert cells["candidate_annotation_matched"].all()
        assert int(states["cells"].sum()) == 5378
