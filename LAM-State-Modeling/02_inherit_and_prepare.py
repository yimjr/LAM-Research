#!/usr/bin/env python3
"""Inherit upstream AnnData/annotations and establish canonical metadata."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import anndata as ad
import pandas as pd
import yaml

from state_modeling_utils import (
    PROJECT_ROOT,
    apply_registry_mapping,
    attach_candidate_annotation,
    attach_upstream_cell_tables,
    canonicalize_gene_aliases,
    discover_annotation_files,
    ensure_counts_layer,
    load_config,
    resolve_dataset_h5ad,
    resolve_shared,
    recreate_log_normalized_x,
    safe_column,
    write_json,
)


CANONICAL_COLUMNS = {
    "cell_id", "source_cell_id", "dataset", "sample_id", "source_sample",
    "specimen_id", "patient_id", "donor_id", "independence_group", "condition",
    "tissue", "assay", "cell_type", "qc_pass", "doublet_score",
    "doublet_predicted", "mapping_status", "mapping_source", "lam_candidate",
    "boundary", "analysis_role", "source_accession", "state_model_qc_pass",
}


def read_yaml(path: Path) -> dict:
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def attach_original_obs(adata: ad.AnnData, original: pd.DataFrame, dataset: str) -> None:
    for col in original.columns:
        if col in CANONICAL_COLUMNS:
            continue
        target = f"upstream_{dataset}_adata_{safe_column(col)}"
        if target not in adata.obs:
            adata.obs[target] = original[col].to_numpy()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    parser.add_argument("--datasets", nargs="*", default=None)
    args = parser.parse_args()
    config = load_config(args.config)
    registry_path = resolve_shared(config, config["annotation_files"]["shared"]["donor_registry"])
    if registry_path is None:
        print("BLOCKED_INPUT: donor_registry.yaml is missing", file=sys.stderr)
        return 2
    registry = read_yaml(registry_path)
    selected = args.datasets or list(config["datasets"])
    inherited_dir = PROJECT_ROOT / "data/interim/inherited"
    inherited_dir.mkdir(parents=True, exist_ok=True)
    audits: list[dict] = []
    errors: list[str] = []

    for dataset in selected:
        spec = config["datasets"][dataset]
        source = resolve_dataset_h5ad(config, dataset)
        annotation_dir = None
        if source is None:
            errors.append(f"{dataset}: AnnData is missing")
            continue
        try:
            adata = ad.read_h5ad(source)
            if adata.obs_names.duplicated().any():
                raise ValueError("duplicate AnnData observation names cannot be joined safely")
            original_obs = adata.obs.copy()
            adata.obs["source_cell_id"] = adata.obs_names.astype(str).to_numpy()
            adata.obs["cell_id"] = adata.obs["source_cell_id"].to_numpy()
            adata.obs["dataset"] = dataset
            adata.obs["source_accession"] = dataset
            if "condition" not in adata.obs:
                adata.obs["condition"] = "unknown"
            if "assay" not in adata.obs:
                adata.obs["assay"] = "unknown"

            # Keep only canonical metadata in the unprefixed namespace.  All
            # inherited cluster/program/state/QC columns are reattached below
            # as upstream_<dataset>_adata_* so they cannot be mistaken for new
            # model inputs or labels.
            keep = [
                col for col in adata.obs.columns
                if col in CANONICAL_COLUMNS
                and not (col == "qc_pass" and spec.get("kind") != "core")
            ]
            adata.obs = adata.obs[keep].copy()
            if "cell_type" not in adata.obs:
                adata.obs["cell_type"] = "unknown"
            if "doublet_score" not in adata.obs:
                adata.obs["doublet_score"] = float("nan")
            if "doublet_predicted" not in adata.obs:
                adata.obs["doublet_predicted"] = False

            adata, alias_audit = canonicalize_gene_aliases(adata, {"FIGF": "VEGFD"})
            counts_audit, copied_counts = ensure_counts_layer(adata)
            if alias_audit["requires_recomputed_log_x"]:
                recreate_log_normalized_x(adata)
            attach_original_obs(adata, original_obs, dataset)

            mapping_audit = apply_registry_mapping(adata, dataset, registry)
            if spec.get("kind") == "core" and "qc_pass" not in adata.obs:
                raise ValueError("core AnnData must inherit qc_pass from LAM-Cell-Research")
            if "qc_pass" in adata.obs and spec.get("kind") == "core":
                from state_modeling_utils import as_bool
                adata.obs["state_model_qc_pass"] = adata.obs["qc_pass"].map(as_bool).to_numpy()
            else:
                adata.obs["state_model_qc_pass"] = True

            from state_modeling_utils import annotation_directory
            annotation_dir = annotation_directory(config, dataset)
            if annotation_dir is None:
                raise ValueError(f"{dataset}: candidate annotation directory is missing")
            candidate_path = annotation_dir / "candidate_pool_labels.csv"
            if not candidate_path.exists():
                raise ValueError(f"{dataset}: candidate_pool_labels.csv is missing")
            candidate_audit = attach_candidate_annotation(adata, candidate_path, dataset)
            upstream_cell_audits = attach_upstream_cell_tables(adata, annotation_dir, dataset)
            upstream_file_manifest = discover_annotation_files(config, dataset)

            unresolved = int((adata.obs["mapping_status"].astype(str) != "registry").sum())
            conflicts = mapping_audit["n_conflicts"]
            if unresolved or conflicts:
                raise ValueError(f"{dataset}: mapping unresolved={unresolved}, conflicts={conflicts}")

            # Store the audit without putting an AnnData object or other
            # non-HDF5 value into uns.
            adata.uns["state_model_inheritance"] = {
                "dataset": dataset,
                "source_h5ad": str(source),
                "source_kind": spec.get("kind"),
                "alias_audit": alias_audit,
                "counts_audit": counts_audit,
                "counts_copied_from_x": copied_counts,
                "mapping_audit": {key: value for key, value in mapping_audit.items() if key != "conflicts"},
                "candidate_audit": candidate_audit,
                # HDF5 cannot serialize a heterogeneous list of mappings in
                # ``uns``; keep the full audit losslessly as JSON instead.
                "upstream_cell_tables_json": json.dumps(upstream_cell_audits, ensure_ascii=False),
                "upstream_file_manifest_json": json.dumps(upstream_file_manifest, ensure_ascii=False),
            }
            output = inherited_dir / f"{dataset}.h5ad"
            adata.write_h5ad(output, compression="gzip")
            audits.append({
                "dataset": dataset,
                "source_h5ad": str(source),
                "output_h5ad": str(output),
                "n_cells": int(adata.n_obs),
                "n_genes": int(adata.n_vars),
                "n_high_confidence": int(adata.obs["lam_candidate"].sum()),
                "n_boundary": int(adata.obs["boundary"].sum()),
                "n_unrestricted": int(adata.obs["upstream_pool_unrestricted_lam"].sum()),
                "mapping_status": "ready",
            })
            print(f"{dataset}: {adata.shape} high={audits[-1]['n_high_confidence']} boundary={audits[-1]['n_boundary']}")
        except Exception as exc:
            errors.append(f"{dataset}: {type(exc).__name__}: {exc}")

    output_manifest = PROJECT_ROOT / "results/stage1_6/upstream_inheritance.json"
    write_json(output_manifest, {"datasets": audits, "errors": errors, "selection_contract": config["selection"]})
    if errors:
        print("BLOCKED_INPUT", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
