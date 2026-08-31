#!/usr/bin/env python3
"""Apply only the new project's QC, harmonize genes, and build prepared AnnData."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from state_modeling_utils import (
    PROJECT_ROOT,
    apply_registry_mapping,
    calculate_count_qc,
    ensure_counts_layer,
    load_config,
    canonicalize_gene_aliases,
    recreate_log_normalized_x,
    resolve_path,
    resolve_shared,
    threshold_for_assay,
    write_json,
)


def external_qc(adata: ad.AnnData, config: dict) -> pd.DataFrame:
    calculate_count_qc(adata)
    qc = config["qc"]
    assays = adata.obs.get("assay", pd.Series("unknown", index=adata.obs.index))
    mt_thresholds = np.asarray([threshold_for_assay(value, qc) for value in assays], dtype=float)
    adata.obs["qc_fail_low_genes_state_model"] = adata.obs["n_genes_by_counts"].to_numpy() < int(qc["min_genes"])
    adata.obs["qc_fail_low_counts_state_model"] = adata.obs["total_counts"].to_numpy() < int(qc["min_counts"])
    adata.obs["qc_fail_mt_state_model"] = adata.obs["pct_counts_mt"].to_numpy() > mt_thresholds
    adata.obs["state_model_qc_pass"] = ~adata.obs[[
        "qc_fail_low_genes_state_model",
        "qc_fail_low_counts_state_model",
        "qc_fail_mt_state_model",
    ]].any(axis=1).to_numpy()
    adata.obs["qc_pass"] = adata.obs["state_model_qc_pass"].to_numpy()
    adata.obs["state_model_qc_source"] = "state_modeling_external_qc"
    return pd.DataFrame({
        "dataset": adata.obs["dataset"].astype(str).to_numpy(),
        "assay": assays.astype(str).to_numpy(),
        "n_cells": 1,
        "qc_pass": adata.obs["state_model_qc_pass"].astype(bool).to_numpy(),
    }).groupby(["dataset", "assay"], as_index=False).agg(cells=("n_cells", "sum"), qc_pass=("qc_pass", "sum"))


def inherited_path(dataset: str) -> Path:
    return PROJECT_ROOT / "data/interim/inherited" / f"{dataset}.h5ad"


def load_optional_normal(config: dict) -> tuple[ad.AnnData | None, dict]:
    normal_cfg = config.get("normal_reference", {})
    path = resolve_path(normal_cfg.get("h5ad_candidates", []), config.get("input_roots", []))
    if path is None:
        return None, {"status": "not_available"}
    try:
        obj = ad.read_h5ad(path)
        obj.obs["source_cell_id"] = obj.obs_names.astype(str).to_numpy()
        obj.obs["cell_id"] = obj.obs["source_cell_id"].to_numpy()
        obj.obs["dataset"] = "GSE122960"
        obj.obs["source_accession"] = "GSE122960"
        obj, _ = canonicalize_gene_aliases(obj, {"FIGF": "VEGFD"})
        ensure_counts_layer(obj)
        registry_path = resolve_shared(config, config["annotation_files"]["shared"]["donor_registry"])
        if registry_path is not None:
            import yaml
            registry = yaml.safe_load(registry_path.read_text(encoding="utf-8")) or {}
            mapping = apply_registry_mapping(obj, "GSE122960", registry)
            if mapping["n_unresolved"] or mapping["n_conflicts"]:
                return None, {"status": "not_available", "path": str(path), "reason": "normal mapping unresolved or conflicting", "mapping": mapping}
        obj.obs["upstream_pool_high_confidence"] = False
        obj.obs["upstream_pool_broad_lam_like"] = False
        obj.obs["upstream_pool_unrestricted_lam"] = False
        obj.obs["lam_candidate"] = False
        obj.obs["boundary"] = False
        obj.obs["analysis_role"] = "normal_reference"
        obj.obs["condition"] = "control"
        if "cell_type" not in obj.obs:
            obj.obs["cell_type"] = "unknown"
        if "doublet_score" not in obj.obs:
            obj.obs["doublet_score"] = float("nan")
        if "doublet_predicted" not in obj.obs:
            obj.obs["doublet_predicted"] = False
        obj.obs["state_model_qc_source"] = "state_modeling_optional_normal_qc"
        calculate_count_qc(obj)
        qc = config["qc"]
        assays = obj.obs.get("assay", pd.Series("scRNA", index=obj.obs.index))
        thresholds = np.asarray([threshold_for_assay(x, qc) for x in assays], dtype=float)
        obj.obs["state_model_qc_pass"] = (
            (obj.obs["n_genes_by_counts"].to_numpy() >= int(qc["min_genes"]))
            & (obj.obs["total_counts"].to_numpy() >= int(qc["min_counts"]))
            & (obj.obs["pct_counts_mt"].to_numpy() <= thresholds)
        )
        obj.obs["qc_pass"] = obj.obs["state_model_qc_pass"].to_numpy()
        obj = obj[obj.obs["state_model_qc_pass"].to_numpy()].copy()
        recreate_log_normalized_x(obj, float(config["preprocess"]["target_sum"]))
        return obj, {"status": "available", "path": str(path), "n_cells": int(obj.n_obs)}
    except Exception as exc:
        return None, {"status": "not_available", "path": str(path), "reason": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    objects: list[ad.AnnData] = []
    qc_rows: list[pd.DataFrame] = []
    errors: list[str] = []

    for dataset, spec in config["datasets"].items():
        path = inherited_path(dataset)
        if not path.exists():
            errors.append(f"{dataset}: run 02_inherit_and_prepare.py first or input is blocked")
            continue
        try:
            obj = ad.read_h5ad(path)
            if spec.get("kind") == "core":
                if "qc_pass" not in obj.obs:
                    raise ValueError("inherited core QC field qc_pass is missing")
                obj.obs["state_model_qc_pass"] = obj.obs["qc_pass"].map(lambda value: str(value).strip().lower() in {"1", "true", "yes", "y"}).to_numpy()
                obj.obs["state_model_qc_source"] = "inherited_LAM-Cell-Research_qc_pass"
                # Do not recalculate or change core QC.  Only fill metrics when
                # absent, for downstream reporting.
                if not {"n_genes_by_counts", "total_counts", "pct_counts_mt"}.issubset(obj.obs.columns):
                    calculate_count_qc(obj)
            else:
                qc_rows.append(external_qc(obj, config))
            obj = obj[obj.obs["state_model_qc_pass"].astype(bool).to_numpy()].copy()
            recreate_log_normalized_x(obj, float(config["preprocess"]["target_sum"]))
            objects.append(obj)
        except Exception as exc:
            errors.append(f"{dataset}: {type(exc).__name__}: {exc}")

    normal, normal_audit = load_optional_normal(config)
    if normal is not None:
        objects.append(normal)

    if errors:
        print("BLOCKED_INPUT", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 2
    high_count = sum(int(obj.obs.get("lam_candidate", pd.Series(False, index=obj.obs.index)).astype(bool).sum()) for obj in objects)
    if high_count == 0:
        print("BLOCKED_INPUT: no pool_high_confidence cells after inherited QC", file=sys.stderr)
        return 2

    common_genes = list(objects[0].var_names.astype(str))
    for obj in objects[1:]:
        present = set(obj.var_names.astype(str))
        common_genes = [gene for gene in common_genes if gene in present]
    if not common_genes:
        print("BLOCKED_INPUT: no common canonical genes across datasets", file=sys.stderr)
        return 2

    prepared: list[ad.AnnData] = []
    for obj in objects:
        obj = obj[:, common_genes].copy()
        dataset = str(obj.obs["dataset"].iloc[0])
        obj.obs["analysis_cell_id"] = [f"{dataset}::{cell}" for cell in obj.obs["source_cell_id"].astype(str)]
        obj.obs_names = obj.obs["analysis_cell_id"].astype(str).to_numpy()
        prepared.append(obj)

    # Variables have already been reduced to the canonical intersection.
    # Use an outer join for obs so dataset-specific upstream annotations are
    # retained instead of silently dropping them during concatenation.
    combined = ad.concat(prepared, join="outer", merge="first", uns_merge="first", index_unique=None)
    combined.var_names = pd.Index(common_genes, dtype=str)
    combined.obs["included_in_prepared"] = True
    combined.uns["state_model_preparation"] = {
        "datasets": [str(x) for x in combined.obs["dataset"].astype(str).unique()],
        "n_common_genes": len(common_genes),
        "common_gene_policy": "canonicalized_intersection",
        "gse135851_qc": "inherited_without_recalculation",
        "external_qc": "recomputed_for_GSE190260_GSE217108_GSE302356",
        "normal_reference": normal_audit,
        "selection": config["selection"],
    }
    output = PROJECT_ROOT / config["outputs"]["prepared_h5ad"]
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.write_h5ad(output, compression="gzip")
    qc_output = PROJECT_ROOT / "results/stage1_6/qc_summary.csv"
    qc_output.parent.mkdir(parents=True, exist_ok=True)
    if qc_rows:
        pd.concat(qc_rows, ignore_index=True).to_csv(qc_output, index=False)
    else:
        pd.DataFrame(columns=["dataset", "assay", "cells", "qc_pass"]).to_csv(qc_output, index=False)
    write_json(PROJECT_ROOT / "results/stage1_6/preparation_summary.json", {
        "output": str(output),
        "shape": list(combined.shape),
        "high_confidence_cells": int(combined.obs["lam_candidate"].astype(bool).sum()),
        "boundary_cells": int(combined.obs["boundary"].astype(bool).sum()),
        "unrestricted_audit_cells": int(combined.obs["upstream_pool_unrestricted_lam"].astype(bool).sum()),
        "normal_reference": normal_audit,
    })
    print(f"prepared={output} shape={combined.shape} high={high_count} genes={len(common_genes)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
