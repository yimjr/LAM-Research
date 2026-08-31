#!/usr/bin/env python3
"""Inventory inherited AnnData and upstream annotations without recomputing them."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import pandas as pd

from state_modeling_utils import (
    PROJECT_ROOT,
    annotation_directory,
    discover_annotation_files,
    load_config,
    resolve_dataset_h5ad,
    resolve_shared,
    write_json,
)


def candidate_summary(path: Path | None) -> dict:
    if path is None or not path.exists():
        return {"available": False}
    frame = pd.read_csv(path, usecols=lambda column: column in {
        "cell_id", "pool_high_confidence", "pool_broad_lam_like", "pool_unrestricted_lam"
    })
    result = {"available": True, "path": str(path), "n_rows": int(len(frame))}
    for col in ["pool_high_confidence", "pool_broad_lam_like", "pool_unrestricted_lam"]:
        if col in frame:
            result[col] = int(frame[col].astype(str).str.strip().str.lower().isin({"1", "true", "yes", "y"}).sum())
        else:
            result[col] = None
    return result


def inspect_h5ad(path: Path | None) -> dict:
    if path is None:
        return {"available": False}
    try:
        obj = ad.read_h5ad(path, backed="r")
        return {
            "available": True,
            "path": str(path),
            "n_cells": int(obj.n_obs),
            "n_genes": int(obj.n_vars),
            "layers": sorted(obj.layers.keys()),
            "has_counts_layer": "counts" in obj.layers,
            "obs_columns": [str(x) for x in obj.obs.columns],
            "var_columns": [str(x) for x in obj.var.columns],
        }
    except Exception as exc:  # inventory should report bad inputs rather than hide them
        return {"available": False, "path": str(path), "error": f"{type(exc).__name__}: {exc}"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    parser.add_argument("--strict", action="store_true", help="exit non-zero when core input blockers are found")
    args = parser.parse_args()
    config = load_config(args.config)
    rows: list[dict] = []
    annotation_rows: list[dict] = []
    datasets: dict[str, dict] = {}
    blockers: list[str] = []

    for dataset, spec in config.get("datasets", {}).items():
        h5ad_path = resolve_dataset_h5ad(config, dataset)
        directory = annotation_directory(config, dataset)
        candidate_path = directory / "candidate_pool_labels.csv" if directory else None
        h5ad_info = inspect_h5ad(h5ad_path)
        candidate_info = candidate_summary(candidate_path)
        annotations = discover_annotation_files(config, dataset)
        annotation_rows.extend(annotations)
        dataset_info = {
            "dataset": dataset,
            "kind": spec.get("kind"),
            "h5ad": h5ad_info,
            "candidate": candidate_info,
            "annotation_directory": str(directory) if directory else None,
            "annotation_files": annotations,
        }
        datasets[dataset] = dataset_info
        if not h5ad_info.get("available"):
            blockers.append(f"{dataset}: missing or unreadable AnnData")
        if not candidate_info.get("available"):
            blockers.append(f"{dataset}: missing candidate_pool_labels.csv")
        elif any(candidate_info.get(col) is None for col in ["pool_high_confidence", "pool_broad_lam_like", "pool_unrestricted_lam"]):
            blockers.append(f"{dataset}: candidate pool schema is incomplete")
        rows.append({
            "dataset": dataset,
            "kind": spec.get("kind"),
            "h5ad_path": h5ad_info.get("path", ""),
            "h5ad_status": "ready" if h5ad_info.get("available") else "missing_or_invalid",
            "n_cells": h5ad_info.get("n_cells", ""),
            "n_genes": h5ad_info.get("n_genes", ""),
            "counts_layer": h5ad_info.get("has_counts_layer", False),
            "candidate_path": candidate_info.get("path", ""),
            "candidate_status": "ready" if candidate_info.get("available") else "missing",
            "pool_high_confidence": candidate_info.get("pool_high_confidence", ""),
            "pool_broad_lam_like": candidate_info.get("pool_broad_lam_like", ""),
            "pool_unrestricted_lam": candidate_info.get("pool_unrestricted_lam", ""),
        })

    shared = {}
    for name, relpath in config.get("annotation_files", {}).get("shared", {}).items():
        path = resolve_shared(config, relpath)
        shared[name] = {"path": str(path) if path else None, "available": path is not None}
        if name in {"donor_registry", "external_modalities", "known_programs"} and path is None:
            blockers.append(f"shared input missing: {name} ({relpath})")

    normal_cfg = config.get("normal_reference", {})
    normal_path = None
    if normal_cfg.get("enabled", True):
        from state_modeling_utils import resolve_path as _resolve_path
        normal_path = _resolve_path(normal_cfg.get("h5ad_candidates", []), config.get("input_roots", []))
    normal = inspect_h5ad(normal_path)
    normal["status"] = "available" if normal.get("available") else "not_available"

    # A candidate file with zero high-confidence cells is an input blocker for
    # the core question; broad/unrestricted pools are intentionally not fallbacks.
    total_high = sum(int(info.get("candidate", {}).get("pool_high_confidence") or 0) for info in datasets.values())
    if total_high == 0:
        blockers.append("no pool_high_confidence cells found across datasets")

    status = "BLOCKED_INPUT" if blockers else "READY"
    payload = {
        "status": status,
        "blockers": blockers,
        "selection_contract": {
            "lam_candidate": "pool_high_confidence",
            "boundary": "pool_broad_lam_like AND NOT pool_high_confidence",
            "unrestricted": "audit_only",
        },
        "datasets": datasets,
        "shared_inputs": shared,
        "normal_reference": normal,
    }
    out_json = PROJECT_ROOT / config["outputs"]["inventory_json"]
    out_csv = PROJECT_ROOT / config["outputs"]["inventory_csv"]
    write_json(out_json, payload)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    pd.DataFrame(annotation_rows).to_csv(
        PROJECT_ROOT / config["outputs"]["upstream_manifest"], index=False
    )
    print(pd.DataFrame(rows).to_string(index=False))
    print(f"inventory_status={status}")
    if blockers:
        print("blockers:")
        for blocker in blockers:
            print(f"- {blocker}")
    return 2 if args.strict and blockers else 0


if __name__ == "__main__":
    sys.exit(main())
