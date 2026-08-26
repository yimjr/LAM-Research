"""Build the Phase 1 baseline without default doublet removal.

The input is a processed GEO count matrix. This script therefore records
downstream QC that remains recoverable from the matrix, while explicitly not
claiming to reproduce FASTQ/Cell Ranger/barcode-calling QC.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import yaml


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/reproduction_core.yaml")
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / args.config).read_text())
    seed = int(config["random_seed"])
    np.random.seed(seed)
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    input_path = ROOT / config["input"]
    output_path = ROOT / config["output"]
    result_dir = ROOT / "results" / "reproduction_core"
    table_dir = result_dir / "tables"
    figure_dir = result_dir / "figures"
    for directory in (output_path.parent, table_dir, figure_dir):
        directory.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(input_path)

    adata = ad.read_h5ad(input_path)
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    adata.layers["counts"] = adata.X.copy()
    if np.issubdtype(adata.layers["counts"].dtype, np.floating):
        values = adata.layers["counts"].data if hasattr(adata.layers["counts"], "data") else np.asarray(adata.layers["counts"])
        if np.any(np.asarray(values) < 0) or not np.allclose(np.asarray(values), np.round(np.asarray(values))):
            raise ValueError("Input counts are not non-negative integers")

    symbols = adata.var.get("gene_symbol", pd.Series(adata.var_names, index=adata.var_names)).astype(str)
    upper = symbols.str.upper()
    adata.var["gene_symbol_upper"] = upper.values
    adata.var["mt"] = upper.str.startswith("MT-").values
    adata.var["ribo"] = upper.str.startswith(("RPS", "RPL")).values
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=False)

    qc = config["qc"]
    mt_limits = qc["mt_pct_by_assay"]
    adata.obs["mt_pct_limit"] = adata.obs["assay"].map(mt_limits).fillna(float(mt_limits["unknown"])).astype(float)
    adata.obs["qc_fail_low_genes"] = adata.obs["n_genes_by_counts"] < int(qc["min_genes"])
    adata.obs["qc_fail_low_counts"] = adata.obs["total_counts"] < int(qc["min_counts"])
    adata.obs["qc_fail_mt"] = adata.obs["pct_counts_mt"] > adata.obs["mt_pct_limit"]
    adata.obs["qc_pass"] = ~adata.obs[["qc_fail_low_genes", "qc_fail_low_counts", "qc_fail_mt"]].any(axis=1)
    adata.obs["doublet_score"] = np.nan
    adata.obs["doublet_predicted"] = False
    doublet_status = "disabled"
    if qc["doublet"].get("enabled", True):
        subset = adata[adata.obs["qc_pass"]].copy()
        nonzero = np.asarray(subset.X.sum(axis=0)).ravel() > 0
        subset = subset[:, nonzero].copy()
        try:
            sc.pp.scrublet(
                subset,
                batch_key="sample_id",
                expected_doublet_rate=float(qc["doublet"].get("expected_doublet_rate", 0.06)),
                random_state=seed,
            )
            adata.obs.loc[subset.obs_names, "doublet_score"] = subset.obs["doublet_score"].astype(float)
            adata.obs.loc[subset.obs_names, "doublet_predicted"] = subset.obs["predicted_doublet"].astype(bool)
            doublet_status = "completed_record_only"
        except Exception as exc:
            doublet_status = f"failed: {type(exc).__name__}: {exc}"

    # This is the explicit Phase 1 rule: doublets are recorded, not removed.
    adata.obs["phase1_included"] = adata.obs["qc_pass"].astype(bool)
    adata.obs["analysis_pass"] = adata.obs["phase1_included"]
    adata.obs["doublet_exclusion_applied"] = False
    adata.obs.insert(0, "cell_id", adata.obs.index)
    adata.obs.to_csv(table_dir / "cell_qc_metrics_phase1.csv", index=False)

    summary = (
        adata.obs.groupby(["sample_id", "donor_id", "tissue", "assay", "condition"], observed=True)
        .agg(
            cells_before=("phase1_included", "size"),
            cells_qc_pass=("qc_pass", "sum"),
            cells_phase1_included=("phase1_included", "sum"),
            predicted_doublets=("doublet_predicted", "sum"),
        )
        .reset_index()
    )
    summary.to_csv(table_dir / "qc_summary_by_sample_phase1.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    adata.obs.boxplot(column="n_genes_by_counts", by="sample_id", ax=axes[0], rot=45)
    axes[0].set_title("Detected genes")
    adata.obs.boxplot(column="total_counts", by="sample_id", ax=axes[1], rot=45)
    axes[1].set_title("Total counts")
    adata.obs.boxplot(column="pct_counts_mt", by="sample_id", ax=axes[2], rot=45)
    axes[2].set_title("Mitochondrial percentage")
    fig.suptitle("")
    fig.tight_layout()
    fig.savefig(figure_dir / "qc_overview_phase1.png", bbox_inches="tight")
    plt.close(fig)

    filtered = adata[adata.obs["phase1_included"]].copy()
    if filtered.n_obs < 100:
        raise RuntimeError(f"Too few cells remain after recoverable QC: {filtered.n_obs}")
    filtered.X = filtered.X.astype(np.float32)
    sc.pp.normalize_total(filtered, target_sum=float(config["preprocess"]["target_sum"]))
    sc.pp.log1p(filtered)
    if not np.isfinite(filtered.X.data).all():
        raise RuntimeError("Non-finite values after normalize_total/log1p")
    filtered.raw = filtered.copy()
    sc.pp.highly_variable_genes(
        filtered,
        n_top_genes=int(config["preprocess"]["n_top_genes"]),
        flavor="cell_ranger",
        batch_key="sample_id",
    )
    work = filtered[:, filtered.var["highly_variable"]].copy()
    sc.pp.scale(work, zero_center=False, max_value=10)
    n_pcs = min(int(config["preprocess"]["n_pcs"]), max(2, work.n_vars - 1))
    sc.tl.pca(work, n_comps=n_pcs, random_state=seed, svd_solver="arpack", zero_center=False)
    sc.pp.neighbors(work, n_neighbors=int(config["preprocess"]["n_neighbors"]), use_rep="X_pca", random_state=seed)
    sc.tl.umap(work, random_state=seed)
    sc.tl.leiden(work, resolution=0.8, random_state=seed, key_added="leiden_r0_8")
    filtered.obs["leiden_r0_8"] = work.obs["leiden_r0_8"].astype(str).values
    for key in work.obsm.keys():
        filtered.obsm[key] = work.obsm[key]
    for key in work.obsp.keys():
        filtered.obsp[key] = work.obsp[key]
    for key in ("neighbors", "umap", "pca"):
        if key in work.uns:
            filtered.uns[key] = work.uns[key]

    filtered.uns["reproduction_baseline"] = {
        "phase": "phase1_core_reproduction",
        "qc_boundary": "downstream QC recoverable from processed matrices; upstream FASTQ/Cell Ranger/barcode-calling QC not reconstructed",
        "doublet_status": doublet_status,
        "doublet_exclusion_applied": False,
        "input": str(input_path.relative_to(ROOT)),
        "config": config,
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    filtered.uns["analysis_config_json"] = json.dumps(config, sort_keys=True)
    filtered.write_h5ad(output_path, compression="gzip")

    manifest_path = ROOT / "manifests" / "run_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("steps", [])
    manifest["steps"] = [step for step in manifest["steps"] if step.get("name") != "phase1_reproduction_baseline"]
    manifest["steps"].append(
        {
            "name": "phase1_reproduction_baseline",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "input": str(input_path.relative_to(ROOT)),
            "output": str(output_path.relative_to(ROOT)),
            "doublet_status": doublet_status,
            "doublet_exclusion_applied": False,
        }
    )
    manifest["status"] = "phase1_baseline_ready_parallel_phase2_phase3"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(json.dumps({"output": str(output_path), "cells": filtered.n_obs, "genes": filtered.n_vars, "doublet_status": doublet_status}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"build_reproduction_baseline.py failed: {exc}", file=sys.stderr)
        raise
