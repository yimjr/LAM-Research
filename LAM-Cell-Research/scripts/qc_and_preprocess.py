"""Run sample-aware QC, normalization, dimensionality reduction and clustering."""

from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
import yaml
import matplotlib.pyplot as plt


def key_for_resolution(resolution: float) -> str:
    return f"leiden_r{str(resolution).replace('.', '_')}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis.yaml")
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / args.config).read_text())
    seed = int(config["random_seed"])
    np.random.seed(seed)
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    accession = config["accession"]
    input_path = ROOT / config["paths"]["interim"] / f"{accession}_combined_raw.h5ad"
    output_dir = ROOT / config["paths"]["processed"]
    result_dir = ROOT / config["paths"]["results"]
    figure_dir = result_dir / "figures"
    table_dir = result_dir / "tables"
    for directory in (output_dir, figure_dir, table_dir):
        directory.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Missing input AnnData: {input_path}; run prepare_matrix.py first")

    adata = ad.read_h5ad(input_path)
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names = adata.obs_names.astype(str)
    adata.layers["counts"] = adata.X.copy()

    symbols = adata.var.get("gene_symbol", pd.Series(adata.var_names, index=adata.var_names)).astype(str)
    upper_symbols = symbols.str.upper()
    adata.var["gene_symbol_upper"] = upper_symbols.values
    adata.var["mt"] = upper_symbols.str.startswith("MT-").values
    adata.var["ribo"] = upper_symbols.str.startswith(("RPS", "RPL")).values
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True, log1p=False)

    qc_config = config["qc"]
    min_genes = int(qc_config["min_genes"])
    min_counts = int(qc_config["min_counts"])
    mt_limits = qc_config["mt_pct_by_assay"]
    assay_limit = adata.obs["assay"].map(mt_limits).fillna(float(mt_limits["unknown"]))
    adata.obs["mt_pct_limit"] = assay_limit.astype(float)
    adata.obs["qc_fail_low_genes"] = adata.obs["n_genes_by_counts"] < min_genes
    adata.obs["qc_fail_low_counts"] = adata.obs["total_counts"] < min_counts
    adata.obs["qc_fail_mt"] = adata.obs["pct_counts_mt"] > adata.obs["mt_pct_limit"]
    adata.obs["qc_pass"] = ~adata.obs[["qc_fail_low_genes", "qc_fail_low_counts", "qc_fail_mt"]].any(axis=1)
    adata.obs["doublet_score"] = np.nan
    adata.obs["doublet_predicted"] = False
    doublet_status = "disabled"

    if qc_config["doublet"].get("enabled", True):
        qc_subset = adata[adata.obs["qc_pass"]].copy()
        # The GEO supplement is a full-barcode matrix. Remove features that
        # are empty in the QC-passed subset before Scrublet to avoid singular
        # sparse decompositions caused by the union of sample feature sets.
        nonzero_features = np.asarray(qc_subset.X.sum(axis=0)).ravel() > 0
        qc_subset = qc_subset[:, nonzero_features].copy()
        try:
            sc.pp.scrublet(
                qc_subset,
                batch_key="sample_id",
                expected_doublet_rate=float(qc_config["doublet"].get("expected_doublet_rate", 0.06)),
                random_state=seed,
            )
            adata.obs.loc[qc_subset.obs_names, "doublet_score"] = qc_subset.obs["doublet_score"].astype(float)
            adata.obs.loc[qc_subset.obs_names, "doublet_predicted"] = qc_subset.obs[
                "predicted_doublet"
            ].astype(bool)
            doublet_status = "completed"
        except Exception as exc:  # diagnostic is retained; pipeline does not silently claim completion
            doublet_status = f"failed: {type(exc).__name__}: {exc}"

    adata.obs["analysis_pass"] = adata.obs["qc_pass"] & ~adata.obs["doublet_predicted"].astype(bool)
    qc_table = adata.obs.copy()
    qc_table.insert(0, "cell_id", qc_table.index)
    qc_table.to_csv(table_dir / "cell_qc_metrics.csv", index=False)
    summary = (
        adata.obs.groupby(["sample_id", "donor_id", "assay", "condition"], observed=True)
        .agg(
            cells_before=("analysis_pass", "size"),
            cells_qc_pass=("qc_pass", "sum"),
            cells_analysis_pass=("analysis_pass", "sum"),
            predicted_doublets=("doublet_predicted", "sum"),
        )
        .reset_index()
    )
    qc_medians = (
        adata.obs.loc[adata.obs["qc_pass"]]
        .groupby(["sample_id", "donor_id", "assay", "condition"], observed=True)
        .agg(
            median_genes=("n_genes_by_counts", "median"),
            median_counts=("total_counts", "median"),
            median_mt_pct=("pct_counts_mt", "median"),
        )
        .reset_index()
    )
    summary = summary.merge(
        qc_medians,
        on=["sample_id", "donor_id", "assay", "condition"],
        how="left",
        validate="one_to_one",
    )
    summary.to_csv(table_dir / "qc_summary_by_sample.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    adata.obs.boxplot(column="n_genes_by_counts", by="sample_id", ax=axes[0], rot=45)
    axes[0].set_title("Detected genes")
    adata.obs.boxplot(column="total_counts", by="sample_id", ax=axes[1], rot=45)
    axes[1].set_title("Total counts")
    adata.obs.boxplot(column="pct_counts_mt", by="sample_id", ax=axes[2], rot=45)
    axes[2].set_title("Mitochondrial percentage")
    fig.suptitle("")
    fig.tight_layout()
    fig.savefig(figure_dir / "qc_overview.png", bbox_inches="tight")
    plt.close(fig)

    filtered = adata[adata.obs["analysis_pass"]].copy()
    if filtered.n_obs < 100:
        raise RuntimeError(f"Too few cells remain after QC: {filtered.n_obs}")
    filtered.uns["qc"] = {
        "min_genes": min_genes,
        "min_counts": min_counts,
        "mt_pct_by_assay": mt_limits,
        "doublet_status": doublet_status,
        "n_obs_before": int(adata.n_obs),
        "n_obs_after": int(filtered.n_obs),
    }
    # h5ad cannot serialize arbitrary nested Python lists/dicts in ``uns``;
    # retain the exact configuration as JSON text for reproducibility.
    filtered.uns["analysis_config_json"] = json.dumps(config, sort_keys=True)
    filtered.uns["source_input"] = str(input_path.relative_to(ROOT))

    # Keep raw counts in layers["counts"], and use float32 for the analysis
    # matrix to reduce memory pressure on the external-disk workstation.
    filtered.X = filtered.X.astype(np.float32)
    sc.pp.normalize_total(filtered, target_sum=float(config["preprocess"]["target_sum"]))
    sc.pp.log1p(filtered)
    if not np.isfinite(filtered.X.data).all():
        raise RuntimeError("Non-finite values detected after normalize_total/log1p")
    filtered.raw = filtered.copy()
    sc.pp.highly_variable_genes(
        filtered,
        n_top_genes=int(config["preprocess"]["n_top_genes"]),
        flavor="cell_ranger",
        batch_key="sample_id",
    )
    # Keep the full log-normalized matrix and all-gene counts in the public
    # AnnData contract. Use a separate HVG view for scaled embeddings.
    work = filtered[:, filtered.var["highly_variable"]].copy()
    # Keep the matrix sparse. A centered 30k x 3k matrix would be densified
    # and can create avoidable numerical and memory failures on this dataset.
    sc.pp.scale(work, zero_center=False, max_value=10)
    n_pcs = min(int(config["preprocess"]["n_pcs"]), max(2, work.n_vars - 1))
    sc.tl.pca(work, n_comps=n_pcs, random_state=seed, svd_solver="arpack", zero_center=False)
    if not np.isfinite(work.obsm["X_pca"]).all():
        raise RuntimeError("Non-finite values detected in PCA representation")

    representation = "X_pca"
    harmony_status = "not_requested"
    if config["preprocess"].get("use_harmony", True):
        try:
            # harmonypy >=0.2 exposes Z_corr as cells x PCs, whereas older
            # Scanpy wrappers expect PCs x cells. Normalize the orientation
            # explicitly for compatibility across both APIs.
            harmony_out = harmonypy.run_harmony(
                work.obsm["X_pca"].astype(np.float64),
                work.obs,
                "sample_id",
                random_state=seed,
            )
            harmony = np.asarray(harmony_out.Z_corr)
            if harmony.shape != work.obsm["X_pca"].shape:
                harmony = harmony.T
            if harmony.shape != work.obsm["X_pca"].shape or not np.isfinite(harmony).all():
                raise RuntimeError(f"Unexpected Harmony output shape or values: {harmony.shape}")
            work.obsm["X_harmony"] = harmony.astype(np.float32)
            representation = "X_harmony"
            harmony_status = "completed"
        except Exception as exc:  # retain fallback state explicitly
            harmony_status = f"failed: {type(exc).__name__}: {exc}"

    sc.pp.neighbors(
        work,
        n_neighbors=int(config["preprocess"]["n_neighbors"]),
        use_rep=representation,
        random_state=seed,
    )
    sc.tl.umap(work, random_state=seed)
    for resolution in config["preprocess"]["leiden_resolutions"]:
        sc.tl.leiden(
            work,
            resolution=float(resolution),
            random_state=seed,
            key_added=key_for_resolution(float(resolution)),
        )

    for key in work.obs.columns:
        if key.startswith("leiden_"):
            filtered.obs[key] = work.obs[key].astype(str).values
    for key in work.obsm.keys():
        filtered.obsm[key] = work.obsm[key]
    for key in work.obsp.keys():
        filtered.obsp[key] = work.obsp[key]
    for key in ("neighbors", "umap", "pca"):
        if key in work.uns:
            filtered.uns[key] = work.uns[key]
    if "X_harmony" in work.obsm:
        filtered.obsm["X_harmony"] = work.obsm["X_harmony"]

    filtered.uns["integration"] = {"representation": representation, "harmony_status": harmony_status}
    output = output_dir / f"{accession}_preprocessed.h5ad"
    filtered.write_h5ad(output, compression="gzip")
    run_manifest = ROOT / config["paths"]["manifests"] / "run_manifest.yaml"
    run_manifest.write_text(
        yaml.safe_dump(
            {
                "status": "preprocessed",
                "project": config["project"],
                "python": sys.version,
                "steps": [
                    {
                        "name": "qc_and_preprocess",
                        "completed_at": datetime.now(timezone.utc).isoformat(),
                        "input": str(input_path.relative_to(ROOT)),
                        "output": str(output.relative_to(ROOT)),
                        "doublet_status": doublet_status,
                        "harmony_status": harmony_status,
                    }
                ],
            },
            sort_keys=False,
        )
    )
    print(json.dumps({"output": str(output), "cells": filtered.n_obs, "genes": filtered.n_vars}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"qc_and_preprocess.py failed: {exc}", file=sys.stderr)
        raise
