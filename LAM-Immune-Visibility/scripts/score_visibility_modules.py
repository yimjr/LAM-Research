"""Score immune-visibility modules without modifying source h5ad files."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd

from common import (
    PROJECT_ROOT,
    add_metadata,
    ensure_output_path,
    expression_intervals,
    gene_state_matrix,
    load_expression,
    load_project_config,
    load_signatures,
    load_source_manifest,
    project_relative,
    resolve_source,
    score_modules,
    quality_status,
    write_json,
)


def score_one_dataset(dataset: str, spec: dict, manifest: dict, config: dict, signatures: dict) -> dict:
    matrix_path = resolve_source(spec["matrix"], manifest["source_root"])
    if not matrix_path.exists():
        raise FileNotFoundError(matrix_path)
    target_genes = sorted({gene.upper() for module in signatures.values() for gene in module.get("genes", [])})
    expr = load_expression(matrix_path, target_genes)
    pool_path = resolve_source(spec["pools"], manifest["source_root"]) if spec.get("pools") else None
    obs = add_metadata(expr, dataset, spec, pool_path if pool_path and pool_path.exists() else None)
    obs["technical_quality_status"] = quality_status(obs, config)
    obs["raw_count_basis"] = "counts_layer" if expr.raw_counts_available else "normalized_X_proxy"

    module_scores, availability = score_modules(expr, signatures, config["modules"]["min_available_genes"])

    state_matrix = np.full((expr.obs.shape[0], len(expr.target_genes)), "not_assayed", dtype=object)
    interval_tables = []
    assay_values = obs["assay_label"].astype(str).to_numpy()
    for assay in sorted(pd.unique(assay_values)):
        rows = np.flatnonzero(assay_values == assay)
        sub = replace(
            expr,
            obs=expr.obs.iloc[rows].copy(),
            raw_counts=expr.raw_counts[rows, :],
            normalized=expr.normalized[rows, :],
        )
        intervals = expression_intervals(
            sub,
            config["expression"]["positive_expression_quantile"],
            config["expression"]["min_positive_observations_for_interval"],
        )
        intervals.insert(0, "assay", assay)
        intervals.insert(0, "dataset", dataset)
        interval_tables.append(intervals)
        state_matrix[rows, :] = gene_state_matrix(sub, intervals)

    state_columns = {
        f"state__{gene}": state_matrix[:, index]
        for index, gene in enumerate(expr.target_genes)
    }
    state_frame = pd.DataFrame(state_columns, index=obs.index)
    antigen_frame = pd.DataFrame(index=obs.index)
    for gene in signatures["antigen_associated"]["genes"]:
        gene = gene.upper()
        if gene in expr.gene_to_col:
            antigen_frame[f"count__{gene}"] = expr.raw_counts[:, expr.gene_to_col[gene]]
            antigen_frame[f"norm__{gene}"] = expr.normalized[:, expr.gene_to_col[gene]]
        else:
            antigen_frame[f"count__{gene}"] = np.nan
            antigen_frame[f"norm__{gene}"] = np.nan
    result = pd.concat([obs, module_scores, state_frame, antigen_frame], axis=1)
    result.index.name = "cell_id"

    output = PROJECT_ROOT / "results" / "cell_scores" / f"{dataset}_cell_visibility_scores.csv"
    intervals_output = PROJECT_ROOT / "manifests" / f"{dataset}_expression_intervals.csv"
    availability_output = PROJECT_ROOT / "manifests" / f"{dataset}_module_availability.csv"
    ensure_output_path(output)
    result.to_csv(output)
    pd.concat(interval_tables, ignore_index=True).to_csv(ensure_output_path(intervals_output), index=False)
    availability.insert(0, "dataset", dataset)
    availability.to_csv(ensure_output_path(availability_output), index=False)

    source_stat = matrix_path.stat()
    run_manifest = {
        "dataset": dataset,
        "source_matrix": project_relative(matrix_path),
        "source_size_bytes": source_stat.st_size,
        "source_mtime_ns": source_stat.st_mtime_ns,
        "shape": [int(expr.obs.shape[0]), int(len(expr.target_genes))],
        "available_gene_count": len(expr.available_genes),
        "target_gene_count": len(expr.target_genes),
        "raw_counts_available": expr.raw_counts_available,
        "state_rule": "raw count > 0 means detected; low/high describe expression intervals among nonzero observations",
        "output": project_relative(output),
        "intervals": project_relative(intervals_output),
    }
    write_json(PROJECT_ROOT / "manifests" / f"{dataset}_score_manifest.json", run_manifest)
    return run_manifest


def main() -> None:
    config = load_project_config()
    manifest = load_source_manifest()
    signatures = load_signatures()
    results = []
    for dataset, spec in manifest["datasets"].items():
        results.append(score_one_dataset(dataset, spec, manifest, config, signatures))
    write_json(PROJECT_ROOT / "manifests" / "score_run_manifest.json", {
        "project": config["project"],
        "datasets": results,
        "read_only_sources": True,
        "state_rule": "detected_low means nonzero count and a low-expression interval after detection",
    })
    print(pd.DataFrame(results)[["dataset", "shape", "available_gene_count", "raw_counts_available"]].to_string(index=False))


if __name__ == "__main__":
    main()
