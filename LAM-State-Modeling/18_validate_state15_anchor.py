#!/usr/bin/env python3
"""Validate the frozen consensus State 15 as a possible LAM-core anchor.

This is a read-only validation stage.  It does not change candidate gates,
recluster cells, or retrain scVI.  State 15 is taken exactly from the current
consensus CSV and is compared with the requested consensus states, boundary
cells, and normal/control cells.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import importlib.util
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy.stats import fisher_exact, mannwhitneyu
from sklearn.neighbors import NearestNeighbors


PROJECT_ROOT = Path(__file__).resolve().parent
DATASETS = ["GSE135851", "GSE190260", "GSE217108", "GSE302356"]
TARGET_STATE = "15"
COMPARATORS = ["18", "20", "12", "7", "5"]
ALIASES = {"FIGF": "VEGFD"}
MARKER_GENES = [
    "PMEL", "MLANA", "MITF", "ACTA2", "ACTG2", "MYH11",
    "VEGFD", "CTSK", "EMX2", "HOXA11", "ESR1",
]


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def load_stage16_module() -> Any:
    path = PROJECT_ROOT / "16_rebuild_lam_identity_gate.py"
    spec = importlib.util.spec_from_file_location("stage16_gate", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Stage 16 module: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def as_bool(value: Any) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return False
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def bool_series(values: pd.Series) -> pd.Series:
    return values.map(as_bool).astype(bool)


def canonical_gene(gene: str) -> str:
    upper = str(gene).strip().upper()
    return ALIASES.get(upper, upper)


def unique_genes(genes: list[str]) -> list[str]:
    result: list[str] = []
    for gene in genes:
        canonical = canonical_gene(gene)
        if canonical and canonical not in result:
            result.append(canonical)
    return result


def resolve_formal_signature(config: dict[str, Any]) -> tuple[list[str], dict[str, Any]]:
    roots = [PROJECT_ROOT / str(root) for root in config.get("input_roots", [])]
    roots.extend(
        [
            Path("/mnt/e/LAM-Research/data-temp"),
            Path("/mnt/e/LAM-Research/LAM-Cell-Research"),
            PROJECT_ROOT / "data/upstream",
        ]
    )
    candidates: list[Path] = []
    for root in roots:
        candidates.extend(
            [
                root / "data/raw/reference/LAM_core_signature_genes.csv",
                root / "LAM_core_signature_genes.csv",
            ]
        )
    seen: set[str] = set()
    for path in candidates:
        key = str(path.resolve())
        if key in seen or not path.exists():
            continue
        seen.add(key)
        table = pd.read_csv(path)
        table.columns = [str(column).strip() for column in table.columns]
        gene_column = next(
            (column for column in ["Gene", "gene", "gene_symbol", "symbol"] if column in table.columns),
            str(table.columns[0]),
        )
        genes = unique_genes(table[gene_column].dropna().astype(str).tolist())
        return genes, {
            "status": "available" if genes else "empty",
            "path": str(path),
            "gene_column": gene_column,
            "n_genes": len(genes),
        }
    return [], {"status": "not_available", "path": "", "n_genes": 0}


def load_programs(config: dict[str, Any]) -> tuple[dict[str, list[str]], dict[str, Any]]:
    roots = [PROJECT_ROOT / str(root) for root in config.get("input_roots", [])]
    roots.extend(
        [
            Path("/mnt/e/LAM-Research/LAM-Cell-Research"),
            Path("/mnt/e/LAM-Research/data-temp"),
            PROJECT_ROOT / "data/upstream",
        ]
    )
    seen: set[str] = set()
    for root in roots:
        path = root / "config/known_lam_programs.yaml"
        key = str(path.resolve())
        if key in seen or not path.exists():
            continue
        seen.add(key)
        payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        programs: dict[str, list[str]] = {}
        for entry in payload.get("programs", []):
            name = str(entry.get("program_name", ""))
            if name:
                programs[name] = unique_genes([str(gene) for gene in entry.get("genes", [])])
        return programs, {"status": "available", "path": str(path), "programs": programs}
    return {}, {"status": "not_available", "path": "", "programs": {}}


def prepared_obs(prepared: ad.AnnData) -> pd.DataFrame:
    requested = [
        "cell_id", "sample_id", "donor_id", "tissue", "condition", "assay",
        "source_accession", "source_sample", "dataset", "cell_type", "patient_id",
        "specimen_id", "independence_group", "boundary", "lam_candidate", "analysis_role",
        "total_counts", "n_genes_by_counts",
    ]
    present = [column for column in requested if column in prepared.obs.columns]
    obs = prepared.obs[present].copy().reset_index(names="obs_name")
    if "analysis_cell_id" in prepared.obs.columns:
        obs["analysis_cell_id"] = prepared.obs["analysis_cell_id"].astype(str).to_numpy()
    else:
        obs["analysis_cell_id"] = obs["obs_name"].astype(str)
    for column in ["dataset", "condition", "patient_id", "assay", "cell_type"]:
        if column in obs:
            obs[column] = obs[column].astype(str)
    return obs


def selected_score_table(
    prepared: ad.AnnData,
    selected: pd.DataFrame,
    modules: dict[str, list[str]],
    block_size: int,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    obs_names = pd.Index(prepared.obs["analysis_cell_id"].astype(str)) if "analysis_cell_id" in prepared.obs else pd.Index(prepared.obs_names.astype(str))
    ids = selected["analysis_cell_id"].astype(str).tolist()
    positions = obs_names.get_indexer(ids)
    missing = [ids[index] for index, value in enumerate(positions) if value < 0]
    if missing:
        raise ValueError(f"Selected IDs missing from prepared AnnData: {len(missing)}")
    var_names = pd.Index(prepared.var_names.astype(str))
    actual_by_canonical: dict[str, list[str]] = {}
    for actual in var_names:
        actual_by_canonical.setdefault(canonical_gene(actual), []).append(str(actual))
    canonical_union: list[str] = []
    for genes in modules.values():
        for gene in unique_genes(genes):
            if gene in actual_by_canonical and gene not in canonical_union:
                canonical_union.append(gene)
    actual_union: list[str] = []
    canonical_indices: dict[str, list[int]] = {}
    for gene in canonical_union:
        canonical_indices[gene] = []
        for actual in actual_by_canonical[gene]:
            canonical_indices[gene].append(len(actual_union))
            actual_union.append(actual)

    output = selected[["analysis_cell_id"]].copy().reset_index(drop=True)
    for start in range(0, len(selected), block_size):
        stop = min(start + block_size, len(selected))
        rows = np.asarray(positions[start:stop], dtype=np.int64)
        block = prepared[rows, actual_union].to_memory()
        values = block.X.toarray() if hasattr(block.X, "toarray") else np.asarray(block.X)
        result_block: dict[str, np.ndarray] = {}
        for name, genes in modules.items():
            present = [gene for gene in unique_genes(genes) if gene in canonical_indices]
            if present:
                result_block[name] = np.column_stack(
                    [values[:, canonical_indices[gene]].sum(axis=1) for gene in present]
                ).mean(axis=1).astype(np.float32)
            else:
                result_block[name] = np.full(stop - start, np.nan, dtype=np.float32)
        for name, block_values in result_block.items():
            if name not in output:
                output[name] = np.nan
            output.loc[start:stop - 1, name] = block_values
    return output, {
        "n_selected_cells": len(selected),
        "n_resolved_canonical_genes": len(canonical_union),
        "resolved_canonical_genes": canonical_union,
        "expression_source": "prepared AnnData X (library-size normalized log1p)",
        "block_size": block_size,
    }


def extract_pseudobulk_counts(
    prepared: ad.AnnData,
    selected: pd.DataFrame,
    genes: list[str],
) -> tuple[pd.DataFrame, list[str]]:
    obs_names = pd.Index(prepared.obs["analysis_cell_id"].astype(str)) if "analysis_cell_id" in prepared.obs else pd.Index(prepared.obs_names.astype(str))
    positions = obs_names.get_indexer(selected["analysis_cell_id"].astype(str).tolist())
    if (positions < 0).any():
        raise ValueError("Pseudobulk IDs missing from prepared AnnData")
    var_names = pd.Index(prepared.var_names.astype(str))
    actual_by_canonical: dict[str, list[str]] = {}
    for actual in var_names:
        actual_by_canonical.setdefault(canonical_gene(actual), []).append(str(actual))
    canonical_genes = [gene for gene in unique_genes(genes) if gene in actual_by_canonical]
    actual_union: list[str] = []
    indices: dict[str, list[int]] = {}
    for gene in canonical_genes:
        indices[gene] = []
        for actual in actual_by_canonical[gene]:
            indices[gene].append(len(actual_union))
            actual_union.append(actual)
    subset = prepared[np.asarray(positions, dtype=np.int64), actual_union].to_memory()
    values = subset.layers["counts"]
    values = values.toarray() if hasattr(values, "toarray") else np.asarray(values)
    result = np.zeros((len(selected), len(canonical_genes)), dtype=np.float64)
    for index, gene in enumerate(canonical_genes):
        result[:, index] = values[:, indices[gene]].sum(axis=1)
    return pd.DataFrame(result, columns=canonical_genes), canonical_genes


def add_metadata_to_scores(score_table: pd.DataFrame, selected: pd.DataFrame) -> pd.DataFrame:
    metadata_columns = [column for column in selected.columns if column != "analysis_cell_id"]
    return selected[["analysis_cell_id", *metadata_columns]].merge(
        score_table, on="analysis_cell_id", how="left", validate="one_to_one"
    )


def summary_rows(table: pd.DataFrame, feature_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for cohort, sub in table.groupby("cohort", observed=True):
        for feature in feature_columns:
            values = pd.to_numeric(sub[feature], errors="coerce").dropna().to_numpy(dtype=float)
            rows.append(
                {
                    "cohort": str(cohort),
                    "cohort_type": "consensus_state" if str(cohort).startswith("State_") else str(cohort),
                    "n_cells": len(sub),
                    "n_datasets": sub["dataset"].nunique() if "dataset" in sub else np.nan,
                    "n_patients": sub["patient_id"].nunique() if "patient_id" in sub else np.nan,
                    "feature": feature,
                    "n": len(values),
                    "mean": float(np.mean(values)) if len(values) else np.nan,
                    "median": float(np.median(values)) if len(values) else np.nan,
                    "q05": float(np.quantile(values, 0.05)) if len(values) else np.nan,
                    "q25": float(np.quantile(values, 0.25)) if len(values) else np.nan,
                    "q75": float(np.quantile(values, 0.75)) if len(values) else np.nan,
                    "q95": float(np.quantile(values, 0.95)) if len(values) else np.nan,
                    "detection_fraction_score_gt_0": float(np.mean(values > 0)) if len(values) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def fisher_enrichment(consensus: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    consensus = consensus.copy()
    consensus["state15"] = consensus["consensus_state"].astype(str).eq(TARGET_STATE)
    consensus["author"] = bool_series(consensus["source_author_style"]) if "source_author_style" in consensus else False
    strata = [("overall", "all")]
    strata.extend(("dataset", str(value)) for value in sorted(consensus["dataset"].astype(str).unique()))
    strata.extend(("patient_id", str(value)) for value in sorted(consensus["patient_id"].astype(str).unique()))
    for level, value in strata:
        sub = consensus if value == "all" else consensus[consensus[level].astype(str).eq(value)]
        a = int((sub["state15"] & sub["author"]).sum())
        b = int((sub["state15"] & ~sub["author"]).sum())
        c = int((~sub["state15"] & sub["author"]).sum())
        d = int((~sub["state15"] & ~sub["author"]).sum())
        odds_ratio, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
        state_fraction = a / (a + b) if a + b else np.nan
        other_fraction = c / (c + d) if c + d else np.nan
        fold = state_fraction / other_fraction if other_fraction and np.isfinite(other_fraction) else np.inf if state_fraction > 0 else np.nan
        rows.append(
            {
                "stratum_level": "patient" if level == "patient_id" else level,
                "stratum": value,
                "state15_cells": a + b,
                "state15_author_style": a,
                "other_cells": c + d,
                "other_author_style": c,
                "state15_author_fraction": state_fraction,
                "other_author_fraction": other_fraction,
                "enrichment_fold": fold,
                "fisher_odds_ratio": float(odds_ratio),
                "fisher_pvalue_greater": float(pvalue),
            }
        )
    return pd.DataFrame(rows)


def compare_states(table: pd.DataFrame, feature_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for comparator in COMPARATORS:
        target = table[table["consensus_state"].astype(str).eq(TARGET_STATE)]
        other = table[table["consensus_state"].astype(str).eq(comparator)]
        for feature in feature_columns:
            x = pd.to_numeric(target[feature], errors="coerce").dropna().to_numpy(dtype=float)
            y = pd.to_numeric(other[feature], errors="coerce").dropna().to_numpy(dtype=float)
            if len(x) and len(y):
                statistic, pvalue = mannwhitneyu(x, y, alternative="two-sided")
            else:
                statistic, pvalue = np.nan, np.nan
            rows.append(
                {
                    "comparison": f"State_{TARGET_STATE}_vs_State_{comparator}",
                    "target_state": f"State_{TARGET_STATE}",
                    "comparator_state": f"State_{comparator}",
                    "feature": feature,
                    "n_target": len(x),
                    "n_comparator": len(y),
                    "target_median": float(np.median(x)) if len(x) else np.nan,
                    "comparator_median": float(np.median(y)) if len(y) else np.nan,
                    "median_difference": float(np.median(x) - np.median(y)) if len(x) and len(y) else np.nan,
                    "mannwhitney_u": float(statistic) if np.isfinite(statistic) else np.nan,
                    "mannwhitney_pvalue": float(pvalue) if np.isfinite(pvalue) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def patient_pseudobulk(
    prepared: ad.AnnData,
    consensus: pd.DataFrame,
    modules: dict[str, list[str]],
) -> pd.DataFrame:
    groups = [TARGET_STATE, *COMPARATORS]
    selected = consensus[consensus["consensus_state"].astype(str).isin(groups)].copy()
    selected = selected.reset_index(drop=True)
    selected["group"] = selected["consensus_state"].astype(str).map(lambda value: f"State_{value}")
    genes = unique_genes([gene for values in modules.values() for gene in values])
    counts, present_genes = extract_pseudobulk_counts(prepared, selected, genes)
    counts.insert(0, "patient_id", selected["patient_id"].astype(str).to_numpy())
    counts.insert(1, "group", selected["group"].to_numpy())
    counts.insert(2, "dataset", selected["dataset"].astype(str).to_numpy())
    if "total_counts" in selected:
        total_counts = pd.to_numeric(selected["total_counts"], errors="coerce").to_numpy()
    elif "total_counts" in prepared.obs:
        prepared_ids = pd.Index(prepared.obs["analysis_cell_id"].astype(str)) if "analysis_cell_id" in prepared.obs else pd.Index(prepared.obs_names.astype(str))
        prepared_totals = pd.Series(pd.to_numeric(prepared.obs["total_counts"], errors="coerce").to_numpy(), index=prepared_ids)
        total_counts = selected["analysis_cell_id"].astype(str).map(prepared_totals).to_numpy()
    else:
        total_counts = np.full(len(selected), np.nan)
    counts["total_counts_cell_metadata"] = total_counts
    rows: list[dict[str, Any]] = []
    for (patient, group), index in counts.groupby(["patient_id", "group"], observed=True).groups.items():
        sub = counts.loc[index]
        total_umi = float(sub["total_counts_cell_metadata"].sum())
        if not np.isfinite(total_umi) or total_umi <= 0:
            total_umi = float(sub[present_genes].to_numpy(dtype=float).sum())
        summed = sub[present_genes].sum(axis=0).to_numpy(dtype=float)
        normalized = np.log1p(summed / max(total_umi, 1.0) * 10000.0)
        source = selected.loc[index]
        for module, module_genes in modules.items():
            present = [gene for gene in unique_genes(module_genes) if gene in present_genes]
            rows.append(
                {
                    "patient_id": str(patient),
                    "group": str(group),
                    "cells": int(len(index)),
                    "dataset_count": int(source["dataset"].astype(str).nunique()),
                    "total_umi": total_umi,
                    "module": module,
                    "n_genes_present": len(present),
                    "pseudobulk_log1p_score": float(np.mean([normalized[present_genes.index(gene)] for gene in present])) if present else np.nan,
                }
            )
    return pd.DataFrame(rows)


def latent_outputs(
    config: dict[str, Any],
    consensus: pd.DataFrame,
    score_table: pd.DataFrame,
    output_dir: Path,
) -> dict[str, Any]:
    scvi_path = PROJECT_ROOT / str(config["outputs"]["scvi_h5ad"])
    if not scvi_path.exists():
        return {"status": "not_available", "path": str(scvi_path)}
    latent_obj = ad.read_h5ad(scvi_path, backed="r")
    if "X_scVI" not in latent_obj.obsm:
        latent_obj.file.close()
        return {"status": "missing_X_scVI", "path": str(scvi_path)}
    obs = latent_obj.obs.copy().reset_index(names="obs_name")
    if "analysis_cell_id" in obs:
        obs["analysis_cell_id"] = obs["analysis_cell_id"].astype(str)
    else:
        obs["analysis_cell_id"] = obs["obs_name"].astype(str)
    consensus_map = dict(zip(consensus["analysis_cell_id"].astype(str), consensus["consensus_state"].astype(str)))
    obs["consensus_state"] = obs["analysis_cell_id"].map(consensus_map)
    boundary = bool_series(obs["boundary"]) if "boundary" in obs else pd.Series(False, index=obs.index)
    condition = obs["condition"].astype(str) if "condition" in obs else pd.Series("", index=obs.index)
    obs["latent_cohort"] = np.select(
        [
            obs["consensus_state"].notna(),
            boundary & condition.eq("LAM"),
            ~condition.eq("LAM"),
        ],
        [
            obs["consensus_state"].map(lambda value: f"State_{value}"),
            "boundary",
            "normal",
        ],
        default="other",
    )
    keep = obs["latent_cohort"].ne("other")
    obs = obs.loc[keep].reset_index(drop=True)
    all_positions = latent_obj.obs_names.astype(str).tolist()
    position_map = {value: index for index, value in enumerate(all_positions)}
    positions = np.asarray([position_map[value] for value in obs["obs_name"].astype(str)], dtype=np.int64)
    latent = np.asarray(latent_obj.obsm["X_scVI"][positions], dtype=np.float32)
    target_mask = obs["latent_cohort"].eq(f"State_{TARGET_STATE}").to_numpy()
    if int(target_mask.sum()) == 0:
        latent_obj.file.close()
        return {"status": "target_state_missing", "path": str(scvi_path)}
    state15_latent = latent[target_mask]
    k = min(31, len(latent))
    neighbors = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(latent)
    distances, indices = neighbors.kneighbors(state15_latent)
    edge_rows: list[dict[str, Any]] = []
    for source_index in range(len(state15_latent)):
        source_id = obs.loc[obs.index[target_mask][source_index], "analysis_cell_id"]
        for rank in range(1, k):
            neighbor_row = obs.iloc[indices[source_index, rank]]
            edge_rows.append(
                {
                    "source_state15_cell": source_id,
                    "neighbor_rank": rank,
                    "neighbor_cell": neighbor_row["analysis_cell_id"],
                    "neighbor_cohort": neighbor_row["latent_cohort"],
                    "neighbor_state": neighbor_row["consensus_state"] if pd.notna(neighbor_row["consensus_state"]) else "",
                    "neighbor_dataset": neighbor_row.get("dataset", ""),
                    "neighbor_patient_id": neighbor_row.get("patient_id", ""),
                    "distance": float(distances[source_index, rank]),
                }
            )
    edges = pd.DataFrame(edge_rows)
    edges.to_csv(output_dir / "state15_latent_neighbor_edges.csv", index=False)
    neighbor_summary = (
        edges.groupby(["neighbor_cohort", "neighbor_state"], dropna=False, observed=True)
        .agg(
            source_state15_cells=("source_state15_cell", "nunique"),
            neighbor_edges=("neighbor_cell", "size"),
            mean_distance=("distance", "mean"),
            median_distance=("distance", "median"),
        )
        .reset_index()
    )
    neighbor_summary["edge_fraction"] = neighbor_summary["neighbor_edges"] / max(len(edges), 1)
    neighbor_summary.to_csv(output_dir / "state15_latent_neighbors.csv", index=False)

    nearest = NearestNeighbors(n_neighbors=1, metric="euclidean").fit(state15_latent)
    nearest_distance = nearest.kneighbors(latent, return_distance=True)[0].ravel()
    latent_cell = obs[["analysis_cell_id", "latent_cohort", "consensus_state"]].copy()
    for column in ["dataset", "patient_id", "condition", "boundary"]:
        if column in obs:
            latent_cell[column] = obs[column].to_numpy()
    latent_cell["nearest_state15_distance"] = nearest_distance
    score_indexed = score_table.set_index("analysis_cell_id")
    feature_columns = [
        column for column in [
            "LAMCORE_777", "melanocytic_identity", "contractile_marker_panel",
            "VEGFD_CTSK_panel", "CORE1", "CORE2", "CORE3_identity",
            "ciliated", "AT2", "macrophage", "endothelial", "fibroblast",
            "mesothelial", "pericyte_VSMC",
        ]
        if column in score_indexed
    ]
    for feature in feature_columns:
        latent_cell[feature] = latent_cell["analysis_cell_id"].map(score_indexed[feature])
    latent_cell.to_csv(output_dir / "state15_latent_distance_by_cell.csv", index=False)
    non_target = latent_cell[latent_cell["latent_cohort"] != f"State_{TARGET_STATE}"].copy()
    if len(non_target):
        non_target["distance_bin"] = pd.qcut(
            non_target["nearest_state15_distance"], q=10, duplicates="drop", labels=False
        )
    gradient_rows: list[dict[str, Any]] = []
    for bin_value, sub in non_target.groupby("distance_bin", observed=True):
        for cohort, cohort_sub in sub.groupby("latent_cohort", observed=True):
            row = {
                "distance_bin": int(bin_value),
                "cohort": str(cohort),
                "n_cells": len(cohort_sub),
                "distance_min": float(cohort_sub["nearest_state15_distance"].min()),
                "distance_median": float(cohort_sub["nearest_state15_distance"].median()),
                "distance_max": float(cohort_sub["nearest_state15_distance"].max()),
            }
            for feature in feature_columns:
                row[f"{feature}_mean"] = float(pd.to_numeric(cohort_sub[feature], errors="coerce").mean())
                row[f"{feature}_median"] = float(pd.to_numeric(cohort_sub[feature], errors="coerce").median())
            gradient_rows.append(row)
    pd.DataFrame(gradient_rows).to_csv(output_dir / "state15_latent_distance_gradient.csv", index=False)
    latent_obj.file.close()
    return {
        "status": "ok",
        "path": str(scvi_path),
        "n_latent_cells": int(len(obs)),
        "n_state15_cells": int(target_mask.sum()),
        "n_neighbors_per_state15_cell": int(k - 1),
        "state15_neighbor_summary": neighbor_summary.to_dict(orient="records"),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/state_modeling.yaml"))
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "results/stage18"))
    parser.add_argument("--block-size", type=int, default=2048)
    args = parser.parse_args()
    config = load_config(Path(args.config).resolve())
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    consensus_path = PROJECT_ROOT / str(config["outputs"]["consensus_upstream_annotations"])
    prepared_path = PROJECT_ROOT / str(config["outputs"]["prepared_h5ad"])
    if not consensus_path.exists():
        raise FileNotFoundError(f"Consensus annotation not found: {consensus_path}")
    if not prepared_path.exists():
        raise FileNotFoundError(f"Prepared AnnData not found: {prepared_path}")
    consensus = pd.read_csv(consensus_path)
    consensus["analysis_cell_id"] = consensus["analysis_cell_id"].astype(str)
    consensus["consensus_state"] = consensus["consensus_state"].astype(str)
    state15_ids = sorted(consensus.loc[consensus["consensus_state"].eq(TARGET_STATE), "analysis_cell_id"].tolist())
    if len(state15_ids) != 200:
        raise ValueError(f"Frozen State 15 expected 200 cells, found {len(state15_ids)}")
    state15_hash = hashlib.sha256("\n".join(state15_ids).encode("utf-8")).hexdigest()

    stage16_module = load_stage16_module()
    programs, program_manifest = load_programs(config)
    formal_genes, formal_manifest = resolve_formal_signature(config)
    competing = {name: list(genes) for name, genes in stage16_module.COMPETING_GENES.items()}
    modules: dict[str, list[str]] = {
        "melanocytic_identity": ["PMEL", "MLANA", "MITF"],
        "contractile_marker_panel": ["ACTA2", "ACTG2", "MYH11"],
        "VEGFD_CTSK_panel": ["VEGFD", "CTSK"],
        "ESR1_hormone_marker": ["ESR1"],
        "HOX_PBX_markers": ["EMX2", "HOXA11"],
    }
    for name in ["CORE1", "CORE2", "CORE3_identity", "LAM_myogenic_contractile", "ECM_remodeling", "mTOR_translation", "hormone_related", "hypoxia_stress", "protease_ECM_niche", "HOX_PBX", "normal_lung_interstitial"]:
        if name in programs:
            modules[name] = programs[name]
    if formal_genes:
        modules["LAMCORE_777"] = formal_genes
    modules.update(competing)
    modules = {name: unique_genes(genes) for name, genes in modules.items()}
    for gene in MARKER_GENES:
        modules[f"gene_{gene}"] = [gene]

    prepared = ad.read_h5ad(prepared_path, backed="r")
    obs = prepared_obs(prepared)
    obs_indexed = obs.set_index("analysis_cell_id", drop=False)
    consensus_meta = consensus.copy()
    consensus_meta["cohort"] = consensus_meta["consensus_state"].map(lambda value: f"State_{value}")
    consensus_meta["condition"] = consensus_meta.get("condition", "LAM")
    required_meta = ["analysis_cell_id", "dataset", "patient_id", "condition", "assay", "cell_type", "cohort", "consensus_state"]
    for column in required_meta:
        if column not in consensus_meta:
            consensus_meta[column] = obs_indexed.reindex(consensus_meta["analysis_cell_id"])[column].to_numpy() if column in obs_indexed else ""
    consensus_meta = consensus_meta[required_meta].copy()
    boundary = obs[(obs.get("condition", "").astype(str).eq("LAM")) & bool_series(obs.get("boundary", pd.Series(False, index=obs.index)))]
    boundary = boundary[~boundary["analysis_cell_id"].isin(consensus_meta["analysis_cell_id"])] if len(boundary) else boundary
    boundary_meta = boundary.copy()
    boundary_meta["cohort"] = "boundary"
    boundary_meta["consensus_state"] = ""
    normal = obs[~obs.get("condition", "").astype(str).eq("LAM")].copy()
    normal = normal[~normal["analysis_cell_id"].isin(consensus_meta["analysis_cell_id"])] if len(normal) else normal
    normal_meta = normal.copy()
    normal_meta["cohort"] = "normal"
    normal_meta["consensus_state"] = ""
    selected = pd.concat([consensus_meta, boundary_meta, normal_meta], ignore_index=True, sort=False)
    selected = selected.drop_duplicates("analysis_cell_id").reset_index(drop=True)

    score_values, score_manifest = selected_score_table(prepared, selected, modules, int(args.block_size))
    selected_scores = add_metadata_to_scores(score_values, selected)
    prepared.file.close()

    feature_columns = [name for name in modules if name in selected_scores.columns]
    profile = summary_rows(selected_scores, feature_columns)
    profile.to_csv(output_dir / "state15_marker_profile.csv", index=False)
    lamcore_profile = profile[profile["feature"].eq("LAMCORE_777")].copy() if formal_genes else pd.DataFrame()
    lamcore_profile.to_csv(output_dir / "state15_lamcore_summary.csv", index=False)

    enrichment = fisher_enrichment(consensus)
    enrichment.to_csv(output_dir / "state15_author_enrichment.csv", index=False)
    consensus_score_table = selected_scores[selected_scores["cohort"].astype(str).str.startswith("State_")].copy()
    comparison_features = [
            name for name in [
            "gene_PMEL", "gene_MLANA", "gene_MITF", "gene_ACTA2", "gene_ACTG2",
            "gene_MYH11", "gene_VEGFD", "gene_CTSK", "gene_EMX2", "gene_HOXA11", "gene_ESR1",
            "melanocytic_identity", "contractile_marker_panel", "VEGFD_CTSK_panel", "ESR1_hormone_marker",
            "CORE1", "CORE2", "CORE3_identity", "HOX_PBX", "LAMCORE_777",
            "ECM_remodeling", "mTOR_translation", "hormone_related", "hypoxia_stress", "protease_ECM_niche",
            "ciliated", "AT2", "macrophage", "endothelial", "fibroblast", "mesothelial", "pericyte_VSMC",
        ]
        if name in consensus_score_table.columns
    ]
    comparisons = compare_states(consensus_score_table, comparison_features)
    comparisons.to_csv(output_dir / "state15_vs_comparators.csv", index=False)

    patient_bulk = ad.read_h5ad(prepared_path, backed="r")
    patient_pseudobulk_table = patient_pseudobulk(patient_bulk, consensus, modules)
    patient_bulk.file.close()
    patient_pseudobulk_table.to_csv(output_dir / "state15_patient_pseudobulk.csv", index=False)

    patient_rows: list[dict[str, Any]] = []
    all_patients = sorted(consensus["patient_id"].astype(str).unique())
    state15_scores = consensus_score_table[consensus_score_table["consensus_state"].eq(TARGET_STATE)]
    for patient in all_patients:
        sub = state15_scores[state15_scores["patient_id"].astype(str).eq(patient)]
        row: dict[str, Any] = {
            "patient_id": patient,
            "state15_present": bool(len(sub)),
            "state15_cells": int(len(sub)),
            "state15_fraction_of_consensus": float(len(sub) / max((consensus["patient_id"].astype(str) == patient).sum(), 1)),
            "dataset_count": int(sub["dataset"].astype(str).nunique()) if len(sub) else 0,
            "author_style_cells": int(bool_series(sub["source_author_style"]).sum()) if len(sub) and "source_author_style" in sub else 0,
        }
        for feature in [
            "LAMCORE_777", "gene_PMEL", "gene_MLANA", "gene_MITF", "gene_VEGFD", "gene_CTSK",
            "CORE1", "CORE2", "CORE3_identity", "contractile_marker_panel",
        ]:
            if feature in sub:
                row[f"{feature}_mean"] = float(pd.to_numeric(sub[feature], errors="coerce").mean())
                row[f"{feature}_median"] = float(pd.to_numeric(sub[feature], errors="coerce").median())
        patient_rows.append(row)
    patient_consistency = pd.DataFrame(patient_rows)
    patient_consistency.to_csv(output_dir / "state15_patient_consistency.csv", index=False)

    target_lamcore = profile[(profile["cohort"] == "State_15") & (profile["feature"] == "LAMCORE_777")]
    boundary_lamcore = profile[(profile["cohort"] == "boundary") & (profile["feature"] == "LAMCORE_777")]
    normal_lamcore = profile[(profile["cohort"] == "normal") & (profile["feature"] == "LAMCORE_777")]
    comparator_lamcore = profile[
        profile["cohort"].isin([f"State_{state}" for state in COMPARATORS])
        & profile["feature"].eq("LAMCORE_777")
    ]
    target_median = float(target_lamcore["median"].iloc[0]) if len(target_lamcore) else np.nan
    boundary_median = float(boundary_lamcore["median"].iloc[0]) if len(boundary_lamcore) else np.nan
    normal_median = float(normal_lamcore["median"].iloc[0]) if len(normal_lamcore) else np.nan
    comparator_medians = comparator_lamcore.set_index("cohort")["median"].to_dict()
    overall_author = enrichment[
        enrichment["stratum_level"].eq("overall") & enrichment["stratum"].eq("all")
    ].iloc[0]

    existing_stability = PROJECT_ROOT / "results/stage7/state_stability_summary.csv"
    existing_reproducibility = PROJECT_ROOT / "results/stage11/state_reproducibility_summary.csv"
    stability_metrics: dict[str, Any] = {}
    if existing_stability.exists():
        stability_table = pd.read_csv(existing_stability)
        stability_row = stability_table[stability_table["consensus_state"].astype(str).eq(TARGET_STATE)]
        if len(stability_row):
            stability_metrics = {
                "mean_within_coassignment": float(stability_row["mean_within_coassignment"].iloc[0]),
                "mean_margin": float(stability_row["mean_margin"].iloc[0]),
            }
    if existing_reproducibility.exists():
        reproducibility_table = pd.read_csv(existing_reproducibility)
        reproducibility_row = reproducibility_table[reproducibility_table["state_id"].astype(str).eq(TARGET_STATE)]
        if len(reproducibility_row):
            for column in [
                "patient_coverage", "dataset_coverage", "structural_stability",
                "biological_reproducibility", "mean_loo_recovery", "mean_loo_additional_loss",
            ]:
                if column in reproducibility_row:
                    value = pd.to_numeric(reproducibility_row[column], errors="coerce").iloc[0]
                    stability_metrics[column] = float(value) if pd.notna(value) else np.nan

    patient_count_with_state15 = int(patient_consistency["state15_present"].sum())
    dominant_patient_cells = int(patient_consistency["state15_cells"].max()) if len(patient_consistency) else 0
    anchor_decision = "provisional_reference_candidate_not_formally_upgraded"
    anchor_decision_basis = {
        "formal_lamcore_elevated_vs_boundary": bool(np.isfinite(target_median) and np.isfinite(boundary_median) and target_median > boundary_median),
        "formal_lamcore_elevated_vs_normal": bool(np.isfinite(target_median) and np.isfinite(normal_median) and target_median > normal_median),
        "formal_lamcore_above_all_requested_comparators": bool(
            len(comparator_medians) == len(COMPARATORS)
            and all(target_median > float(value) for value in comparator_medians.values())
        ),
        "author_style_overall_enriched": bool(float(overall_author["enrichment_fold"]) > 1 and float(overall_author["fisher_pvalue_greater"]) < 0.05),
        "state15_patient_coverage": patient_count_with_state15,
        "consensus_patient_count": int(len(all_patients)),
        "dominant_patient": str(patient_consistency.loc[patient_consistency["state15_cells"].idxmax(), "patient_id"]) if len(patient_consistency) else "",
        "dominant_patient_cells": dominant_patient_cells,
        "dominant_patient_fraction_of_state15": float(dominant_patient_cells / max(len(state15_ids), 1)),
        "author_style_support_datasets": sorted(
            consensus.loc[
                bool_series(consensus["source_author_style"]) & consensus["consensus_state"].eq(TARGET_STATE), "dataset"
            ].astype(str).unique()
        ) if "source_author_style" in consensus else [],
        "existing_consensus_and_loo_metrics": stability_metrics,
        "decision_note": "Formal score and comparator separation support State 15 as a useful provisional reference, but patient concentration and author-label concentration prevent formal upgrade to a cross-patient anchor in this validation stage.",
    }

    latent_manifest = latent_outputs(config, consensus, selected_scores, output_dir)
    anchor_summary = {
        "target_state": TARGET_STATE,
        "target_cell_count": len(state15_ids),
        "frozen_state15_cell_id_sha256": state15_hash,
        "consensus_cell_count": len(consensus),
        "state15_dataset_count": int(consensus.loc[consensus["consensus_state"].eq(TARGET_STATE), "dataset"].astype(str).nunique()),
        "state15_patient_count": int(consensus.loc[consensus["consensus_state"].eq(TARGET_STATE), "patient_id"].astype(str).nunique()),
        "state15_dataset_counts": consensus.loc[consensus["consensus_state"].eq(TARGET_STATE), "dataset"].astype(str).value_counts().to_dict(),
        "state15_patient_counts": consensus.loc[consensus["consensus_state"].eq(TARGET_STATE), "patient_id"].astype(str).value_counts().to_dict(),
        "formal_signature_manifest": formal_manifest,
        "program_manifest": program_manifest,
        "score_manifest": score_manifest,
        "latent_manifest": latent_manifest,
        "no_reclustering": True,
        "no_scvi_training": True,
        "formal_score_is_independent_of_stage16": True,
        "outputs_are_diagnostic_only": True,
        "anchor_decision": anchor_decision,
        "anchor_decision_basis": anchor_decision_basis,
    }
    (output_dir / "state15_anchor_summary.json").write_text(json.dumps(anchor_summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    overall = enrichment[(enrichment["stratum_level"] == "overall") & (enrichment["stratum"] == "all")].iloc[0]
    target_profile = target_lamcore
    normal_profile = normal_lamcore
    boundary_profile = boundary_lamcore
    report = [
        "# State 15 LAM-core reference anchor validation",
        "",
        "State 15 was frozen from the existing consensus annotation. This stage did not recluster, retrain scVI, or modify any candidate gate.",
        "",
        f"- Frozen State 15 cells: {len(state15_ids)}",
        f"- Frozen State 15 ID SHA-256: {state15_hash}",
        f"- State 15 dataset coverage: {anchor_summary['state15_dataset_count']} datasets",
        f"- State 15 patient coverage: {anchor_summary['state15_patient_count']} patients",
        f"- Formal signature: {formal_manifest.get('status')} ({formal_manifest.get('n_genes', 0)} genes)",
        "",
        "## Formal LAMCORE and author-label evidence",
        "",
    ]
    if len(target_profile):
        report.append(f"- State 15 LAMCORE median: {float(target_profile['median'].iloc[0]):.4f}")
    if len(boundary_profile):
        report.append(f"- Boundary LAMCORE median: {float(boundary_profile['median'].iloc[0]):.4f}")
    if len(normal_profile):
        report.append(f"- Normal/control LAMCORE median: {float(normal_profile['median'].iloc[0]):.4f}")
    report.extend(
        [
            f"- Overall author-style enrichment fold: {float(overall['enrichment_fold']):.4f}",
            f"- Overall Fisher exact one-sided p-value: {float(overall['fisher_pvalue_greater']):.4g}",
            "",
            "## Anchor decision",
            "",
            f"- Decision: `{anchor_decision}`",
            f"- State 15 is present in {patient_count_with_state15}/{len(all_patients)} consensus patients; the largest contribution is {anchor_decision_basis['dominant_patient']} ({dominant_patient_cells}/{len(state15_ids)} cells).",
            f"- Author-style support is present in datasets: {', '.join(anchor_decision_basis['author_style_support_datasets']) or 'none'}.",
            "- Formal LAMCORE/comparator separation is supportive, but the current patient and author-label concentration is not sufficient to promote State 15 to a formally cross-patient reference anchor.",
            "- State 15 remains frozen as a provisional reference candidate for a later, explicitly designed expansion analysis; no gate or existing state artifact is changed here.",
            "",
            "## Interpretation",
            "",
            "The tables retain State 15 as a fixed validation object. Formal LAMCORE elevation, author-label enrichment, patient-level consistency, comparator separation, and latent-neighborhood continuity must be considered together; no single score or cell count is promoted to an automatic gate.",
            "",
            "## Outputs",
            "",
            "- state15_lamcore_summary.csv",
            "- state15_marker_profile.csv",
            "- state15_author_enrichment.csv",
            "- state15_vs_comparators.csv",
            "- state15_patient_pseudobulk.csv",
            "- state15_patient_consistency.csv",
            "- state15_latent_neighbors.csv",
            "- state15_latent_neighbor_edges.csv",
            "- state15_latent_distance_by_cell.csv",
            "- state15_latent_distance_gradient.csv",
        ]
    )
    (output_dir / "state15_anchor_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Frozen State 15: {len(state15_ids)} cells")
    print(f"Scored selected cohorts: {len(selected_scores)} cells")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
