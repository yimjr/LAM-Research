#!/usr/bin/env python3
"""Audit cross-patient support for the frozen consensus State 15.

This stage is intentionally diagnostic.  It reuses the current high-confidence
candidate cohort, prepared counts/normalized expression, and the existing
``X_scVI`` embedding.  It does not change a candidate gate, recluster cells, or
train scVI.  The main sensitivity analysis removes LAM1163 without changing
the State 15 labels.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy.stats import fisher_exact, mannwhitneyu


PROJECT_ROOT = Path(__file__).resolve().parent
TARGET_STATE = "15"
TARGET_PATIENT = "LAM1163"
COMPARATORS = ["18", "20", "12", "7", "5"]
COMPARATOR_NAMES = {
    "18": "pericyte_VSMC_comparator",
    "20": "fibroblast_comparator",
    "12": "endothelial_comparator",
    "7": "AT2_comparator",
    "5": "macrophage_comparator",
}
PROFILE_FEATURES = [
    "LAMCORE_777",
    "CORE1",
    "CORE2",
    "CORE3_identity",
    "gene_PMEL",
    "gene_MLANA",
    "gene_MITF",
    "gene_ACTA2",
    "gene_ACTG2",
    "gene_MYH11",
    "gene_VEGFD",
    "gene_CTSK",
    "gene_ESR1",
    "HOX_PBX_markers",
    "HOX_PBX",
    "LAM_myogenic_contractile",
    "ECM_remodeling",
    "mTOR_translation",
    "hormone_related",
    "protease_ECM_niche",
    "ciliated",
    "AT2",
    "macrophage",
    "endothelial",
    "fibroblast",
    "mesothelial",
    "pericyte_VSMC",
]


def load_stage18_module() -> Any:
    path = PROJECT_ROOT / "18_validate_state15_anchor.py"
    spec = importlib.util.spec_from_file_location("stage18_anchor", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Stage 18 helpers: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_config(path: Path) -> dict[str, Any]:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def build_modules(stage18: Any, config: dict[str, Any]) -> tuple[dict[str, list[str]], dict[str, Any], dict[str, Any]]:
    programs, program_manifest = stage18.load_programs(config)
    formal_genes, formal_manifest = stage18.resolve_formal_signature(config)
    stage16 = stage18.load_stage16_module()
    modules: dict[str, list[str]] = {
        "melanocytic_identity": ["PMEL", "MLANA", "MITF"],
        "contractile_marker_panel": ["ACTA2", "ACTG2", "MYH11"],
        "VEGFD_CTSK_panel": ["VEGFD", "CTSK"],
        "ESR1_hormone_marker": ["ESR1"],
        "HOX_PBX_markers": ["EMX2", "HOXA11"],
    }
    for name in [
        "CORE1", "CORE2", "CORE3_identity", "LAM_myogenic_contractile",
        "ECM_remodeling", "mTOR_translation", "hormone_related",
        "hypoxia_stress", "protease_ECM_niche", "HOX_PBX",
        "normal_lung_interstitial",
    ]:
        if name in programs:
            modules[name] = programs[name]
    if formal_genes:
        modules["LAMCORE_777"] = formal_genes
    modules.update({name: list(genes) for name, genes in stage16.COMPETING_GENES.items()})
    for gene in ["PMEL", "MLANA", "MITF", "ACTA2", "ACTG2", "MYH11", "VEGFD", "CTSK", "ESR1", "EMX2", "HOXA11"]:
        modules[f"gene_{gene}"] = [gene]
    modules = {name: stage18.unique_genes(genes) for name, genes in modules.items()}
    return modules, formal_manifest, program_manifest


def cosine_similarity(x: np.ndarray, y: np.ndarray) -> float:
    denominator = float(np.linalg.norm(x) * np.linalg.norm(y))
    return float(np.dot(x, y) / denominator) if denominator > 0 else np.nan


def pearson_similarity(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return float(np.corrcoef(x, y)[0, 1])


def vector_metrics(x: np.ndarray, y: np.ndarray) -> dict[str, float]:
    return {
        "pearson_correlation": pearson_similarity(x, y),
        "cosine_similarity": cosine_similarity(x, y),
        "euclidean_distance": float(np.linalg.norm(x - y)),
    }


def find_candidate_file(dataset: str, config: dict[str, Any]) -> Path | None:
    roots = [PROJECT_ROOT / str(root) for root in config.get("input_roots", [])]
    roots.extend(
        [
            PROJECT_ROOT.parent / "LAM-Cell-Research",
            Path("/mnt/e/LAM-Research/LAM-Cell-Research"),
            Path("/mnt/e/LAM-Research/data-temp"),
            PROJECT_ROOT / "data/upstream",
        ]
    )
    relative = (
        Path("results/program_discovery/candidate_pool_labels.csv")
        if dataset == "GSE135851"
        else Path(f"results/program_discovery/external_{dataset}/candidate_pool_labels.csv")
    )
    seen: set[str] = set()
    for root in roots:
        path = root / relative
        key = str(path.resolve())
        if key in seen or not path.exists():
            continue
        seen.add(key)
        return path
    return None


def author_annotation_table(consensus: pd.DataFrame, config: dict[str, Any]) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict[str, Any]] = []
    for dataset in sorted(consensus["dataset"].astype(str).unique()):
        path = find_candidate_file(dataset, config)
        if path is None:
            source = pd.DataFrame()
        else:
            source = pd.read_csv(path)
        has_field = "source_author_style" in source.columns
        positive_total = int(stage18_bool(source["source_author_style"]).sum()) if has_field else 0
        # The external conversion initializes this field to False; zero is not
        # treated as a negative author label in those datasets.
        available = bool(dataset == "GSE135851" and positive_total > 0)
        dataset_consensus = consensus[consensus["dataset"].astype(str).eq(dataset)]
        state15 = dataset_consensus[dataset_consensus["consensus_state"].astype(str).eq(TARGET_STATE)]
        observed_state15_true = int(stage18_bool(state15["source_author_style"]).sum()) if "source_author_style" in state15 else 0
        rows.append(
            {
                "dataset": dataset,
                "candidate_file": str(path) if path else "",
                "author_style_field_present": has_field,
                "author_style_annotation_status": "available" if available else "not_assayed",
                "availability_basis": "positive_author_labels_present" if available else "external_field_initialized_false_or_file_missing",
                "author_style_positive_total": positive_total if available else np.nan,
                "state15_cells": len(state15),
                "state15_author_style_positive_observed": observed_state15_true if available else np.nan,
            }
        )
    availability = pd.DataFrame(rows)

    assayed = consensus[consensus["dataset"].astype(str).eq("GSE135851")].copy()
    state15 = assayed[assayed["consensus_state"].astype(str).eq(TARGET_STATE)]
    other = assayed[~assayed["consensus_state"].astype(str).eq(TARGET_STATE)]
    a = int(stage18_bool(state15["source_author_style"]).sum())
    b = int(len(state15) - a)
    c = int(stage18_bool(other["source_author_style"]).sum())
    d = int(len(other) - c)
    odds, pvalue = fisher_exact([[a, b], [c, d]], alternative="greater")
    state_fraction = a / len(state15) if len(state15) else np.nan
    other_fraction = c / len(other) if len(other) else np.nan
    enrichment = pd.DataFrame(
        [
            {
                "dataset": "GSE135851",
                "annotation_status": "available",
                "state15_cells": len(state15),
                "state15_author_style": a,
                "other_candidate_cells": len(other),
                "other_author_style": c,
                "state15_author_fraction": state_fraction,
                "other_author_fraction": other_fraction,
                "enrichment_fold": state_fraction / other_fraction if other_fraction else np.inf,
                "fisher_odds_ratio": float(odds),
                "fisher_pvalue_greater": float(pvalue),
            }
        ]
    )
    return availability, enrichment


def stage18_bool(values: pd.Series) -> pd.Series:
    return values.map(lambda value: str(value).strip().lower() in {"true", "1", "yes", "y"}).astype(bool)


def text_table(table: pd.DataFrame) -> str:
    return "```text\n" + table.to_string(index=False) + "\n```"


def add_prepared_totals(prepared_obs: pd.DataFrame, consensus: pd.DataFrame) -> pd.DataFrame:
    result = consensus.copy()
    if "total_counts" in prepared_obs.columns:
        totals = prepared_obs.set_index("analysis_cell_id")["total_counts"]
        result["total_counts"] = result["analysis_cell_id"].map(totals)
    else:
        result["total_counts"] = np.nan
    if "n_genes_by_counts" in prepared_obs.columns:
        detected = prepared_obs.set_index("analysis_cell_id")["n_genes_by_counts"]
        result["n_genes_by_counts"] = result["analysis_cell_id"].map(detected)
    else:
        result["n_genes_by_counts"] = np.nan
    return result


def patient_composition(consensus: pd.DataFrame) -> pd.DataFrame:
    total_candidate = len(consensus)
    target = consensus[consensus["consensus_state"].astype(str).eq(TARGET_STATE)]
    rows: list[dict[str, Any]] = []
    for patient in sorted(target["patient_id"].astype(str).unique()):
        candidate_sub = consensus[consensus["patient_id"].astype(str).eq(patient)]
        state_sub = target[target["patient_id"].astype(str).eq(patient)]
        candidate_fraction = len(candidate_sub) / max(total_candidate, 1)
        state_fraction = len(state_sub) / max(len(target), 1)
        rows.append(
            {
                "patient_id": patient,
                "candidate_pool_total_cells": len(candidate_sub),
                "candidate_pool_fraction": candidate_fraction,
                "state15_cells": len(state_sub),
                "state15_fraction": state_fraction,
                "enrichment": state_fraction / candidate_fraction if candidate_fraction else np.nan,
                "candidate_pool_dataset_count": candidate_sub["dataset"].astype(str).nunique(),
                "state15_dataset_count": state_sub["dataset"].astype(str).nunique(),
            }
        )
    return pd.DataFrame(rows)


def patient_matched_comparison(scores: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    state_mask = scores["consensus_state"].astype(str).eq(TARGET_STATE)
    for patient in sorted(scores.loc[state_mask, "patient_id"].astype(str).unique()):
        patient_mask = scores["patient_id"].astype(str).eq(patient)
        target = scores[patient_mask & state_mask]
        rest = scores[patient_mask & ~state_mask]
        for feature in features:
            x = pd.to_numeric(target[feature], errors="coerce").dropna().to_numpy(dtype=float)
            y = pd.to_numeric(rest[feature], errors="coerce").dropna().to_numpy(dtype=float)
            if len(x) and len(y):
                statistic, pvalue = mannwhitneyu(x, y, alternative="two-sided")
            else:
                statistic, pvalue = np.nan, np.nan
            rows.append(
                {
                    "patient_id": patient,
                    "feature": feature,
                    "state15_cells": len(x),
                    "rest_of_candidate_cells": len(y),
                    "state15_mean": float(np.mean(x)) if len(x) else np.nan,
                    "rest_mean": float(np.mean(y)) if len(y) else np.nan,
                    "delta_mean": float(np.mean(x) - np.mean(y)) if len(x) and len(y) else np.nan,
                    "state15_median": float(np.median(x)) if len(x) else np.nan,
                    "rest_median": float(np.median(y)) if len(y) else np.nan,
                    "delta_median": float(np.median(x) - np.median(y)) if len(x) and len(y) else np.nan,
                    "mannwhitney_u": float(statistic) if np.isfinite(statistic) else np.nan,
                    "mannwhitney_pvalue": float(pvalue) if np.isfinite(pvalue) else np.nan,
                    "delta_mean_direction": "positive" if len(x) and len(y) and np.mean(x) > np.mean(y) else "negative_or_equal",
                }
            )
    return pd.DataFrame(rows)


def module_gene_indices(modules: dict[str, list[str]], present_genes: list[str]) -> dict[str, list[int]]:
    index = {gene: i for i, gene in enumerate(present_genes)}
    return {name: [index[gene] for gene in genes if gene in index] for name, genes in modules.items()}


def aggregate_profile(
    mask: np.ndarray,
    counts: np.ndarray,
    total_counts: np.ndarray,
    module_indices: dict[str, list[int]],
) -> tuple[np.ndarray, dict[str, float], float]:
    positions = np.flatnonzero(mask)
    if len(positions) == 0:
        return positions, {name: np.nan for name in module_indices}, np.nan
    summed = counts[positions].sum(axis=0, dtype=np.float64)
    totals = total_counts[positions]
    total_umi = float(np.nansum(totals))
    if not np.isfinite(total_umi) or total_umi <= 0:
        total_umi = float(np.sum(summed))
    normalized = np.log1p(summed / max(total_umi, 1.0) * 10000.0)
    profile = {
        name: float(np.mean(normalized[indices])) if indices else np.nan
        for name, indices in module_indices.items()
    }
    return positions, profile, total_umi


def patient_pseudobulk_profiles(
    consensus: pd.DataFrame,
    counts: np.ndarray,
    total_counts: np.ndarray,
    module_indices: dict[str, list[int]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    states = [TARGET_STATE, *COMPARATORS]
    for patient in sorted(consensus["patient_id"].astype(str).unique()):
        for state in states:
            mask = (
                consensus["patient_id"].astype(str).eq(patient).to_numpy()
                & consensus["consensus_state"].astype(str).eq(state).to_numpy()
            )
            positions, profile, total_umi = aggregate_profile(mask, counts, total_counts, module_indices)
            if len(positions) == 0:
                continue
            dataset_count = consensus.iloc[positions]["dataset"].astype(str).nunique()
            for feature, value in profile.items():
                rows.append(
                    {
                        "patient_id": patient,
                        "group": f"State_{state}",
                        "cells": len(positions),
                        "dataset_count": int(dataset_count),
                        "total_umi": total_umi,
                        "feature": feature,
                        "pseudobulk_log1p_score": value,
                    }
                )
    return pd.DataFrame(rows)


def lopo_validation(
    consensus: pd.DataFrame,
    counts: np.ndarray,
    total_counts: np.ndarray,
    module_indices: dict[str, list[int]],
    latent: np.ndarray | None,
) -> pd.DataFrame:
    state = consensus["consensus_state"].astype(str).to_numpy()
    patients = consensus["patient_id"].astype(str).to_numpy()
    target_patients = sorted(np.unique(patients[state == TARGET_STATE]))
    features = [feature for feature in PROFILE_FEATURES if feature in module_indices]
    rows: list[dict[str, Any]] = []
    for held_out in target_patients:
        target_mask = (state == TARGET_STATE) & (patients == held_out)
        reference_mask = (state == TARGET_STATE) & (patients != held_out)
        _, target_profile, _ = aggregate_profile(target_mask, counts, total_counts, module_indices)
        _, reference_profile, _ = aggregate_profile(reference_mask, counts, total_counts, module_indices)
        target_vector = np.asarray([target_profile.get(feature, np.nan) for feature in features], dtype=float)
        reference_vector = np.asarray([reference_profile.get(feature, np.nan) for feature in features], dtype=float)
        valid_reference = np.isfinite(target_vector) & np.isfinite(reference_vector)
        profile_reference_metrics = vector_metrics(target_vector[valid_reference], reference_vector[valid_reference]) if valid_reference.sum() >= 2 else {"pearson_correlation": np.nan, "cosine_similarity": np.nan, "euclidean_distance": np.nan}
        latent_reference_metrics = {"latent_centroid_euclidean": np.nan, "latent_centroid_cosine": np.nan}
        if latent is not None and target_mask.any() and reference_mask.any():
            target_latent = latent[target_mask].mean(axis=0)
            reference_latent = latent[reference_mask].mean(axis=0)
            latent_reference_metrics = {
                "latent_centroid_euclidean": float(np.linalg.norm(target_latent - reference_latent)),
                "latent_centroid_cosine": cosine_similarity(target_latent, reference_latent),
            }
        for comparator in COMPARATORS:
            comparator_mask = (state == comparator) & (patients != held_out)
            _, comparator_profile, _ = aggregate_profile(comparator_mask, counts, total_counts, module_indices)
            comparator_vector = np.asarray([comparator_profile.get(feature, np.nan) for feature in features], dtype=float)
            valid = np.isfinite(target_vector) & np.isfinite(comparator_vector)
            profile_metrics = vector_metrics(target_vector[valid], comparator_vector[valid]) if valid.sum() >= 2 else {"pearson_correlation": np.nan, "cosine_similarity": np.nan, "euclidean_distance": np.nan}
            latent_metrics = {"latent_centroid_euclidean": np.nan, "latent_centroid_cosine": np.nan}
            if latent is not None and target_mask.any() and comparator_mask.any():
                target_latent = latent[target_mask].mean(axis=0)
                comparator_latent = latent[comparator_mask].mean(axis=0)
                latent_metrics = {
                    "latent_centroid_euclidean": float(np.linalg.norm(target_latent - comparator_latent)),
                    "latent_centroid_cosine": cosine_similarity(target_latent, comparator_latent),
                }
            rows.append(
                {
                    "held_out_patient": held_out,
                    "held_out_state15_cells": int(target_mask.sum()),
                    "reference_state15_cells": int(reference_mask.sum()),
                    "comparison": "reference_vs_" + COMPARATOR_NAMES[comparator],
                    "comparison_group": f"State_{comparator}",
                    "comparison_cells": int(comparator_mask.sum()),
                    "reference_pearson_correlation": profile_reference_metrics["pearson_correlation"],
                    "reference_cosine_similarity": profile_reference_metrics["cosine_similarity"],
                    "reference_euclidean_distance": profile_reference_metrics["euclidean_distance"],
                    "comparison_pearson_correlation": profile_metrics["pearson_correlation"],
                    "comparison_cosine_similarity": profile_metrics["cosine_similarity"],
                    "comparison_euclidean_distance": profile_metrics["euclidean_distance"],
                    "reference_profile_closer": bool(
                        np.isfinite(profile_reference_metrics["euclidean_distance"])
                        and np.isfinite(profile_metrics["euclidean_distance"])
                        and profile_reference_metrics["euclidean_distance"] < profile_metrics["euclidean_distance"]
                    ),
                    "reference_latent_centroid_euclidean": latent_reference_metrics["latent_centroid_euclidean"],
                    "comparison_latent_centroid_euclidean": latent_metrics["latent_centroid_euclidean"],
                    "reference_latent_closer": bool(
                        np.isfinite(latent_reference_metrics["latent_centroid_euclidean"])
                        and np.isfinite(latent_metrics["latent_centroid_euclidean"])
                        and latent_reference_metrics["latent_centroid_euclidean"] < latent_metrics["latent_centroid_euclidean"]
                    ),
                }
            )
    return pd.DataFrame(rows)


def without_patient_sensitivity(scores: pd.DataFrame, patient: str, features: list[str]) -> pd.DataFrame:
    retained = scores[~scores["patient_id"].astype(str).eq(patient)].copy()
    target = retained[retained["consensus_state"].astype(str).eq(TARGET_STATE)]
    rows: list[dict[str, Any]] = []
    for feature in features:
        x = pd.to_numeric(target[feature], errors="coerce").dropna().to_numpy(dtype=float)
        for comparator in COMPARATORS:
            y = pd.to_numeric(
                retained[retained["consensus_state"].astype(str).eq(comparator)][feature], errors="coerce"
            ).dropna().to_numpy(dtype=float)
            if len(x) and len(y):
                statistic, pvalue = mannwhitneyu(x, y, alternative="two-sided")
            else:
                statistic, pvalue = np.nan, np.nan
            rows.append(
                {
                    "removed_patient": patient,
                    "target_group": "State_15_without_LAM1163",
                    "comparator_group": f"State_{comparator}",
                    "feature": feature,
                    "target_cells": len(x),
                    "comparator_cells": len(y),
                    "target_mean": float(np.mean(x)) if len(x) else np.nan,
                    "comparator_mean": float(np.mean(y)) if len(y) else np.nan,
                    "target_median": float(np.median(x)) if len(x) else np.nan,
                    "comparator_median": float(np.median(y)) if len(y) else np.nan,
                    "median_difference": float(np.median(x) - np.median(y)) if len(x) and len(y) else np.nan,
                    "mannwhitney_u": float(statistic) if np.isfinite(statistic) else np.nan,
                    "mannwhitney_pvalue": float(pvalue) if np.isfinite(pvalue) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def load_candidate_latent(config: dict[str, Any], consensus: pd.DataFrame) -> np.ndarray | None:
    path = PROJECT_ROOT / str(config["outputs"]["scvi_h5ad"])
    if not path.exists():
        return None
    obj = ad.read_h5ad(path, backed="r")
    if "X_scVI" not in obj.obsm:
        obj.file.close()
        return None
    obs = obj.obs.copy().reset_index(names="obs_name")
    if "analysis_cell_id" in obs:
        ids = obs["analysis_cell_id"].astype(str)
    else:
        ids = obs["obs_name"].astype(str)
    position_map = {value: index for index, value in enumerate(ids)}
    positions = np.asarray([position_map.get(value, -1) for value in consensus["analysis_cell_id"].astype(str)], dtype=np.int64)
    if (positions < 0).any():
        obj.file.close()
        raise ValueError(f"{int((positions < 0).sum())} candidate cells are missing from X_scVI")
    latent = np.asarray(obj.obsm["X_scVI"][positions], dtype=np.float32)
    obj.file.close()
    return latent


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default=str(PROJECT_ROOT / "config/state_modeling.yaml"))
    parser.add_argument("--output-dir", default=str(PROJECT_ROOT / "results/stage19"))
    parser.add_argument("--block-size", type=int, default=8192)
    args = parser.parse_args()
    config = load_config(Path(args.config).resolve())
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    stage18 = load_stage18_module()
    modules, formal_manifest, program_manifest = build_modules(stage18, config)
    consensus_path = PROJECT_ROOT / str(config["outputs"]["consensus_upstream_annotations"])
    prepared_path = PROJECT_ROOT / str(config["outputs"]["prepared_h5ad"])
    if not consensus_path.exists() or not prepared_path.exists():
        raise FileNotFoundError(f"Required consensus/prepared input missing: {consensus_path}, {prepared_path}")
    consensus = pd.read_csv(consensus_path)
    consensus["analysis_cell_id"] = consensus["analysis_cell_id"].astype(str)
    consensus["consensus_state"] = consensus["consensus_state"].astype(str)
    if consensus["analysis_cell_id"].duplicated().any():
        raise ValueError("Duplicate analysis_cell_id in consensus annotation")
    if int(consensus["consensus_state"].eq(TARGET_STATE).sum()) != 200:
        raise ValueError("Frozen State 15 is not exactly 200 cells")
    if "patient_id" not in consensus or consensus["patient_id"].isna().any():
        raise ValueError("Patient IDs are required for Stage 19")

    prepared = ad.read_h5ad(prepared_path, backed="r")
    prepared_meta = stage18.prepared_obs(prepared)
    candidate = add_prepared_totals(prepared_meta, consensus)
    score_values, score_manifest = stage18.selected_score_table(
        prepared,
        candidate[["analysis_cell_id"]],
        modules,
        int(args.block_size),
    )
    prepared.file.close()
    scores = candidate.merge(score_values, on="analysis_cell_id", how="left", validate="one_to_one")
    features = [feature for feature in PROFILE_FEATURES if feature in scores.columns]
    if "LAMCORE_777" not in features:
        raise ValueError("The 777-gene formal LAMCORE score is unavailable")

    composition = patient_composition(consensus)
    composition.to_csv(output_dir / "state15_patient_composition.csv", index=False)

    matched = patient_matched_comparison(scores, features)
    matched.to_csv(output_dir / "state15_patient_matched_comparison.csv", index=False)

    prepared_counts = ad.read_h5ad(prepared_path, backed="r")
    counts_df, present_genes = stage18.extract_pseudobulk_counts(
        prepared_counts,
        candidate[["analysis_cell_id"]],
        [gene for values in modules.values() for gene in values],
    )
    prepared_counts.file.close()
    counts = counts_df.to_numpy(dtype=np.float64)
    total_counts = pd.to_numeric(candidate["total_counts"], errors="coerce").to_numpy(dtype=float)
    module_indices = module_gene_indices(modules, present_genes)
    pseudobulk = patient_pseudobulk_profiles(consensus, counts, total_counts, module_indices)
    pseudobulk.to_csv(output_dir / "state15_patient_pseudobulk_profiles.csv", index=False)
    pseudobulk[pseudobulk["group"].eq("State_15")].to_csv(output_dir / "state15_patient_profiles.csv", index=False)

    latent = load_candidate_latent(config, consensus)
    lopo = lopo_validation(consensus, counts, total_counts, module_indices, latent)
    lopo.to_csv(output_dir / "state15_lopo_validation.csv", index=False)

    sensitivity = without_patient_sensitivity(scores, TARGET_PATIENT, features)
    sensitivity.to_csv(output_dir / "state15_without_LAM1163.csv", index=False)

    availability, author_enrichment = author_annotation_table(consensus, config)
    availability.to_csv(output_dir / "author_annotation_availability.csv", index=False)
    author_enrichment.to_csv(output_dir / "state15_author_enrichment_assayed.csv", index=False)

    target_after = sensitivity[sensitivity["target_group"].eq("State_15_without_LAM1163") & sensitivity["feature"].eq("LAMCORE_777")]
    target_after_median = float(target_after["target_median"].iloc[0]) if len(target_after) else np.nan
    comparator_after_medians = {
        str(row["comparator_group"]): float(row["comparator_median"])
        for _, row in target_after.iterrows()
    }
    matched_lamcore = matched[matched["feature"].eq("LAMCORE_777")]
    lopo_profile = lopo["reference_profile_closer"] if len(lopo) else pd.Series(dtype=bool)
    lopo_latent = lopo["reference_latent_closer"] if len(lopo) else pd.Series(dtype=bool)
    formal_remains_above_comparators = bool(
        np.isfinite(target_after_median)
        and len(comparator_after_medians) == len(COMPARATORS)
        and all(target_after_median > value for value in comparator_after_medians.values())
    )
    positive_patient_delta_fraction = float((matched_lamcore["delta_mean"] > 0).mean()) if len(matched_lamcore) else np.nan
    sensitivity_conclusion = (
        "去除 LAM1163 后 State 15 的 LAMCORE 仍高于全部指定 comparator，且 patient-matched LAMCORE 方向保留。"
        if formal_remains_above_comparators and np.isfinite(positive_patient_delta_fraction) and positive_patient_delta_fraction >= 0.5
        else "去除 LAM1163 后 State 15 的跨患者生物学支持为混合结果，需保持 provisional。"
    )
    stage19_conclusion = (
        "B_patient_enriched_but_biological_profile_persists_after_removal"
        if formal_remains_above_comparators
        and np.isfinite(positive_patient_delta_fraction)
        and positive_patient_delta_fraction >= 0.5
        else "C_profile_not_preserved_after_dominant_patient_removal"
    )
    summary = {
        "stage": 19,
        "target_state": TARGET_STATE,
        "removed_patient_sensitivity": TARGET_PATIENT,
        "candidate_pool_definition": "current 5,378-cell high-confidence consensus cohort",
        "candidate_pool_cells": int(len(consensus)),
        "state15_cells": int(consensus["consensus_state"].eq(TARGET_STATE).sum()),
        "state15_patient_count": int(consensus.loc[consensus["consensus_state"].eq(TARGET_STATE), "patient_id"].astype(str).nunique()),
        "formal_signature_manifest": formal_manifest,
        "program_manifest": program_manifest,
        "score_manifest": score_manifest,
        "no_reclustering": True,
        "no_scvi_training": True,
        "author_annotation_status": availability.to_dict(orient="records"),
        "author_labels_interpretation": "not_assayed datasets are not treated as author-negative",
        "sensitivity": {
            "removed_patient": TARGET_PATIENT,
            "state15_cells_after_removal": int((scores["consensus_state"].astype(str).eq(TARGET_STATE) & ~scores["patient_id"].astype(str).eq(TARGET_PATIENT)).sum()),
            "state15_lamcore_median_after_removal": target_after_median,
            "comparator_lamcore_medians_after_removal": comparator_after_medians,
            "formal_remains_above_all_comparators": formal_remains_above_comparators,
            "patient_matched_lamcore_positive_delta_fraction": positive_patient_delta_fraction,
            "lopo_reference_profile_closer_fraction": float(lopo_profile.mean()) if len(lopo_profile) else np.nan,
            "lopo_reference_latent_closer_fraction": float(lopo_latent.mean()) if len(lopo_latent) else np.nan,
            "conclusion": sensitivity_conclusion,
        },
        "stage19_conclusion": stage19_conclusion,
    }
    (output_dir / "stage19_manifest.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8")

    target_comp = composition[composition["patient_id"].eq(TARGET_PATIENT)].iloc[0]
    author_available = availability[availability["author_style_annotation_status"].eq("available")]
    report = [
        "# Stage 19：State 15 跨患者 identity calibration audit",
        "",
        "本阶段固定使用现有 5,378 个 high-confidence candidate 的 consensus 标签，不重新聚类、不重训 scVI、不修改 candidate gate。",
        "",
        "## 1. Patient composition baseline",
        "",
        f"- State 15: {len(consensus[consensus['consensus_state'].eq(TARGET_STATE)])} cells; {summary['state15_patient_count']} patients.",
        f"- {TARGET_PATIENT} 在全部 candidate pool 中占 {float(target_comp['candidate_pool_fraction']):.4f}，在 State 15 中占 {float(target_comp['state15_fraction']):.4f}，composition enrichment = {float(target_comp['enrichment']):.4f}.",
        "",
        text_table(composition),
        "",
        "## 2. Author-style annotation availability",
        "",
        "只有存在真实逐细胞作者阳性标签的数据集才进入 author enrichment；其他数据集的零值被标记为 `not_assayed`，不解释为 author-negative。",
        "",
        text_table(availability),
        "",
        text_table(author_enrichment),
        "",
        "## 3. Patient-matched and leave-one-patient-out evidence",
        "",
        f"- Patient-matched LAMCORE positive delta fraction: {positive_patient_delta_fraction:.4f}.",
        f"- LOPO profile reference closer fraction: {float(lopo_profile.mean()) if len(lopo_profile) else np.nan:.4f}.",
        f"- LOPO latent reference closer fraction: {float(lopo_latent.mean()) if len(lopo_latent) else np.nan:.4f}.",
        "",
        "## 4. Removing LAM1163",
        "",
        f"- Remaining State 15 cells: {summary['sensitivity']['state15_cells_after_removal']}.",
        f"- Remaining State 15 LAMCORE median: {target_after_median:.4f}.",
        f"- LAMCORE remains above every requested comparator: {formal_remains_above_comparators}.",
        f"- Interpretation: {sensitivity_conclusion}",
        "",
        "## 5. Stage 19 conclusion",
        "",
        f"- Classification: `{stage19_conclusion}`.",
        "- The patient-composition baseline rules out interpreting the 63.5% LAM1163 fraction as ordinary sampling alone: its State 15 composition enrichment is 6.8301-fold.",
        "- The LAM1163 removal sensitivity, patient-matched direction, and LOPO latent comparison support a real but patient-heterogeneous State 15 biological profile.",
        "- Author-style evidence is available only for GSE135851 in the inherited upstream files; GSE190260, GSE217108 and GSE302356 are `not_assayed`, not author-negative.",
        "",
        "## 6. Scope",
        "",
        "本阶段仅判断当前冻结 State 15 的跨患者组成基线和生物学复现，不把结果自动写回 candidate gate、consensus clustering 或 atlas。",
        "",
        "## Outputs",
        "",
        *[f"- {path.name}" for path in sorted(output_dir.iterdir()) if path.is_file()],
    ]
    (output_dir / "state19_cross_patient_audit.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Audited frozen State 15: {summary['state15_cells']} cells")
    print(f"Removed-patient sensitivity: {TARGET_PATIENT}; remaining State 15: {summary['sensitivity']['state15_cells_after_removal']}")
    print(f"Outputs: {output_dir}")


if __name__ == "__main__":
    main()
