"""Release-aware LINCS drug x gene program analysis.

This is a downstream analysis of the already computed local WTCS/NCS results.
It reads Level 5 GCTX values for the combined plastic-or-hydrogel concordant
compound set and the genetic
perturbations relevant to lestaurtinib and QL-X-138.  Each LINCS release is
summarised independently before any cross-release comparison.

The primary evidence is direction and stability, not a product of two scores:
drug_effect, reversal/mimic/neutral status, and context consistency are kept
separately.  weighted_contribution is retained only as an auxiliary field.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np
import pandas as pd
import yaml
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import fisher_exact, spearmanr, trim_mean
from statsmodels.stats.multitest import multipletests

from common import (
    CANDIDATE_ANALYSIS_REPORTS,
    CANDIDATE_AUDIT,
    CANDIDATE_DRUG_TARGETS,
    CANDIDATE_PROGRAMS,
    CANDIDATE_VALIDATION,
    ROOT,
)


DATASETS = {
    "GSE92742": {
        "gctx": "data/processed/LINCS/gctx/GSE92742_Level5.gctx",
        "sig_info": "data/raw/LINCS/GSE92742/GSE92742_Broad_LINCS_sig_info.txt.gz",
        "gene_info": "data/raw/LINCS/GSE92742/GSE92742_Broad_LINCS_gene_info.txt.gz",
    },
    "GSE70138": {
        "gctx": "data/processed/LINCS/gctx/GSE70138_Level5.gctx",
        "sig_info": "data/raw/LINCS/GSE70138/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
        "gene_info": "data/raw/LINCS/GSE70138/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
    },
}

MIN_VALID_SIGNATURES = 3
MIN_DIRECTION_FRACTION = 0.60
CLUSTER_DISTANCE_THRESHOLD = 0.75
TARGET_DRUGS = ("lestaurtinib", "QL-X-138")
TARGET_PERT_TYPES = ("trt_sh", "trt_sh.cgs", "trt_sh.css", "trt_oe", "trt_lig", "trt_xpr")
OUTPUT_PREFIX = "tsc2_loss_plastic_or_hydrogel_replicated_concordant"


def norm(value: object) -> str:
    return re.sub(r"\s+", " ", str(value).strip().casefold())


def decode(values: Iterable[object]) -> list[str]:
    result = []
    for value in values:
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
        result.append(str(value))
    return result


def load_panel(signature_path: Path, contrast: str, n_genes: int = 150) -> pd.DataFrame:
    signatures = pd.read_csv(signature_path)
    group = signatures.loc[signatures["contrast"].eq(contrast)].copy()
    if group.empty:
        raise ValueError(f"{contrast} is absent from the GSE179044 signature file")
    group["gene"] = group["gene"].astype(str).str.upper().str.strip()
    group = group.drop_duplicates(["direction", "gene"])
    rows = []
    for direction, ascending in (("up", False), ("down", True)):
        selected = (
            group.loc[group["direction"].eq(direction)]
            .sort_values(["moderated_q", "signed_score"], ascending=[True, ascending])
            .head(n_genes)
        )
        rows.append(selected[["gene", "direction", "signed_score", "moderated_q"]])
    panel = pd.concat(rows, ignore_index=True).drop_duplicates("gene")
    panel = panel.rename(columns={"signed_score": "disease_weight", "direction": "disease_direction"})
    panel["disease_weight"] = panel["disease_weight"].astype(float)
    return panel


def load_compound_scope(path: Path) -> pd.DataFrame:
    compounds = pd.read_csv(path)
    required = {"pert_iname", "direction_pattern"}
    missing = required - set(compounds.columns)
    if missing:
        raise ValueError(f"compound list missing columns: {sorted(missing)}")
    result = compounds[["pert_iname", "direction_pattern", "cross_phase_status"]].copy()
    result["pert_iname"] = result["pert_iname"].astype(str)
    result["entity_id"] = "compound::" + result["pert_iname"].map(norm)
    result["scope"] = result["direction_pattern"].map(
        {"reversal_direction": "reversal_only", "mimic_direction": "mimic_only"}
    ).fillna("unclassified")
    return result


def load_target_axes(path: Path) -> dict[str, list[str]]:
    targets = pd.read_csv(path, low_memory=False)
    axes: dict[str, set[str]] = {drug: set() for drug in TARGET_DRUGS}
    for drug in TARGET_DRUGS:
        subset = targets.loc[targets["pert_iname"].map(norm).eq(norm(drug))]
        for value in subset.get("target_gene_symbol", pd.Series(dtype=str)).dropna():
            for gene in str(value).split(";"):
                gene = gene.strip().upper()
                if gene and gene != "NAN":
                    axes[drug].add(gene)
    return {drug: sorted(genes) for drug, genes in axes.items()}


def load_gene_mapping(dataset: str, panel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    cfg = DATASETS[dataset]
    gene_info = pd.read_csv(ROOT / cfg["gene_info"], sep="\t", compression="gzip", low_memory=False)
    bing = gene_info.loc[gene_info["pr_is_bing"].astype(int).eq(1)].copy()
    bing["gene"] = bing["pr_gene_symbol"].astype(str).str.upper().str.strip()
    symbol_to_id = {}
    for _, row in bing.iterrows():
        symbol_to_id.setdefault(row["gene"], str(row["pr_gene_id"]))
    with h5py.File(ROOT / cfg["gctx"], "r") as h5:
        gene_ids = decode(h5["0/META/ROW/id"][:])
        gene_lookup = {gene_id: i for i, gene_id in enumerate(gene_ids)}
    mapped = panel.copy()
    mapped["gene_id"] = mapped["gene"].map(symbol_to_id)
    mapped["matrix_col"] = mapped["gene_id"].map(gene_lookup)
    mapped["gene_available"] = mapped["matrix_col"].notna()
    mapped["matrix_col"] = mapped["matrix_col"].fillna(-1).astype(int)
    mapped = mapped.sort_values("gene").reset_index(drop=True)
    return mapped, bing, gene_ids, sorted(gene_lookup)


def load_sig_info(dataset: str, compound_scope: pd.DataFrame, target_axes: dict[str, list[str]]) -> tuple[pd.DataFrame, dict[str, int]]:
    cfg = DATASETS[dataset]
    sig = pd.read_csv(ROOT / cfg["sig_info"], sep="\t", compression="gzip", low_memory=False)
    with h5py.File(ROOT / cfg["gctx"], "r") as h5:
        sig_ids = decode(h5["0/META/COL/id"][:])
    sig_lookup = {sig_id: i for i, sig_id in enumerate(sig_ids)}
    sig["sig_id"] = sig["sig_id"].astype(str)
    sig["pert_iname"] = sig["pert_iname"].astype(str)
    sig["pert_type"] = sig["pert_type"].astype(str)
    sig["cell_id"] = sig["cell_id"].astype(str)
    sig["name_key"] = sig["pert_iname"].map(norm)
    compound_map = {norm(row.pert_iname): row.entity_id for row in compound_scope.itertuples()}
    compound_scope_map = {row.entity_id: row.scope for row in compound_scope.itertuples()}
    target_map = {gene.casefold(): gene for genes in target_axes.values() for gene in genes}
    is_compound = sig["pert_type"].eq("trt_cp") & sig["name_key"].isin(compound_map)
    is_target = sig["pert_type"].str.startswith("trt_") & sig["name_key"].isin(target_map)
    selected = sig.loc[is_compound | is_target].copy()
    selected["entity_id"] = selected["name_key"].map(compound_map)
    selected["entity_type"] = "compound"
    selected.loc[is_target.loc[selected.index], "entity_id"] = selected.loc[is_target.loc[selected.index], ["name_key", "pert_type"]].apply(
        lambda row: "genetic::" + target_map[row["name_key"]] + "::" + row["pert_type"], axis=1
    )
    selected.loc[is_target.loc[selected.index], "entity_type"] = "genetic"
    selected["scope"] = selected["entity_id"].map(compound_scope_map).fillna("target_validation")
    selected["gctx_row"] = selected["sig_id"].map(sig_lookup)
    selected = selected.loc[selected["gctx_row"].notna()].copy()
    selected["gctx_row"] = selected["gctx_row"].astype(int)
    selected["dose"] = _first_existing(selected, ["pert_idose", "pert_dose"])
    selected["time"] = _first_existing(selected, ["pert_itime", "pert_time"])
    selected["context_id"] = (
        selected["cell_id"].astype(str) + "|" + selected["dose"].astype(str) + "|" + selected["time"].astype(str)
    )
    selected = selected.sort_values("gctx_row").reset_index(drop=True)
    return selected, sig_lookup


def _first_existing(frame: pd.DataFrame, columns: list[str]) -> pd.Series:
    result = pd.Series("", index=frame.index, dtype=object)
    for column in columns:
        if column in frame.columns:
            values = frame[column].fillna("").astype(str)
            result = result.where(result.ne(""), values)
    return result


def summarize_release(
    dataset: str,
    mapped_panel: pd.DataFrame,
    selected: pd.DataFrame,
    entity_defs: pd.DataFrame,
    chunk_size: int = 256,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    cfg = DATASETS[dataset]
    mapped = mapped_panel.loc[mapped_panel["gene_available"]].copy()
    genes = mapped["gene"].tolist()
    weights = mapped["disease_weight"].to_numpy(float)
    matrix_cols = mapped["matrix_col"].to_numpy(int)
    n_entities = len(entity_defs)
    n_genes = len(genes)
    entity_index = {entity: i for i, entity in enumerate(entity_defs["entity_id"])}
    selected = selected.loc[selected["entity_id"].isin(entity_index)].copy()
    selected["entity_idx"] = selected["entity_id"].map(entity_index).astype(int)
    cell_levels = {value: i for i, value in enumerate(sorted(selected["cell_id"].unique()))}
    context_levels = {value: i for i, value in enumerate(sorted(selected["context_id"].unique()))}
    selected["cell_idx"] = selected["cell_id"].map(cell_levels).astype(int)
    selected["context_idx"] = selected["context_id"].map(context_levels).astype(int)

    values_acc = [[[] for _ in range(n_genes)] for _ in range(n_entities)]
    rev_counts = np.zeros((n_entities, n_genes), dtype=np.int64)
    mimic_counts = np.zeros((n_entities, n_genes), dtype=np.int64)
    neutral_counts = np.zeros((n_entities, n_genes), dtype=np.int64)
    cell_sum = np.zeros((n_entities, n_genes, len(cell_levels)), dtype=np.float64)
    cell_count = np.zeros((n_entities, n_genes, len(cell_levels)), dtype=np.int64)
    context_sum = np.zeros((n_entities, n_genes, len(context_levels)), dtype=np.float64)
    context_count = np.zeros((n_entities, n_genes, len(context_levels)), dtype=np.int64)

    with h5py.File(ROOT / cfg["gctx"], "r") as h5:
        matrix = h5["0/DATA/0/matrix"]
        for start in range(0, len(selected), chunk_size):
            chunk = selected.iloc[start : start + chunk_size]
            rows = chunk["gctx_row"].to_numpy(int)
            values = np.asarray(matrix[rows, :][:, matrix_cols], dtype=np.float32)
            for i, (_, meta) in enumerate(chunk.iterrows()):
                entity = int(meta["entity_idx"])
                row = values[i].astype(float)
                sign = row * weights
                rev = sign < 0
                mimic = sign > 0
                neutral = ~(rev | mimic)
                rev_counts[entity] += rev.astype(np.int64)
                mimic_counts[entity] += mimic.astype(np.int64)
                neutral_counts[entity] += neutral.astype(np.int64)
                for gene_idx, value in enumerate(row):
                    if np.isfinite(value):
                        values_acc[entity][gene_idx].append(float(value))
                cell = int(meta["cell_idx"])
                context = int(meta["context_idx"])
                cell_sum[entity, :, cell] += np.nan_to_num(row, nan=0.0)
                cell_count[entity, :, cell] += np.isfinite(row).astype(np.int64)
                context_sum[entity, :, context] += np.nan_to_num(row, nan=0.0)
                context_count[entity, :, context] += np.isfinite(row).astype(np.int64)
            if start == 0 or start % (chunk_size * 20) == 0:
                print(f"[{dataset}] gene responses {start + len(chunk):,}/{len(selected):,}", flush=True)

    entity_meta = entity_defs.set_index("entity_id").to_dict("index")
    rows = []
    for entity_id, entity_idx in entity_index.items():
        meta = entity_meta[entity_id]
        for gene_idx, gene in enumerate(genes):
            values = np.asarray(values_acc[entity_idx][gene_idx], dtype=float)
            values = values[np.isfinite(values)]
            n_valid = len(values)
            weight = float(weights[gene_idx])
            rev_fraction = float(rev_counts[entity_idx, gene_idx] / n_valid) if n_valid else np.nan
            mimic_fraction = float(mimic_counts[entity_idx, gene_idx] / n_valid) if n_valid else np.nan
            neutral_fraction = float(neutral_counts[entity_idx, gene_idx] / n_valid) if n_valid else np.nan
            if n_valid < MIN_VALID_SIGNATURES:
                direction = "not_available"
                status = "not_available"
            elif rev_fraction >= MIN_DIRECTION_FRACTION:
                direction = "reversal"
                status = "available"
            elif mimic_fraction >= MIN_DIRECTION_FRACTION:
                direction = "mimic"
                status = "available"
            else:
                direction = "neutral"
                status = "available"
            cell_means = np.divide(
                cell_sum[entity_idx, gene_idx],
                cell_count[entity_idx, gene_idx],
                out=np.full(len(cell_levels), np.nan),
                where=cell_count[entity_idx, gene_idx] > 0,
            )
            context_means = np.divide(
                context_sum[entity_idx, gene_idx],
                context_count[entity_idx, gene_idx],
                out=np.full(len(context_levels), np.nan),
                where=context_count[entity_idx, gene_idx] > 0,
            )
            cell_sign = cell_means * weight
            context_sign = context_means * weight
            cell_n = int(np.isfinite(cell_means).sum())
            context_n = int(np.isfinite(context_means).sum())
            cell_rev = float(np.mean(cell_sign[np.isfinite(cell_sign)] < 0)) if cell_n else np.nan
            cell_mimic = float(np.mean(cell_sign[np.isfinite(cell_sign)] > 0)) if cell_n else np.nan
            context_rev = float(np.mean(context_sign[np.isfinite(context_sign)] < 0)) if context_n else np.nan
            context_mimic = float(np.mean(context_sign[np.isfinite(context_sign)] > 0)) if context_n else np.nan
            rows.append(
                {
                    "dataset": dataset,
                    "entity_id": entity_id,
                    "entity_type": meta["entity_type"],
                    "pert_iname": meta["pert_iname"],
                    "scope": meta["scope"],
                    "gene": gene,
                    "disease_weight": weight,
                    "drug_effect_median": float(np.median(values)) if n_valid else np.nan,
                    "drug_effect_trimmed_mean": float(trim_mean(values, 0.2)) if n_valid >= 5 else (float(values.mean()) if n_valid else np.nan),
                    "drug_effect_q25": float(np.quantile(values, 0.25)) if n_valid else np.nan,
                    "drug_effect_q75": float(np.quantile(values, 0.75)) if n_valid else np.nan,
                    "weighted_contribution_median": float(np.median(-weight * values)) if n_valid else np.nan,
                    "n_valid_signatures": n_valid,
                    "n_cells": cell_n,
                    "n_contexts": context_n,
                    "reversal_fraction": rev_fraction,
                    "mimic_fraction": mimic_fraction,
                    "neutral_fraction": neutral_fraction,
                    "cell_line_reversal_fraction": cell_rev,
                    "cell_line_mimic_fraction": cell_mimic,
                    "cell_line_consistency": max(cell_rev, cell_mimic) if cell_n else np.nan,
                    "context_reversal_fraction": context_rev,
                    "context_mimic_fraction": context_mimic,
                    "context_consistency": max(context_rev, context_mimic) if context_n else np.nan,
                    "direction_status": direction,
                    "data_status": status,
                    "gene_available": True,
                }
            )
    # Add genes absent from this release explicitly as not_available.
    mapped_symbols = set(genes)
    unavailable = mapped_panel.loc[~mapped_panel["gene"].isin(mapped_symbols)]
    for entity_id, entity_idx in entity_index.items():
        meta = entity_meta[entity_id]
        for _, panel_row in unavailable.iterrows():
            rows.append(
                {
                    "dataset": dataset,
                    "entity_id": entity_id,
                    "entity_type": meta["entity_type"],
                    "pert_iname": meta["pert_iname"],
                    "scope": meta["scope"],
                    "gene": panel_row["gene"],
                    "disease_weight": panel_row["disease_weight"],
                    "drug_effect_median": np.nan,
                    "drug_effect_trimmed_mean": np.nan,
                    "drug_effect_q25": np.nan,
                    "drug_effect_q75": np.nan,
                    "weighted_contribution_median": np.nan,
                    "n_valid_signatures": 0,
                    "n_cells": 0,
                    "n_contexts": 0,
                    "reversal_fraction": np.nan,
                    "mimic_fraction": np.nan,
                    "neutral_fraction": np.nan,
                    "cell_line_reversal_fraction": np.nan,
                    "cell_line_mimic_fraction": np.nan,
                    "cell_line_consistency": np.nan,
                    "context_reversal_fraction": np.nan,
                    "context_mimic_fraction": np.nan,
                    "context_consistency": np.nan,
                    "direction_status": "not_available",
                    "data_status": "not_available",
                    "gene_available": False,
                }
            )
    summary = pd.DataFrame(rows)
    summary["direction_stability"] = summary["reversal_fraction"] - summary["mimic_fraction"]
    selected_audit = selected.groupby(["entity_id", "entity_type", "scope"], as_index=False).agg(
        n_selected_signatures=("sig_id", "size"),
        n_selected_cells=("cell_id", "nunique"),
        perturbation_types=("pert_type", lambda x: "|".join(sorted(set(x.astype(str))))),
    )
    summary = summary.merge(
        selected_audit[["entity_id", "perturbation_types"]],
        on="entity_id",
        how="left",
        validate="many_to_one",
    )
    return summary, selected_audit


def add_not_available_entities(summary: pd.DataFrame, entity_defs: pd.DataFrame, panel: pd.DataFrame, dataset: str) -> pd.DataFrame:
    present = set(summary["entity_id"])
    rows = []
    for entity in entity_defs.itertuples():
        if entity.entity_id in present:
            continue
        for panel_row in panel.itertuples():
            rows.append(
                {
                    "dataset": dataset,
                    "entity_id": entity.entity_id,
                    "entity_type": entity.entity_type,
                    "pert_iname": entity.pert_iname,
                    "scope": entity.scope,
                    "gene": panel_row.gene,
                    "disease_weight": panel_row.disease_weight,
                    "direction_status": "not_available",
                    "data_status": "not_available",
                    "gene_available": False,
                    "n_valid_signatures": 0,
                    "n_cells": 0,
                    "n_contexts": 0,
                }
            )
    return pd.concat([summary, pd.DataFrame(rows)], ignore_index=True) if rows else summary


def make_matrix(summary: pd.DataFrame, value: str, entity_type: str, scope: str | None = None) -> pd.DataFrame:
    frame = summary.loc[summary["entity_type"].eq(entity_type)].copy()
    if scope is not None:
        frame = frame.loc[frame["scope"].eq(scope)]
    return frame.pivot_table(index="entity_id", columns="gene", values=value, aggfunc="first")


def hierarchical_clusters(matrix: pd.DataFrame, object_type: str, dataset: str, scope: str) -> pd.DataFrame:
    if matrix.shape[0] < 2:
        return pd.DataFrame(columns=["dataset", "scope", "object_type", "object_id", "cluster_id", "n_features"])
    values = matrix.astype(float).fillna(0.0).to_numpy()
    if np.allclose(values, values[0]):
        labels = np.ones(matrix.shape[0], dtype=int)
    else:
        distances = pdist(values, metric="cosine")
        distances = np.nan_to_num(distances, nan=1.0, posinf=1.0, neginf=1.0)
        tree = linkage(distances, method="average")
        labels = fcluster(tree, t=CLUSTER_DISTANCE_THRESHOLD, criterion="distance")
    return pd.DataFrame(
        {
            "dataset": dataset,
            "scope": scope,
            "object_type": object_type,
            "object_id": matrix.index.astype(str),
            "cluster_id": labels.astype(int),
            "n_features": matrix.shape[1],
        }
    )


def cluster_release(summary: pd.DataFrame, dataset: str, scope: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    direction = make_matrix(summary, "direction_stability", "compound", scope)
    drug_clusters = hierarchical_clusters(direction, "drug", dataset, scope)
    gene_clusters = hierarchical_clusters(direction.T, "gene", dataset, scope)
    return drug_clusters, gene_clusters


def module_matches(gene_clusters: pd.DataFrame, common_genes: set[str]) -> pd.DataFrame:
    if gene_clusters.empty:
        return pd.DataFrame()
    grouped = {
        (dataset, scope, int(cluster)): set(group.loc[group["object_type"].eq("gene"), "object_id"]) & common_genes
        for (dataset, scope, cluster), group in gene_clusters.groupby(["dataset", "scope", "cluster_id"])
    }
    rows = []
    left_keys = [key for key in grouped if key[0] == "GSE92742"]
    right_keys = [key for key in grouped if key[0] == "GSE70138"]
    for left in left_keys:
        for right in right_keys:
            a, b = grouped[left], grouped[right]
            if left[1] != right[1]:
                continue
            union = a | b
            intersection = a & b
            jaccard = len(intersection) / len(union) if union else 0.0
            overlap = len(intersection) / min(len(a), len(b)) if a and b else 0.0
            rows.append(
                {
                    "scope": left[1],
                    "module_GSE92742": left[2],
                    "module_GSE70138": right[2],
                    "n_genes_GSE92742_common_panel": len(a),
                    "n_genes_GSE70138_common_panel": len(b),
                    "n_common_genes": len(intersection),
                    "jaccard": jaccard,
                    "overlap_coefficient": overlap,
                    "same_module_first_pass": bool(jaccard >= 0.5 and len(intersection) >= 3),
                    "common_genes": ";".join(sorted(intersection)),
                }
            )
    return pd.DataFrame(rows)


def posthoc_module_overlap(summary: pd.DataFrame, gene_clusters: pd.DataFrame, config_path: Path) -> pd.DataFrame:
    config = yaml.safe_load((ROOT / config_path).read_text())
    module_sets = {name: {str(g).upper() for g in genes} for name, genes in config.get("module_sets", {}).items()}
    rows = []
    for (dataset, scope, cluster), group in gene_clusters.groupby(["dataset", "scope", "cluster_id"]):
        cluster_genes = set(group["object_id"].astype(str))
        background = set(
            summary.loc[
                summary["dataset"].eq(dataset)
                & summary["scope"].eq(scope)
                & summary["entity_type"].eq("compound")
                & summary["data_status"].eq("available"),
                "gene",
            ]
        )
        # The background is the actual clusterable gene universe for this
        # release/scope, never the whole human genome and never only the
        # current cluster.
        if not background:
            background = cluster_genes
        for module_name, module_genes in module_sets.items():
            a = len(cluster_genes & module_genes)
            b = len(cluster_genes - module_genes)
            c = len((background & module_genes) - cluster_genes)
            d = len(background - cluster_genes - module_genes)
            _, p_value = fisher_exact([[a, b], [c, d]], alternative="greater")
            rows.append(
                {
                    "dataset": dataset,
                    "scope": scope,
                    "cluster_id": int(cluster),
                    "cluster_size": len(cluster_genes),
                    "module": module_name,
                    "overlap_genes": a,
                    "background_size": len(background),
                    "p_value": p_value,
                    "background_definition": "actual_clusterable_genes",
                }
            )
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = multipletests(result["p_value"].fillna(1.0), method="fdr_bh")[1]
        result = result.sort_values(["fdr", "p_value"])
    return result


def target_axis_analysis(summary: pd.DataFrame, axes: dict[str, list[str]]) -> pd.DataFrame:
    rows = []
    for dataset in sorted(summary["dataset"].dropna().unique()):
        dataset_summary = summary.loc[summary["dataset"].eq(dataset)]
        for drug, targets in axes.items():
            drug_id = "compound::" + norm(drug)
            drug_frame = dataset_summary.loc[
                dataset_summary["entity_id"].eq(drug_id) & dataset_summary["data_status"].eq("available")
            ].drop_duplicates("gene").set_index("gene")
            drug_reversal = set(drug_frame.index[drug_frame["direction_status"].eq("reversal")])
            for target in targets:
                target_prefix = "genetic::" + target + "::"
                target_ids = sorted(
                    entity_id
                    for entity_id in dataset_summary["entity_id"].dropna().unique()
                    if str(entity_id).startswith(target_prefix)
                )
                if not target_ids:
                    target_ids = [target_prefix + "not_available"]
                for target_id in target_ids:
                    target_perturbation_type = target_id.split("::", 2)[-1]
                    target_frame = dataset_summary.loc[
                        dataset_summary["entity_id"].eq(target_id) & dataset_summary["data_status"].eq("available")
                    ].drop_duplicates("gene").set_index("gene")
                    if target_frame.empty or drug_frame.empty:
                        rows.append(
                            {
                                "dataset": dataset,
                                "drug": drug,
                                "target_gene": target,
                                "target_perturbation_type": target_perturbation_type,
                                "status": "not_available",
                                "n_shared_stable_reversal": 0,
                            }
                        )
                        continue
                    target_reversal = set(target_frame.index[target_frame["direction_status"].eq("reversal")])
                    shared = drug_reversal & target_reversal
                    union = drug_reversal | target_reversal
                    common = drug_frame.index.intersection(target_frame.index)
                    d_effect = drug_frame.loc[common, "drug_effect_median"].to_numpy(float)
                    t_effect = target_frame.loc[common, "drug_effect_median"].to_numpy(float)
                    valid = np.isfinite(d_effect) & np.isfinite(t_effect)
                    d_effect, t_effect = d_effect[valid], t_effect[valid]
                    direction_concordance = float(np.mean(np.sign(d_effect) == np.sign(t_effect))) if len(d_effect) else np.nan
                    rho = float(spearmanr(d_effect, t_effect).statistic) if len(d_effect) >= 3 else np.nan
                    overlap_jaccard = len(shared) / len(union) if union else 0.0
                    if len(shared) >= 5 and direction_concordance >= 0.6:
                        status = "supportive"
                    elif len(shared) >= 3 or (np.isfinite(direction_concordance) and direction_concordance >= 0.55):
                        status = "weak"
                    elif np.isfinite(direction_concordance) and direction_concordance < 0.4:
                        status = "discordant"
                    else:
                        status = "weak"
                    rows.append(
                        {
                            "dataset": dataset,
                            "drug": drug,
                            "target_gene": target,
                            "target_perturbation_type": target_perturbation_type,
                            "status": status,
                            "n_shared_stable_reversal": len(shared),
                            "drug_stable_reversal_genes": len(drug_reversal),
                            "target_stable_reversal_genes": len(target_reversal),
                            "stable_reversal_jaccard": overlap_jaccard,
                            "n_common_genes_compared": len(d_effect),
                            "direction_concordance": direction_concordance,
                            "spearman_sensitivity": rho,
                            "target_n_valid_signatures": int(target_frame["n_valid_signatures"].max()),
                            "target_perturbation_types": target_perturbation_type,
                        }
                    )
    return pd.DataFrame(rows)


def cross_release_gene_comparison(summary: pd.DataFrame, common_genes: set[str]) -> pd.DataFrame:
    compounds = summary.loc[summary["entity_type"].eq("compound") & summary["gene"].isin(common_genes)].copy()
    left = compounds.loc[compounds.dataset.eq("GSE92742")]
    right = compounds.loc[compounds.dataset.eq("GSE70138")]
    keys = ["entity_id", "gene"]
    merged = left.merge(right, on=keys, how="outer", suffixes=("_GSE92742", "_GSE70138"))
    statuses = []
    for row in merged.itertuples():
        a, b = row.direction_status_GSE92742, row.direction_status_GSE70138
        if a == "not_available" or b == "not_available" or pd.isna(a) or pd.isna(b):
            status = "not_comparable_not_available"
        elif a == b and a in {"reversal", "mimic"}:
            status = "concordant_direction"
        elif {a, b} == {"reversal", "mimic"}:
            status = "discordant_direction"
        else:
            status = "available_but_neutral_or_weak"
        statuses.append(status)
    merged["cross_release_gene_status"] = statuses
    merged["common_gene_panel"] = True
    return merged


def write_report(
    out_dir: Path,
    summary: pd.DataFrame,
    module_matches_df: pd.DataFrame,
    target_df: pd.DataFrame,
    common_genes: set[str],
    contrast: str,
    output_prefix: str,
) -> None:
    compounds = summary.loc[summary["entity_type"].eq("compound")]
    counts = compounds.groupby(["scope", "entity_id"], as_index=False).agg(
        n_reversal_genes=("direction_status", lambda x: int((x == "reversal").sum())),
        n_mimic_genes=("direction_status", lambda x: int((x == "mimic").sum())),
        n_neutral_genes=("direction_status", lambda x: int((x == "neutral").sum())),
        n_not_available_genes=("direction_status", lambda x: int((x == "not_available").sum())),
    )
    lines = [
        "# LINCS drug × gene program analysis",
        "",
        "本报告基于 GSE92742 与 GSE70138 的 Level 5 GCTX，两个 release 先独立汇总，再进行共同基因比较。",
        "方向判断与 drug effect 分开保存；weighted contribution 仅为辅助指标。",
        "",
        f"- primary disease contrast: {contrast} top150 up + top150 down",
        f"- common analyzable genes used for cross-release comparison: {len(common_genes)}",
        f"- compound rows: {len(compounds):,}",
        f"- cross-release candidate module matches: {int(module_matches_df['same_module_first_pass'].sum()) if not module_matches_df.empty else 0}",
        "- GO/Reactome/MSigDB: not run in this local-only pass because no local GMT was present; post-hoc predefined-module overlap was calculated separately.",
        "",
        "## Per-drug stable gene counts",
        "",
        counts.sort_values(["scope", "n_reversal_genes"], ascending=[True, False]).to_string(index=False),
        "",
        "## Lestaurtinib / QL-X-138 target-axis evidence",
        "",
        target_df.to_string(index=False) if not target_df.empty else "No target perturbation was available.",
        "",
        "Interpretation remains hypothesis-generating; target causality and generic cytotoxicity require external/selective perturbation validation.",
    ]
    (out_dir / f"{output_prefix}_LINCS_gene_program_analysis.md").write_text("\n".join(lines))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signature", default="results/signatures/GSE179044_cmap_query_signatures.csv")
    parser.add_argument("--compounds", default="results/candidates/tsc2_loss_plastic_or_hydrogel_replicated_concordant_compounds_unique.csv")
    parser.add_argument("--targets", default="data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_or_hydrogel_replicated_concordant_compound_targets.csv")
    parser.add_argument("--contrast", choices=("tsc2_loss_plastic", "tsc2_loss_hydrogel"), default="tsc2_loss_plastic")
    parser.add_argument("--output-prefix", default=OUTPUT_PREFIX)
    parser.add_argument("--config", default="config/analysis.yaml")
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        assert MIN_VALID_SIGNATURES == 3
        assert MIN_DIRECTION_FRACTION == 0.60
        assert norm("QL-X-138") == "ql-x-138"
        print("self_test: PASS")
        return

    out_dir = CANDIDATE_PROGRAMS
    out_dir.mkdir(parents=True, exist_ok=True)
    CANDIDATE_VALIDATION.mkdir(parents=True, exist_ok=True)
    CANDIDATE_AUDIT.mkdir(parents=True, exist_ok=True)
    CANDIDATE_ANALYSIS_REPORTS.mkdir(parents=True, exist_ok=True)
    panel = load_panel(ROOT / args.signature, args.contrast)
    compound_scope = load_compound_scope(ROOT / args.compounds)
    target_axes = load_target_axes(ROOT / args.targets)
    target_entities = sorted({gene for genes in target_axes.values() for gene in genes})
    entity_rows = []
    for row in compound_scope.itertuples():
        entity_rows.append({"entity_id": row.entity_id, "entity_type": "compound", "pert_iname": row.pert_iname, "scope": row.scope})
    for gene in target_entities:
        for pert_type in TARGET_PERT_TYPES:
            entity_rows.append(
                {
                    "entity_id": "genetic::" + gene + "::" + pert_type,
                    "entity_type": "genetic",
                    "pert_iname": gene,
                    "scope": "target_validation",
                }
            )
    entity_defs = pd.DataFrame(entity_rows)

    all_summaries = []
    audits = []
    mappings = []
    release_mapped: dict[str, set[str]] = {}
    for dataset in DATASETS:
        mapped, _, _, _ = load_gene_mapping(dataset, panel)
        mappings.append(mapped.assign(dataset=dataset))
        release_mapped[dataset] = set(mapped.loc[mapped["gene_available"], "gene"])
        selected, _ = load_sig_info(dataset, compound_scope, target_axes)
        summary, audit = summarize_release(dataset, mapped, selected, entity_defs, args.chunk_size)
        summary = add_not_available_entities(summary, entity_defs, panel, dataset)
        all_summaries.append(summary)
        audits.append(audit.assign(dataset=dataset))
        print(f"[{dataset}] selected signatures={len(selected):,}; summary rows={len(summary):,}", flush=True)

    summary = pd.concat(all_summaries, ignore_index=True)
    common_genes = release_mapped["GSE92742"] & release_mapped["GSE70138"]
    gene_comparison = cross_release_gene_comparison(summary, common_genes)

    clusters = []
    for dataset in DATASETS:
        release_summary = summary.loc[summary.dataset.eq(dataset)]
        for scope in ("reversal_only", "mimic_only"):
            drug_cluster, gene_cluster = cluster_release(release_summary, dataset, scope)
            clusters.extend([drug_cluster, gene_cluster])
    cluster_table = pd.concat([x for x in clusters if not x.empty], ignore_index=True) if any(not x.empty for x in clusters) else pd.DataFrame()
    gene_clusters = cluster_table.loc[cluster_table["object_type"].eq("gene")] if not cluster_table.empty else pd.DataFrame()
    matches = module_matches(gene_clusters, common_genes)
    posthoc = posthoc_module_overlap(summary, gene_clusters, ROOT / args.config)
    target_analysis = target_axis_analysis(summary, target_axes)

    summary.to_csv(out_dir / f"{args.output_prefix}_LINCS_drug_gene_response_summary.csv.gz", index=False, compression="gzip")
    for dataset in DATASETS:
        for scope in ("reversal_only", "mimic_only"):
            matrix = make_matrix(summary.loc[summary.dataset.eq(dataset)], "direction_stability", "compound", scope)
            matrix.to_csv(out_dir / f"{args.output_prefix}_{dataset}_{scope}_drug_gene_direction_stability.csv")
    gene_comparison.to_csv(out_dir / f"{args.output_prefix}_LINCS_cross_release_common_gene_comparison.csv.gz", index=False, compression="gzip")
    cluster_table.to_csv(out_dir / f"{args.output_prefix}_LINCS_drug_gene_clusters.csv", index=False)
    matches.to_csv(out_dir / f"{args.output_prefix}_LINCS_cross_release_module_matches.csv", index=False)
    posthoc.to_csv(out_dir / f"{args.output_prefix}_LINCS_posthoc_module_overlap.csv", index=False)
    target_analysis.to_csv(CANDIDATE_VALIDATION / f"{args.output_prefix}_LINCS_lestaurtinib_qlx138_target_axis_analysis.csv", index=False)
    pd.concat(mappings, ignore_index=True).to_csv(CANDIDATE_AUDIT / f"{args.output_prefix}_LINCS_gene_panel_mapping_audit.csv", index=False)
    pd.concat(audits, ignore_index=True).to_csv(CANDIDATE_AUDIT / f"{args.output_prefix}_LINCS_gene_analysis_signature_audit.csv", index=False)
    manifest = {
        "output_prefix": args.output_prefix,
        "input_compounds": args.compounds,
        "input_targets": args.targets,
        "primary_contrast": args.contrast,
        "primary_panel": f"{args.contrast} top150 up + top150 down",
        "n_primary_panel_genes": int(len(panel)),
        "n_common_genes_GSE92742_GSE70138": int(len(common_genes)),
        "n_unique_compounds": int(len(compound_scope)),
        "n_reversal_only_compounds": int((compound_scope.scope == "reversal_only").sum()),
        "n_mimic_only_compounds": int((compound_scope.scope == "mimic_only").sum()),
        "min_valid_signatures": MIN_VALID_SIGNATURES,
        "min_direction_fraction": MIN_DIRECTION_FRACTION,
        "neutral_definition": "available data but neither reversal nor mimic fraction reaches stability threshold",
        "not_available_definition": "insufficient valid signatures or gene not available in the release",
        "cross_release_gene_panel_definition": "intersection of release-specific analyzable genes",
        "module_match_rule": "Jaccard >= 0.5 and at least 3 shared genes",
        "enrichment_background": "actual clusterable genes per release; common analyzable genes for cross-release comparison",
        "go_reactome_msigdb": "not_run_no_local_gmt",
    }
    (ROOT / "manifests" / f"{args.output_prefix}_LINCS_gene_program_analysis_manifest.json").write_text(json.dumps(manifest, indent=2))
    write_report(CANDIDATE_ANALYSIS_REPORTS, summary, matches, target_analysis, common_genes, args.contrast, args.output_prefix)
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
