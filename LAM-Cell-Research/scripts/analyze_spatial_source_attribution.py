#!/usr/bin/env python3
"""Estimate protease source attribution and proteolytic balance.

The script intentionally does not assign a protease gene to a cell type in
advance. Identity markers are used to define broad source states from actual
single-cell expression; protease genes are then measured independently in
those states. Spatial output is modality-specific and uses a standardized
protease-minus-antiprotease balance, never a score ratio.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from sklearn.neighbors import NearestNeighbors

from analyze_spatial_niche import decode, load_coordinates, read_10x_target_matrix

ROOT = Path(__file__).resolve().parents[1]

PROTEASE_GENES = ["CTSK", "MMP2", "MMP8", "MMP9", "MMP12", "ELANE", "PRTN3", "CTSS", "CTSB"]
ANTIPROTEASE_GENES = ["TIMP1", "TIMP2", "SERPINA1", "SERPINE1", "SLPI", "PI3"]
SOURCE_MARKERS = {
    "LAMCORE_like": ["PMEL", "MLANA", "MITF", "ACTA2", "ESR1", "FIGF", "VEGFD", "CTSK"],
    "LAF_fibroblast": ["PDGFRA", "COL1A1", "COL1A2", "DCN", "LUM", "CFD", "C7", "CXCL14", "COL3A1"],
    "immune": ["PTPRC", "LST1", "AIF1", "TYROBP", "FCER1G", "CD68", "CTSS", "CD3D", "CD3E"],
    "lymphatic_endothelial": ["PDPN", "LYVE1", "FLT4", "CCL21", "PROX1", "KDR", "PECAM1", "EMCN"],
    "vascular_endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2", "CA4"],
}


def zscore(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    mean = np.nanmean(values)
    sd = np.nanstd(values)
    return (values - mean) / sd if np.isfinite(sd) and sd > 0 else np.zeros_like(values, dtype=float)


def gene_symbols(adata: ad.AnnData) -> np.ndarray:
    for col in ["gene_symbol_upper", "gene_symbol"]:
        if col in adata.var:
            return adata.var[col].astype(str).str.upper().to_numpy()
    return adata.var_names.astype(str).str.upper().to_numpy()


def source_labels_from_core(obs: pd.DataFrame, norm: np.ndarray, selected_symbols: np.ndarray) -> pd.Series:
    labels = pd.Series("unclassified", index=obs.index, dtype="object")
    if "lamcore_candidate_author_style" in obs:
        labels.loc[obs["lamcore_candidate_author_style"].astype(bool)] = "LAMCORE_like"
    # Source labels are estimated from the actual expression of non-protease
    # identity markers, not inherited from broad precomputed microenvironment
    # scores whose baselines may differ substantially between sources.
    source_scores = {}
    for source, markers in SOURCE_MARKERS.items():
        rows = [i for i, gene in enumerate(selected_symbols) if gene in markers and gene not in PROTEASE_GENES + ANTIPROTEASE_GENES]
        if rows:
            source_scores[source] = np.nanmean(norm[:, rows], axis=1)
    if source_scores:
        score_matrix = np.column_stack([zscore(values) for values in source_scores.values()])
        names = list(source_scores)
        best = np.nanargmax(np.nan_to_num(score_matrix, nan=-np.inf), axis=1)
        best_values = np.nanmax(np.nan_to_num(score_matrix, nan=-np.inf), axis=1)
        for i, value in enumerate(best_values):
            if labels.iloc[i] == "LAMCORE_like":
                continue
            if np.isfinite(value) and value >= 0.5:
                labels.iloc[i] = names[best[i]]
    return labels


def load_core_source_reference(path: Path, output_dir: Path) -> tuple[pd.DataFrame, dict]:
    adata = ad.read_h5ad(path, backed="r")
    symbols = gene_symbols(adata)
    target_genes = sorted(set(PROTEASE_GENES + ANTIPROTEASE_GENES + sum(SOURCE_MARKERS.values(), [])))
    indices = np.flatnonzero(np.isin(symbols, target_genes))
    if "counts" in adata.layers:
        matrix = adata[:, indices].layers["counts"]
    else:
        matrix = adata[:, indices].X
    raw = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    raw = np.asarray(raw, dtype=float)
    totals = pd.to_numeric(adata.obs.get("total_counts", pd.Series(raw.sum(axis=1), index=adata.obs.index)), errors="coerce").to_numpy(float)
    norm = np.log1p(raw / np.maximum(totals[:, None], 1.0) * 1e4)
    selected_symbols = symbols[indices]
    labels = source_labels_from_core(adata.obs, norm, selected_symbols)
    obs = adata.obs.copy()
    obs["source_state"] = labels.to_numpy()
    obs_index = pd.Index(obs.index)
    cell_ids = obs["cell_id"].astype(str).to_numpy() if "cell_id" in obs else obs.index.astype(str).to_numpy()
    source_rows = []
    donor_rows = []
    for gene in PROTEASE_GENES + ANTIPROTEASE_GENES:
        gene_idx = np.flatnonzero(selected_symbols == gene)
        if len(gene_idx) == 0:
            continue
        expression = raw[:, gene_idx].sum(axis=1)
        for source, group in obs.groupby("source_state", observed=True):
            ix = obs_index.get_indexer(group.index)
            ix = ix[ix >= 0]
            if len(ix) == 0:
                continue
            values = expression[ix]
            source_rows.append({
                "gene": gene,
                "gene_class": "protease" if gene in PROTEASE_GENES else "antiprotease",
                "source_state": source,
                "n_cells": int(len(ix)),
                "mean_counts_per_cell": float(values.mean()),
                "pct_positive": float((values > 0).mean()),
                "transcript_fraction": float(values.sum() / expression.sum()) if expression.sum() > 0 else np.nan,
            })
            for donor_id, donor_group in group.groupby("donor_id", observed=True):
                dix = obs_index.get_indexer(donor_group.index)
                dix = dix[dix >= 0]
                dvalues = expression[dix]
                donor_total = expression[obs_index.get_indexer(obs.loc[obs["donor_id"].eq(donor_id)].index)]
                donor_rows.append({
                    "donor_id": str(donor_id),
                    "gene": gene,
                    "gene_class": "protease" if gene in PROTEASE_GENES else "antiprotease",
                    "source_state": source,
                    "n_cells": int(len(dix)),
                    "mean_counts_per_cell": float(dvalues.mean()),
                    "pct_positive": float((dvalues > 0).mean()),
                    "transcript_fraction_within_donor": float(dvalues.sum() / donor_total.sum()) if donor_total.sum() > 0 else np.nan,
                })
    source_df = pd.DataFrame(source_rows)
    donor_df = pd.DataFrame(donor_rows)
    source_df.to_csv(output_dir / "single_cell_protease_source_attribution.csv", index=False)
    donor_df.to_csv(output_dir / "single_cell_protease_source_attribution_by_donor.csv", index=False)
    # Store source marker expression summaries for spatial projection. Protease
    # and antiprotease genes are excluded from the projection features.
    projection_genes = sorted(set(sum(SOURCE_MARKERS.values(), [])) - set(PROTEASE_GENES + ANTIPROTEASE_GENES))
    profiles = []
    for source in SOURCE_MARKERS:
        marker_idx = [i for i, g in enumerate(selected_symbols) if g in projection_genes and g in SOURCE_MARKERS[source]]
        if not marker_idx:
            continue
        subset = norm[:, marker_idx]
        for gene_i, gene in zip(marker_idx, selected_symbols[marker_idx]):
            for source_state in [source]:
                profiles.append({"source_state": source_state, "gene": gene, "mean_log1p_cpm": float(subset[:, list(selected_symbols[marker_idx]).index(gene)].mean())})
    pd.DataFrame(profiles).to_csv(output_dir / "source_reference_marker_profiles.csv", index=False)
    manifest = {
        "reference": str(path.relative_to(ROOT)),
        "n_cells": int(adata.n_obs),
        "source_label_rule": "LAMCORE from existing candidate label; other broad source states from z-scored non-protease identity-marker expression in the actual single-cell matrix; protease genes were not used to assign labels.",
        "protease_genes": PROTEASE_GENES,
        "antiprotease_genes": ANTIPROTEASE_GENES,
        "external_datasets": "External h5ad objects are retained for later projection; current source labels are not silently inferred where cell-state metadata are absent.",
    }
    (output_dir / "source_reference_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    return obs[["source_state", "donor_id"]].copy(), manifest


def spatial_scores(modality: str, cfg: dict, out: Path, source_reference: pd.DataFrame, max_units: int, seed: int) -> dict:
    all_genes = sorted(set(PROTEASE_GENES + ANTIPROTEASE_GENES + sum(SOURCE_MARKERS.values(), [])))
    feature_names, barcodes, target, totals, detected = read_10x_target_matrix(ROOT / cfg["matrix"], set(all_genes))
    coords = load_coordinates(modality, cfg, barcodes)
    positions = {b: i for i, b in enumerate(barcodes)}
    keep = np.asarray([positions[x] for x in coords.index], dtype=int)
    features = np.asarray([g.upper() for g in feature_names])
    target_original = np.flatnonzero(np.isin(features, all_genes))
    row_lookup = {int(original): i for i, original in enumerate(target_original)}
    norm_rows = {}
    for gene in all_genes:
        original = np.flatnonzero(features == gene)
        rows = [row_lookup[int(i)] for i in original if int(i) in row_lookup]
        if rows:
            raw = np.asarray(target[rows, :].sum(axis=0)).ravel().astype(float)
            norm_rows[gene] = np.log1p(raw / np.maximum(totals, 1.0) * 1e4)
    table = coords.copy()
    table["total_counts"] = totals[keep]
    table["n_genes_detected"] = detected[keep]
    source_activity = {}
    for source, genes in SOURCE_MARKERS.items():
        present = [g for g in genes if g in norm_rows and g not in PROTEASE_GENES + ANTIPROTEASE_GENES]
        if present:
            activity = np.nanmean(np.vstack([norm_rows[g] for g in present]), axis=0)
            source_activity[source] = activity[keep]
            table[f"source_{source}"] = activity[keep]
        else:
            source_activity[source] = np.full(len(barcodes), np.nan)
            table[f"source_{source}"] = np.nan
    protease_gene_scores = []
    antiprotease_gene_scores = []
    for gene in PROTEASE_GENES:
        if gene in norm_rows:
            table[f"gene_{gene}"] = norm_rows[gene][keep]
            protease_gene_scores.append(zscore(norm_rows[gene][keep]))
    for gene in ANTIPROTEASE_GENES:
        if gene in norm_rows:
            table[f"gene_{gene}"] = norm_rows[gene][keep]
            antiprotease_gene_scores.append(zscore(norm_rows[gene][keep]))
    table["protease_activity"] = np.nanmean(np.vstack(protease_gene_scores), axis=0) if protease_gene_scores else np.nan
    table["antiprotease_activity"] = np.nanmean(np.vstack(antiprotease_gene_scores), axis=0) if antiprotease_gene_scores else np.nan
    table["proteolytic_balance_z"] = table["protease_activity"] - table["antiprotease_activity"]
    source_cols = [f"source_{x}" for x in SOURCE_MARKERS]
    source_matrix = table[source_cols].to_numpy(float)
    source_matrix = np.maximum(source_matrix, 0)
    source_denominator = np.nansum(source_matrix, axis=1)
    for i, source in enumerate(SOURCE_MARKERS):
        table[f"source_fraction_{source}"] = np.divide(source_matrix[:, i], source_denominator, out=np.full(len(table), np.nan), where=source_denominator > 0)
    table["modality"] = cfg["assay"]
    table["patient_id"] = cfg["patient_id"]
    table["sample_id"] = cfg["sample_id"]
    table["unit_type"] = cfg["unit"]
    table.to_csv(out / f"{modality}_source_attribution_and_balance.csv", index=False)
    # Candidate local association summaries. These are not pooled across
    # technologies and are not communication claims.
    ids = table.index
    rng = np.random.default_rng(seed)
    if len(ids) > max_units:
        ids = pd.Index(rng.choice(ids.to_numpy(), size=max_units, replace=False))
    local = table.loc[ids].copy()
    local_coords = local[["x", "y"]].to_numpy(float)
    valid = np.isfinite(local_coords).all(axis=1)
    local = local.iloc[np.flatnonzero(valid)]
    local_coords = local_coords[valid]
    rows = []
    if len(local) >= 10:
        nn = NearestNeighbors(n_neighbors=min(7, len(local))).fit(local_coords)
        _, neighbor_idx = nn.kneighbors(local_coords)
        neighbor_idx = neighbor_idx[:, 1:]
        targets = ["protease_activity", "antiprotease_activity", "proteolytic_balance_z"]
        for source in source_cols:
            for target_name in targets:
                x = local[source].to_numpy(float)
                y = local[target_name].to_numpy(float)
                neighbor_values = y[neighbor_idx]
                neighbor_counts = np.isfinite(neighbor_values).sum(axis=1)
                neigh_y = np.full(len(local), np.nan, dtype=float)
                valid_neighbor_rows = neighbor_counts > 0
                if valid_neighbor_rows.any():
                    neigh_y[valid_neighbor_rows] = (
                        np.nansum(neighbor_values[valid_neighbor_rows], axis=1)
                        / neighbor_counts[valid_neighbor_rows]
                    )
                valid_xy = np.isfinite(x) & np.isfinite(y)
                corr = float(pd.Series(x[valid_xy]).corr(pd.Series(y[valid_xy]), method="spearman")) if valid_xy.sum() >= 10 else np.nan
                rows.append({
                    "modality": cfg["assay"],
                    "source_activity": source,
                    "target": target_name,
                    "n_units": int(len(local)),
                    "spearman_same_unit": corr,
                    "spearman_source_vs_neighbor_target": float(pd.Series(x[np.isfinite(neigh_y)]).corr(pd.Series(neigh_y[np.isfinite(neigh_y)]), method="spearman")) if np.isfinite(neigh_y).sum() >= 10 else np.nan,
                    "interpretation": "source attribution and proteolytic-balance association; not proof of communication or lesion causality",
                })
    pd.DataFrame(rows).to_csv(out / f"{modality}_source_balance_associations.csv", index=False)
    return {"modality": modality, "patient_id": cfg["patient_id"], "units": int(len(table)), "source_output": str((out / f"{modality}_source_attribution_and_balance.csv").relative_to(ROOT)), "lesion_edge_status": "unavailable_no_reproducible_lesion_mask"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spatial-config", default="config/spatial_analysis.yaml")
    parser.add_argument("--reference", default="data/processed/reproduction_core/GSE135851_core_reproduction.h5ad")
    parser.add_argument("--output-dir", default="results/spatial/GSE302356")
    parser.add_argument("--max-neighbor-units", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260823)
    args = parser.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load((ROOT / args.spatial_config).read_text())
    source_reference, ref_manifest = load_core_source_reference(ROOT / args.reference, out)
    modality_results = [spatial_scores(name, cfg, out, source_reference, args.max_neighbor_units, args.seed) for name, cfg in config["modalities"].items()]
    manifest = {
        "spatial_config": args.spatial_config,
        "reference": args.reference,
        "reference_manifest": ref_manifest,
        "modalities": modality_results,
        "protease_genes": PROTEASE_GENES,
        "antiprotease_genes": ANTIPROTEASE_GENES,
        "balance_definition": "within-modality standardized protease activity minus standardized antiprotease activity; no score ratio",
        "source_attribution_rule": "source identity is estimated from non-protease markers in actual single-cell data; protease genes are measured after labeling",
        "lesion_rule": "No formal lesion-edge claim without a reproducible lesion/cyst mask.",
        "modality_rule": "Visium, Visium HD and Xenium outputs remain separate; same PatientID is orthogonal evidence, not an independent donor.",
    }
    (out / "source_attribution_analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
