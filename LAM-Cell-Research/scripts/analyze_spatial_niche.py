#!/usr/bin/env python3
"""Analyze the GSE302356 spatial modalities without pooling their units.

This is an evidence-generation pass for the protease/ECM niche question. The
three technologies are deliberately processed separately:

* Visium: spot-level
* Visium HD: bin/segment-level
* Xenium: cell-level, targeted panel

Scores are simple library-size-normalized module summaries. They are not cell
type calls and are not treated as equivalent measurements across assays.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]

PROGRAMS: dict[str, list[str]] = {
    "LAMCORE": ["PMEL", "MLANA", "ACTA2", "FIGF", "VEGFD", "CTSK", "MITF", "ESR1"],
    "LAF_fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "CFD", "C7", "CXCL14", "CCL2", "TGFB1"],
    "lymphatic_endothelial": ["PDPN", "LYVE1", "FLT4", "CCL21", "PROX1", "PECAM1", "EMCN", "KDR"],
    "immune": ["PTPRC", "LST1", "AIF1", "TYROBP", "FCER1G", "CTSS", "CD68", "CD3D", "CD3E", "S100A8", "S100A9"],
    "protease": ["MMP2", "MMP8", "MMP9", "MMP12", "CTSK", "CTSB", "ELANE", "CTSS"],
    "antiprotease": ["TIMP1", "TIMP2", "SERPINA1", "SERPINE1", "SLPI", "PI3"],
    "ECM_remodeling": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "FN1", "VCAN", "TNC"],
}


def decode(values: np.ndarray) -> list[str]:
    return [x.decode() if isinstance(x, (bytes, np.bytes_)) else str(x) for x in values]


def read_10x_target_matrix(path: Path, target_genes: set[str]) -> tuple[list[str], list[str], sp.csc_matrix, np.ndarray, np.ndarray]:
    """Read only the target feature rows while retaining full library sizes."""
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        feature_names = decode(matrix["features"]["name"][:])
        barcodes = decode(matrix["barcodes"][:])
        data = matrix["data"][:]
        indices = matrix["indices"][:]
        indptr = matrix["indptr"][:]
        shape = tuple(int(x) for x in matrix["shape"][:])
        full = sp.csc_matrix((data, indices, indptr), shape=shape)
        gene_upper = np.asarray([g.upper() for g in feature_names])
        selected = np.flatnonzero(np.isin(gene_upper, sorted(target_genes)))
        target = full[selected, :].tocsc()
        totals = np.asarray(full.sum(axis=0)).ravel().astype(float)
        detected = np.asarray(full.getnnz(axis=0)).ravel().astype(int)
    return feature_names, barcodes, target, totals, detected


def load_coordinates(modality: str, config: dict, barcodes: list[str]) -> pd.DataFrame:
    if modality in {"visium", "visium_hd"}:
        path = ROOT / config["positions"]
        coords = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)
        coords["barcode"] = coords["barcode"].astype(str)
        coords = coords[coords["barcode"].isin(set(barcodes))].copy()
        if "in_tissue" in coords.columns:
            coords = coords[coords["in_tissue"].astype(int).eq(1)].copy()
        coords = coords.rename(columns={"barcode": "unit_id", "pxl_col_in_fullres": "x", "pxl_row_in_fullres": "y"})
    else:
        coords = pd.read_parquet(ROOT / config["cells"])
        coords["cell_id"] = coords["cell_id"].astype(str)
        coords = coords[coords["cell_id"].isin(set(barcodes))].copy()
        coords = coords.rename(columns={"cell_id": "unit_id", "x_centroid": "x", "y_centroid": "y"})
    return coords.set_index("unit_id", drop=False)


def module_score(target: sp.csc_matrix, feature_names: list[str], barcodes: list[str], total_counts: np.ndarray, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    gene_upper = np.asarray([g.upper() for g in feature_names])
    present = [g for g in genes if np.any(gene_upper == g)]
    if not present:
        return np.full(len(barcodes), np.nan), []
    rows = np.flatnonzero(np.isin(gene_upper, present))
    # target contains the union of all requested rows. Mapping back to the
    # target row indices keeps this calculation sparse and memory bounded.
    target_genes = gene_upper[np.flatnonzero(np.isin(gene_upper, sorted(set(sum(PROGRAMS.values(), [])))))]
    # The target rows are in the same order as the sorted union selected above.
    selected_all = np.flatnonzero(np.isin(gene_upper, sorted(set(sum(PROGRAMS.values(), [])))))
    row_lookup = {int(original): i for i, original in enumerate(selected_all)}
    target_rows = [row_lookup[int(r)] for r in rows if int(r) in row_lookup]
    if not target_rows:
        return np.full(len(barcodes), np.nan), []
    raw = np.asarray(target[target_rows, :].sum(axis=0)).ravel().astype(float)
    denom = np.maximum(total_counts, 1.0) * len(target_rows)
    score = np.log1p(raw / denom * 1e4)
    return score, present


def pairwise_spatial_summary(scores: pd.DataFrame, coords: pd.DataFrame, modality: str, max_units: int, seed: int) -> pd.DataFrame:
    program_names = [x for x in PROGRAMS if x in scores.columns]
    if len(program_names) < 2 or len(scores) < 10:
        return pd.DataFrame()
    rng = np.random.default_rng(seed)
    ids = scores.index.intersection(coords.index)
    if len(ids) > max_units:
        ids = pd.Index(rng.choice(ids.to_numpy(), size=max_units, replace=False))
    local_scores = scores.loc[ids, program_names]
    local_coords = coords.loc[ids, ["x", "y"]].astype(float)
    valid = np.isfinite(local_coords.to_numpy()).all(axis=1)
    local_scores = local_scores.iloc[np.flatnonzero(valid)]
    local_coords = local_coords.iloc[np.flatnonzero(valid)]
    if len(local_scores) < 10:
        return pd.DataFrame()
    neighbors = NearestNeighbors(n_neighbors=min(7, len(local_scores)), algorithm="auto").fit(local_coords.to_numpy())
    _, indices = neighbors.kneighbors(local_coords.to_numpy())
    indices = indices[:, 1:]
    rows = []
    for source in program_names:
        source_values = local_scores[source].to_numpy(float)
        source_high = source_values >= np.nanquantile(source_values, 0.9)
        for target in program_names:
            if source == target:
                continue
            target_values = local_scores[target].to_numpy(float)
            finite = np.isfinite(source_values) & np.isfinite(target_values)
            corr = float(pd.Series(source_values[finite]).corr(pd.Series(target_values[finite]), method="spearman")) if finite.sum() >= 10 else np.nan
            neigh_mean = np.nanmean(target_values[indices], axis=1)
            target_high = target_values >= np.nanquantile(target_values, 0.9)
            source_high_valid = source_high & np.isfinite(neigh_mean)
            observed = float(target_high[source_high_valid].mean()) if source_high_valid.any() else np.nan
            global_rate = float(target_high[np.isfinite(target_values)].mean()) if np.isfinite(target_values).any() else np.nan
            rows.append({
                "modality": modality,
                "source_program": source,
                "target_program": target,
                "n_units_analyzed": int(len(local_scores)),
                "n_neighbors": int(indices.shape[1]),
                "spearman_unit_score": corr,
                "source_top_decile_target_neighbor_top_decile_rate": observed,
                "target_top_decile_global_rate": global_rate,
                "neighbor_enrichment_ratio": observed / global_rate if global_rate and np.isfinite(observed) else np.nan,
                "interpretation": "spatial co-localization candidate; not proof of cellular communication",
            })
    return pd.DataFrame(rows)


def analyze_modality(modality: str, config: dict, out: Path, max_units: int, seed: int) -> dict:
    matrix_path = ROOT / config["matrix"]
    all_genes, barcodes, target, totals, detected = read_10x_target_matrix(matrix_path, set(sum(PROGRAMS.values(), [])))
    coords = load_coordinates(modality, config, barcodes)
    barcode_pos = {b: i for i, b in enumerate(barcodes)}
    keep = np.asarray([barcode_pos[x] for x in coords.index], dtype=int)
    unit_scores = coords.copy()
    unit_scores["total_counts"] = totals[keep]
    unit_scores["n_genes_detected"] = detected[keep]
    availability: dict[str, list[str]] = {}
    for name, genes in PROGRAMS.items():
        score, present = module_score(target, all_genes, barcodes, totals, genes)
        unit_scores[name] = score[keep]
        availability[name] = present
    unit_scores["modality"] = config["assay"]
    unit_scores["patient_id"] = config["patient_id"]
    unit_scores["sample_id"] = config["sample_id"]
    unit_scores["unit_type"] = config["unit"]
    unit_scores.to_csv(out / f"{modality}_unit_scores.csv", index=False)
    pairwise = pairwise_spatial_summary(unit_scores.set_index("unit_id"), coords, config["assay"], max_units, seed)
    pairwise.to_csv(out / f"{modality}_co_localization.csv", index=False)
    manifest = {
        "dataset": "GSE302356",
        "modality": modality,
        "sample_id": config["sample_id"],
        "patient_id": config["patient_id"],
        "assay": config["assay"],
        "unit": config["unit"],
        "n_matrix_units": len(barcodes),
        "n_units_analyzed_after_position_filter": len(unit_scores),
        "program_genes_present": availability,
        "panel_unobserved_rule": "A missing gene in Xenium is panel-unobserved, not a biological negative.",
        "cross_modality_rule": "Do not pool raw units, scores, or p-values across Visium, Visium HD, and Xenium.",
        "co_localization_note": "Nearest-neighbor enrichment and score correlation are candidate niche evidence only.",
        "max_units_for_neighbor_summary": max_units,
        "seed": seed,
    }
    (out / f"{modality}_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    return {"modality": modality, "units": len(unit_scores), "matrix_units": len(barcodes), "co_localization_rows": len(pairwise)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/spatial_analysis.yaml")
    parser.add_argument("--output-dir", default="results/spatial/GSE302356")
    parser.add_argument("--max-neighbor-units", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20260823)
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / args.config).read_text())
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    results = [analyze_modality(name, cfg, out, args.max_neighbor_units, args.seed) for name, cfg in config["modalities"].items()]
    run_manifest = {
        "config": args.config,
        "programs": PROGRAMS,
        "modalities": results,
        "analysis_rule": config["analysis_rules"],
        "interpretation": "The modalities support evidence-level triangulation, not pooled quantitative inference.",
    }
    (out / "spatial_analysis_manifest.json").write_text(json.dumps(run_manifest, ensure_ascii=False, indent=2))
    print(json.dumps(run_manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
