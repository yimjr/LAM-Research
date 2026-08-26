"""Score GSE302356 raw 10x samples against LAM residual programs.

This is a modality-aware preliminary human validation. The downloaded subset
contains LAM3/LAM4 scRNA-seq and LAM18/LAM20 spatial matrices, but no formal
cell labels. Paper-derived operational marker panels are scored as enrichment;
they are not treated as formal state assignments, especially for LAMCORE3.
"""

from __future__ import annotations

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml

from common import ROOT, write_json


SAMPLE_META = {
    "LAM3": {"modality": "scRNA-seq", "sample_accession": "GSM9102552"},
    "LAM4": {"modality": "scRNA-seq", "sample_accession": "GSM9102553"},
    "LAM18": {"modality": "Visium HD", "sample_accession": "GSM9226654"},
    "LAM20": {"modality": "Visium", "sample_accession": "GSM9226656"},
}


def load_programs() -> tuple[dict[str, dict[str, list[str]]], dict[str, list[str]], dict[str, list[str]]]:
    signatures = pd.read_csv(ROOT / "results/signatures/GSE179044_cmap_query_signatures.csv")
    signed: dict[str, dict[str, list[str]]] = {}
    for contrast, group in signatures.groupby("contrast"):
        signed[contrast] = {
            "up": group.loc[group.direction.eq("up"), "gene"].astype(str).str.upper().unique().tolist(),
            "down": group.loc[group.direction.eq("down"), "gene"].astype(str).str.upper().unique().tolist(),
        }
    config = yaml.safe_load((ROOT / "config/analysis.yaml").read_text())
    modules = {name: [str(gene).upper() for gene in genes] for name, genes in config["module_sets"].items()}
    state_path = ROOT / "data/processed/GSE302356/paper_state_marker_panels.csv"
    state_table = pd.read_csv(state_path)
    state_panels = {
        state: group["marker"].astype(str).str.upper().unique().tolist()
        for state, group in state_table.groupby("state", sort=False)
    }
    return signed, modules, state_panels


def gene_index(adata: ad.AnnData) -> dict[str, int]:
    result: dict[str, int] = {}
    for i, gene in enumerate(adata.var_names.astype(str)):
        result.setdefault(gene.upper(), i)
    return result


def standardize_scores(x: sp.csr_matrix, indices: list[int], signs: np.ndarray | None = None) -> np.ndarray:
    if not indices:
        return np.full(x.shape[0], np.nan)
    subset = x[:, indices].tocsr()
    mean = np.asarray(subset.mean(axis=0)).ravel()
    second = np.asarray(subset.multiply(subset).mean(axis=0)).ravel()
    std = np.sqrt(np.maximum(second - mean * mean, 1e-8))
    if signs is None:
        signs = np.ones(len(indices), dtype=float)
    weights = signs.astype(float) / std / len(indices)
    score = np.asarray(subset @ weights).ravel()
    return score - float(np.sum(mean * weights))


def signed_signature_score(x: sp.csr_matrix, lookup: dict[str, int], up: list[str], down: list[str]) -> tuple[np.ndarray, int, int]:
    up_indices = [lookup[g] for g in up if g in lookup]
    down_indices = [lookup[g] for g in down if g in lookup]
    up_score = standardize_scores(x, up_indices)
    down_score = standardize_scores(x, down_indices)
    return up_score - down_score, len(up_indices), len(down_indices)


def raw_signature_score(x: sp.csr_matrix, lookup: dict[str, int], up: list[str], down: list[str]) -> tuple[np.ndarray, int, int]:
    up_indices = [lookup[g] for g in up if g in lookup]
    down_indices = [lookup[g] for g in down if g in lookup]
    up_score = np.asarray(x[:, up_indices].mean(axis=1)).ravel() if up_indices else np.zeros(x.shape[0])
    down_score = np.asarray(x[:, down_indices].mean(axis=1)).ravel() if down_indices else np.zeros(x.shape[0])
    return up_score - down_score, len(up_indices), len(down_indices)


def process_sample(
    path: Path,
    signed: dict[str, dict[str, list[str]]],
    modules: dict[str, list[str]],
    state_panels: dict[str, list[str]],
) -> tuple[pd.DataFrame, dict]:
    sample_id = path.parent.name
    adata = sc.read_10x_h5(path)
    adata.var_names = adata.var_names.astype(str)
    x = adata.X.tocsr().astype(np.float32)
    total_counts = np.asarray(x.sum(axis=1)).ravel()
    detected = np.asarray((x > 0).sum(axis=1)).ravel()
    valid = total_counts > 0
    tissue_filter = "not_applicable"
    if SAMPLE_META[sample_id]["modality"].startswith("Visium"):
        position_path = path.parent / f"{sample_id}_tissue_positions.parquet"
        if not position_path.exists():
            position_path = path.parent / f"{sample_id}_tissue_positions.csv"
        positions = pd.read_parquet(position_path) if position_path.suffix == ".parquet" else pd.read_csv(position_path)
        tissue_barcodes = set(positions.loc[positions["in_tissue"].astype(int).eq(1), "barcode"].astype(str))
        valid &= np.array([str(barcode) in tissue_barcodes for barcode in adata.obs_names])
        tissue_filter = f"in_tissue==1 from {position_path.name}"
    n_before_filter = int(valid.size)
    x = x[valid]
    total_counts = total_counts[valid]
    detected = detected[valid]
    x = sp.diags(1e4 / total_counts) @ x
    x = x.tocsr()
    x.data = np.log1p(x.data)
    lookup = gene_index(adata)

    scores = pd.DataFrame({
        "cell_id": adata.obs_names.astype(str)[valid],
        "sample_id": sample_id,
        "modality": SAMPLE_META[sample_id]["modality"],
        "n_counts": total_counts,
        "n_genes": detected,
    })
    availability = {}
    for contrast, sets in signed.items():
        score, n_up, n_down = raw_signature_score(x, lookup, sets["up"], sets["down"])
        z_score, _, _ = signed_signature_score(x, lookup, sets["up"], sets["down"])
        scores[f"signature_{contrast}"] = score
        scores[f"z_signature_{contrast}"] = z_score
        availability[contrast] = {"up": n_up, "down": n_down}
    for module, genes in modules.items():
        indices = [lookup[g] for g in genes if g in lookup]
        scores[f"module_{module}"] = np.asarray(x[:, indices].mean(axis=1)).ravel() if indices else np.full(x.shape[0], np.nan)
        scores[f"z_module_{module}"] = standardize_scores(x, indices)
    state_availability = {}
    for state, genes in state_panels.items():
        indices = [lookup[g] for g in genes if g in lookup]
        scores[f"state_{state}"] = np.asarray(x[:, indices].mean(axis=1)).ravel() if indices else np.full(x.shape[0], np.nan)
        scores[f"z_state_{state}"] = standardize_scores(x, indices)
        state_availability[state] = {"n_panel_genes": len(genes), "n_detected": len(indices)}

    metadata = {
        "sample_id": sample_id,
        "sample_accession": SAMPLE_META[sample_id]["sample_accession"],
        "modality": SAMPLE_META[sample_id]["modality"],
        "n_cells_or_spots": int(len(scores)),
        "n_before_tissue_filter": n_before_filter,
        "n_genes": int(adata.n_vars),
        "median_counts": float(np.median(total_counts)),
        "median_genes": float(np.median(detected)),
        "tissue_filter": tissue_filter,
        "program_gene_availability": availability,
        "state_marker_availability": state_availability,
    }
    return scores, metadata


def main() -> None:
    signed, modules, state_panels = load_programs()
    input_root = ROOT / "data/raw/GSE302356/unpacked"
    paths = [input_root / sample / f"{sample}_filtered_feature_bc_matrix.h5" for sample in SAMPLE_META]
    rows = []
    metadata = []
    for path in paths:
        scores, info = process_sample(path, signed, modules, state_panels)
        rows.append(scores)
        metadata.append(info)
    all_scores = pd.concat(rows, ignore_index=True)
    out_dir = ROOT / "results/human_mapping"
    out_dir.mkdir(parents=True, exist_ok=True)
    all_scores.to_csv(out_dir / "GSE302356_raw_cell_scores.csv.gz", index=False, compression="gzip")

    score_columns = [
        column for column in all_scores.columns
        if column.startswith("signature_") or column.startswith("module_") or column.startswith("state_") or column.startswith("z_state_")
    ]
    summary_rows = []
    for (sample_id, modality), group in all_scores.groupby(["sample_id", "modality"], sort=False):
        row = {"sample_id": sample_id, "modality": modality, "n_cells_or_spots": len(group)}
        for column in score_columns:
            row[f"{column}_mean"] = float(group[column].mean())
            row[f"{column}_median"] = float(group[column].median())
            row[f"{column}_p90"] = float(group[column].quantile(0.9))
        for state in state_panels:
            z_column = f"z_state_{state}"
            row[f"fraction_{z_column}_gt1"] = float((group[z_column] > 1).mean())
        summary_rows.append(row)
    pd.DataFrame(summary_rows).to_csv(out_dir / "GSE302356_raw_sample_scores.csv", index=False)
    write_json(ROOT / "manifests/GSE302356_raw_analysis.json", {
        "status": "preliminary_modality_level_mapping",
        "samples": metadata,
        "state_labels": "Formal cell labels are not present in these raw archives; state_* columns use paper-derived operational marker panels.",
        "state_panel_source": "data/processed/GSE302356/paper_state_marker_panels.csv; based on DOI 10.1183/13993003.02049-2025 and the accessible preprint text. LAMCORE3 has no unique marker panel in the paper and is represented only by a shared-core plus translation-enriched surrogate.",
        "interpretation": "Raw scores are mean log1p-normalized expression of up genes minus down genes (or mean module/state expression); z scores are within-sample diagnostics. These test whether GSE179044 programs and paper-derived state programs are detectable in human LAM scRNA/spatial samples, not formal patient-level state assignments, and cross-modality results are not same-patient tri-modal validation.",
    })
    print({"dataset": "GSE302356", "samples": [(m["sample_id"], m["modality"], m["n_cells_or_spots"]) for m in metadata]})


if __name__ == "__main__":
    main()
