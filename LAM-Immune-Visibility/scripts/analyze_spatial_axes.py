"""Analyze immune-visibility axes in GSE302356 spatial modalities separately."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from common import PROJECT_ROOT, ensure_output_path, load_project_config, load_signatures, load_source_manifest, load_yaml, resolve_source, write_json


SOURCE_SCRIPTS = resolve_source("../LAM-Cell-Research/scripts", ".")
sys.path.insert(0, str(SOURCE_SCRIPTS))
from analyze_spatial_niche import load_coordinates, read_10x_target_matrix  # noqa: E402


def modality_paths(cfg: dict, source_project_root: Path) -> dict:
    result = dict(cfg)
    for key in ["matrix", "positions", "cells", "transcripts", "cell_boundaries", "image", "scalefactors"]:
        if key in result:
            value = Path(result[key])
            result[key] = str(value if value.is_absolute() else source_project_root / value)
    return result


def score_modality(name: str, cfg: dict, signatures: dict, source_project_root: Path, seed: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cfg = modality_paths(cfg, source_project_root)
    module_genes = {module: [str(g).upper() for g in spec.get("genes", [])] for module, spec in signatures.items()}
    module_genes["source_immune"] = sorted(set(signatures["t_cell_markers"]["genes"] + signatures["nk_markers"]["genes"] + signatures["macrophage_markers"]["genes"]))
    target_genes = sorted(set(sum(module_genes.values(), [])))
    feature_names, barcodes, target, totals, detected = read_10x_target_matrix(Path(cfg["matrix"]), set(target_genes))
    coords = load_coordinates(name, cfg, barcodes)
    barcode_index = {barcode: index for index, barcode in enumerate(barcodes)}
    keep = np.asarray([barcode_index[unit] for unit in coords.index], dtype=int)
    gene_upper = np.asarray([gene.upper() for gene in feature_names])
    selected_original = np.flatnonzero(np.isin(gene_upper, sorted(set(target_genes))))
    target_row = {int(original): index for index, original in enumerate(selected_original)}
    table = coords.copy()
    table["total_counts"] = totals[keep]
    table["n_genes_detected"] = detected[keep]
    availability_rows = []
    unit_vectors: dict[str, np.ndarray] = {}
    for gene in target_genes:
        original = np.flatnonzero(gene_upper == gene)
        rows = [target_row[int(index)] for index in original if int(index) in target_row]
        if rows:
            raw = np.asarray(target[rows, :].sum(axis=0)).ravel().astype(float)
            norm = np.log1p(raw / np.maximum(totals, 1.0) * 1e4)[keep]
            unit_vectors[gene] = norm
            positive = norm > 0
            low_cut = float(np.quantile(norm[positive], 0.25)) if positive.sum() >= 20 else np.nan
            high_cut = float(np.quantile(norm[positive], 0.75)) if positive.sum() >= 20 else np.nan
            availability_rows.append({
                "modality": cfg["assay"],
                "patient_id": cfg["patient_id"],
                "gene": gene,
                "panel_gene_available": True,
                "n_units": len(norm),
                "n_detected_units": int(positive.sum()),
                "not_detected_fraction": float((~positive).mean()),
                "detected_low_fraction": float(((norm > 0) & (norm < low_cut)).mean()) if np.isfinite(low_cut) else np.nan,
                "detected_high_fraction": float((norm >= high_cut).mean()) if np.isfinite(high_cut) else np.nan,
                "low_expression_upper_bound": low_cut,
                "high_expression_lower_bound": high_cut,
                "state_note": "nonzero transcript means detected; low/high are expression intervals",
            })
        else:
            availability_rows.append({"modality": cfg["assay"], "patient_id": cfg["patient_id"], "gene": gene, "panel_gene_available": False, "n_units": len(keep), "state_note": "not_assayed"})
    for module, genes in module_genes.items():
        vectors = [unit_vectors[gene] for gene in genes if gene in unit_vectors]
        table[f"module_{module}"] = np.nanmean(np.vstack(vectors), axis=0) if vectors else np.nan
        table[f"module_{module}__n_available"] = sum(gene in unit_vectors for gene in genes)
        table[f"module_{module}__status"] = "analyzable" if vectors else "not_assayed"

    ids = table.index
    rng = np.random.default_rng(seed)
    if len(ids) > 50000:
        ids = pd.Index(rng.choice(ids.to_numpy(), size=50000, replace=False))
    local = table.loc[ids].copy()
    local_coords = local[["x", "y"]].astype(float).to_numpy()
    valid = np.isfinite(local_coords).all(axis=1)
    local = local.iloc[np.flatnonzero(valid)]
    local_coords = local_coords[valid]
    rows = []
    if len(local) >= 10:
        nn = NearestNeighbors(n_neighbors=min(7, len(local))).fit(local_coords)
        _, neighbor_indices = nn.kneighbors(local_coords)
        neighbor_indices = neighbor_indices[:, 1:]
        source = pd.to_numeric(local["module_identity_protected"], errors="coerce").to_numpy(float)
        for target_module in ["source_immune", "t_cell_exhaustion", "nk_state", "macrophage_suppressive", "immune_evasion"]:
            target_values = pd.to_numeric(local[f"module_{target_module}"], errors="coerce").to_numpy(float)
            source_high = source >= np.nanquantile(source, 0.9)
            target_high = target_values >= np.nanquantile(target_values, 0.9)
            neighbor_mean = np.nanmean(target_values[neighbor_indices], axis=1)
            target_neighbor_high = neighbor_mean >= np.nanquantile(neighbor_mean[np.isfinite(neighbor_mean)], 0.9) if np.isfinite(neighbor_mean).any() else np.zeros(len(neighbor_mean), dtype=bool)
            usable = source_high & np.isfinite(neighbor_mean)
            rows.append({
                "modality": cfg["assay"],
                "patient_id": cfg["patient_id"],
                "source": "identity_protected_LAMCORE_like_top_decile",
                "target": target_module,
                "n_units": len(local),
                "source_units": int(source_high.sum()),
                "target_neighbor_high_rate": float(target_neighbor_high[usable].mean()) if usable.any() else np.nan,
                "target_global_high_rate": float(target_high[np.isfinite(target_values)].mean()) if np.isfinite(target_values).any() else np.nan,
                "interpretation": "spatial candidate association; not proof of communication",
            })
    return table, pd.DataFrame(availability_rows), pd.DataFrame(rows)


def main() -> None:
    config = load_project_config()
    manifest = load_source_manifest()
    spatial_config_path = resolve_source(manifest["spatial"]["config"], manifest["source_root"])
    spatial_config = load_yaml(spatial_config_path)
    source_project_root = spatial_config_path.parent.parent
    signatures = load_signatures()
    output_dir = PROJECT_ROOT / "results" / "spatial"
    all_availability = []
    all_associations = []
    summaries = []
    for name, cfg in spatial_config["modalities"].items():
        units, availability, associations = score_modality(name, cfg, signatures, source_project_root, config["random_seed"])
        units.to_csv(ensure_output_path(output_dir / f"{name}_visibility_unit_scores.csv"), index=False)
        availability.to_csv(ensure_output_path(output_dir / f"{name}_gene_availability.csv"), index=False)
        associations.to_csv(ensure_output_path(output_dir / f"{name}_visibility_associations.csv"), index=False)
        all_availability.append(availability)
        all_associations.append(associations)
        summaries.append({"modality": name, "patient_id": cfg["patient_id"], "assay": cfg["assay"], "n_units": len(units), "n_association_rows": len(associations)})
    pd.concat(all_availability, ignore_index=True).to_csv(ensure_output_path(output_dir / "gene_availability_summary.csv"), index=False)
    pd.concat(all_associations, ignore_index=True).to_csv(ensure_output_path(output_dir / "visibility_association_summary.csv"), index=False)
    write_json(PROJECT_ROOT / "manifests" / "spatial_visibility_manifest.json", {
        "modalities": summaries,
        "same_patient_modalities_are_not_independent": True,
        "missing_panel_rule": "panel-unobserved is not biological negative",
        "interpretation": "Spatial axes are candidate associations only.",
    })
    print(pd.DataFrame(summaries).to_string(index=False))


if __name__ == "__main__":
    main()
