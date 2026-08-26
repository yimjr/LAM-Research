#!/usr/bin/env python3
"""Compare lineage/reference and lung-adaptation programs across tissues.

This pass is deliberately a reference comparison, not a proof of cellular
origin. LAMCORE-like cells are taken from the conservative fallback candidate
label in the core reanalysis; normal lung and mouse uterus are reference
cohorts, not universal negative controls.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

PROGRAMS = {
    "lineage_origin_smooth_muscle": ["ACTA2", "TAGLN", "MYH11", "CNN1", "DES", "CALD1", "MYL9", "TPM2"],
    "lung_adaptation": ["VEGFD", "FIGF", "CCL21", "PDPN", "LYVE1", "FLT4", "CTSK", "MMP2", "MMP9", "CXCL14"],
    "ECM_interaction": ["COL1A1", "COL1A2", "COL3A1", "FN1", "VCAN", "MMP2", "MMP9", "TIMP1", "TIMP2"],
    "survival_stress": ["BCL2", "BCL2L1", "LGALS3", "HIF1A", "VEGFA", "BNIP3", "DDIT4", "NDRG1", "TXNIP"],
    "immune_evasion": ["LGALS3", "CTSS", "TYROBP", "CCL2", "TGFB1", "CD274", "SERPINE1"],
    "hormone_related": ["ESR1", "ESR2", "PGR", "IGF1R", "EGFR", "GREB1"],
}


def gene_symbols(adata: ad.AnnData) -> np.ndarray:
    for column in ["gene_symbol_upper", "gene_symbol"]:
        if column in adata.var:
            return adata.var[column].astype(str).str.upper().to_numpy()
    return adata.var_names.astype(str).str.upper().to_numpy()


def row_totals(adata: ad.AnnData, chunk_size: int = 2000) -> np.ndarray:
    if "total_counts" in adata.obs:
        values = pd.to_numeric(adata.obs["total_counts"], errors="coerce").to_numpy(float)
        if np.isfinite(values).all() and np.all(values > 0):
            return values
    totals = np.zeros(adata.n_obs, dtype=float)
    for start in range(0, adata.n_obs, chunk_size):
        block = adata.X[start:start + chunk_size]
        totals[start:start + block.shape[0]] = np.asarray(block.sum(axis=1)).ravel()
    return totals


def score_dataset(path: Path, cohort: str, candidate_label: str | None = None, candidate_values: set[str] | None = None) -> tuple[pd.DataFrame, dict]:
    adata = ad.read_h5ad(path, backed="r")
    symbols = gene_symbols(adata)
    totals = row_totals(adata)
    out = pd.DataFrame(index=adata.obs_names.astype(str))
    out["cohort"] = cohort
    out["sample_id"] = adata.obs["sample_id"].astype(str).to_numpy() if "sample_id" in adata.obs else cohort
    out["donor_id"] = adata.obs["donor_id"].astype(str).to_numpy() if "donor_id" in adata.obs else cohort
    out["tissue"] = adata.obs["tissue"].astype(str).to_numpy() if "tissue" in adata.obs else cohort
    out["assay"] = adata.obs["assay"].astype(str).to_numpy() if "assay" in adata.obs else "unknown"
    out["total_counts"] = totals
    if candidate_label and candidate_label in adata.obs:
        values = adata.obs[candidate_label].astype(str)
        out["lamcore_candidate"] = values.isin(candidate_values or {"True", "true", "1"}).to_numpy()
        # The normal lung reference donor in the core object is named Donor1;
        # it must not enter the LAMCORE candidate cohort.
        if cohort == "LAM_lung_core_reanalysis":
            out.loc[~out["donor_id"].isin({"LAM1", "LAM2", "LAM3", "LAM4"}), "lamcore_candidate"] = False
    else:
        out["lamcore_candidate"] = False
    available: dict[str, list[str]] = {}
    for name, genes in PROGRAMS.items():
        selected = [g for g in genes if g in set(symbols)]
        available[name] = selected
        if not selected:
            out[name] = np.nan
            continue
        indices = np.flatnonzero(np.isin(symbols, selected))
        block = adata[:, indices]
        if "counts" in adata.layers:
            matrix = block.layers["counts"]
        else:
            matrix = block.X
        raw = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
        raw = np.asarray(raw, dtype=float)
        summed = raw.sum(axis=1)
        out[name] = np.log1p(summed / (np.maximum(totals, 1.0) * len(indices)) * 1e4)
    out["cell_id"] = out.index
    manifest = {"path": str(path.relative_to(ROOT)), "cohort": cohort, "n_cells": len(out), "available_genes": available}
    return out.reset_index(drop=True), manifest


def donor_summary(cells: pd.DataFrame) -> pd.DataFrame:
    program_names = list(PROGRAMS)
    rows = []
    for (cohort, donor_id), group in cells.groupby(["cohort", "donor_id"], observed=True):
        for label, subset in [("all_cells", group), ("lamcore_candidate", group[group["lamcore_candidate"]])]:
            if subset.empty:
                continue
            row = {"cohort": cohort, "donor_id": donor_id, "cell_subset": label, "n_cells": len(subset)}
            for program in program_names:
                row[f"{program}_mean"] = float(subset[program].mean()) if subset[program].notna().any() else np.nan
                row[f"{program}_median"] = float(subset[program].median()) if subset[program].notna().any() else np.nan
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/adaptation")
    args = parser.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    datasets = [
        ("data/processed/GSE135851_lam_states.h5ad", "LAM_lung_core_reanalysis", "lamcore_candidate_fallback", {"True", "true", "1"}),
        ("data/processed/external/GSE122960_normal_lung.h5ad", "normal_lung_reference", None, None),
        ("data/processed/external/GSE118180_wildtype_uterus.h5ad", "mouse_uterus_lineage_reference", None, None),
    ]
    frames = []
    manifests = []
    for path, cohort, label, values in datasets:
        frame, manifest = score_dataset(ROOT / path, cohort, label, values)
        frames.append(frame)
        manifests.append(manifest)
    cells = pd.concat(frames, ignore_index=True)
    cells.to_csv(out / "adaptation_cell_scores.csv", index=False)
    donors = donor_summary(cells)
    donors.to_csv(out / "adaptation_donor_program_scores.csv", index=False)

    # Evidence-level contrasts are descriptive and use donor summaries as the
    # unit. They do not claim cell-matched or causal lineage relationships.
    candidate = donors[donors["cell_subset"].eq("lamcore_candidate")]
    normal_lung = donors[(donors["cohort"].eq("normal_lung_reference")) & donors["cell_subset"].eq("all_cells")]
    uterus = donors[(donors["cohort"].eq("mouse_uterus_lineage_reference")) & donors["cell_subset"].eq("all_cells")]
    contrasts = []
    for program in PROGRAMS:
        col = f"{program}_median"
        lam = float(candidate[col].median()) if col in candidate and candidate[col].notna().any() else np.nan
        lung = float(normal_lung[col].median()) if col in normal_lung and normal_lung[col].notna().any() else np.nan
        uter = float(uterus[col].median()) if col in uterus and uterus[col].notna().any() else np.nan
        contrasts.append({
            "program": program,
            "lamcore_candidate_median_across_donors": lam,
            "normal_lung_reference_median_across_donors": lung,
            "mouse_uterus_reference_median": uter,
            "lam_vs_normal_lung_difference": lam - lung if np.isfinite(lam) and np.isfinite(lung) else np.nan,
            "lam_vs_mouse_uterus_difference": lam - uter if np.isfinite(lam) and np.isfinite(uter) else np.nan,
            "interpretation": "descriptive origin/adaptation comparison; uterus is a lineage reference, not a negative control",
        })
    pd.DataFrame(contrasts).to_csv(out / "adaptation_contrasts.csv", index=False)
    manifest = {
        "programs": PROGRAMS,
        "datasets": manifests,
        "candidate_rule": "LAM1-LAM4 cells with lamcore_candidate_fallback=True; the broad formal candidate was not used as the conservative identity set.",
        "reference_rule": "Normal lung is a disease-control reference; mouse uterus is a lineage reference; neither is a universal exclusion criterion.",
        "species_rule": "Mouse uterus is reported as a symbol-overlap reference; formal ortholog mapping remains a follow-up if this line is upgraded.",
        "interpretation": "Scores are evidence for comparative programs, not proof of cell origin or a lung-specific causal mechanism.",
    }
    (out / "adaptation_analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(json.dumps({"cells": len(cells), "donor_rows": len(donors), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
