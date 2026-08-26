#!/usr/bin/env python3
"""Four-group lung-adaptation comparison using difference-in-differences logic.

Groups:
normal uterus, LAM uterus, normal lung, pulmonary LAM.

The two uterus matrices are public processed Matrix Market files containing
many empty/raw barcodes. They are filtered by total counts before scoring.
Because the uterus references currently have one source specimen each, the
interaction is reported as descriptive population-level evidence, not a
patient-matched causal estimate.
"""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

ROOT = Path(__file__).resolve().parents[1]

PROGRAMS = {
    "lineage_origin_smooth_muscle": ["ACTA2", "TAGLN", "MYH11", "CNN1", "DES", "CALD1", "MYL9", "TPM2"],
    "lung_adaptation": ["VEGFD", "FIGF", "CCL21", "PDPN", "LYVE1", "FLT4", "CTSK", "MMP2", "MMP9", "CXCL14"],
    "ECM_interaction": ["COL1A1", "COL1A2", "COL3A1", "FN1", "VCAN", "MMP2", "MMP9", "TIMP1", "TIMP2"],
    "survival_stress": ["BCL2", "BCL2L1", "LGALS3", "HIF1A", "VEGFA", "BNIP3", "DDIT4", "NDRG1", "TXNIP"],
    "immune_evasion": ["LGALS3", "CTSS", "TYROBP", "CCL2", "TGFB1", "CD274", "SERPINE1"],
    "hormone_related": ["ESR1", "ESR2", "PGR", "IGF1R", "EGFR", "GREB1"],
}


def score_matrix_market(matrix_path: Path, genes_path: Path, group: str, donor_id: str, min_counts: int) -> dict:
    genes_table = pd.read_csv(genes_path, sep="\t", header=None)
    symbols = genes_table.iloc[:, 1].astype(str).str.upper().to_numpy() if genes_table.shape[1] >= 2 else genes_table.iloc[:, 0].astype(str).str.upper().to_numpy()
    matrix = mmread(matrix_path).tocsr()
    totals = np.asarray(matrix.sum(axis=0)).ravel().astype(float)
    keep = totals >= min_counts
    matrix = matrix[:, keep]
    totals = totals[keep]
    summary = {"group": group, "donor_id": donor_id, "n_cells_after_filter": int(keep.sum()), "min_total_counts": min_counts}
    for program, program_genes in PROGRAMS.items():
        present = [g for g in program_genes if g in set(symbols)]
        if not present:
            summary[f"{program}_mean"] = np.nan
            summary[f"{program}_median"] = np.nan
            summary[f"{program}_pseudobulk_cpm"] = np.nan
            summary[f"{program}_genes_present"] = ""
            continue
        rows = np.flatnonzero(np.isin(symbols, present))
        counts = matrix[rows, :].sum(axis=0)
        counts = np.asarray(counts).ravel().astype(float)
        scores = np.log1p(counts / (np.maximum(totals, 1.0) * len(rows)) * 1e4)
        summary[f"{program}_mean"] = float(np.mean(scores))
        summary[f"{program}_median"] = float(np.median(scores))
        summary[f"{program}_pseudobulk_cpm"] = float(counts.sum() / max(matrix.sum(), 1) * 1e6)
        summary[f"{program}_genes_present"] = ";".join(present)
    del matrix
    gc.collect()
    return summary


def h5ad_group_summary(path: Path, group: str, subset: str | None = None) -> pd.DataFrame:
    cells = pd.read_csv(ROOT / "results/adaptation/adaptation_cell_scores.csv")
    cells = cells[cells["cohort"].eq(group)].copy()
    if subset is not None:
        cells = cells.query(subset).copy()
    rows = []
    for donor_id, donor in cells.groupby("donor_id", observed=True):
        row = {"group": group, "donor_id": str(donor_id), "n_cells_after_filter": len(donor)}
        for program in PROGRAMS:
            row[f"{program}_mean"] = float(donor[program].mean())
            row[f"{program}_median"] = float(donor[program].median())
            row[f"{program}_pseudobulk_cpm"] = np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/adaptation")
    parser.add_argument("--min-total-counts", type=int, default=500)
    args = parser.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    uterus_specs = [
        ("LAM_uterus", "GSE135851_GSM4035470", "data/raw/GSE135851/extracted/GSM4035470_LAM_Uterus.mtx.tsv/Uterus_LAM/matrix.mtx", "data/raw/GSE135851/extracted/GSM4035470_LAM_Uterus.mtx.tsv/Uterus_LAM/genes.tsv"),
        ("normal_uterus", "GSE135851_GSM4035471", "data/raw/GSE135851/extracted/GSM4035471_Normal_Uterus.mtx.tsv/Normal_Uterus/matrix.mtx", "data/raw/GSE135851/extracted/GSM4035471_Normal_Uterus.mtx.tsv/Normal_Uterus/genes.tsv"),
    ]
    rows = []
    for group, donor_id, matrix, genes in uterus_specs:
        rows.append(score_matrix_market(ROOT / matrix, ROOT / genes, group, donor_id, args.min_total_counts))
    # Reuse the already computed cell-level scores for lung groups. The
    # conservative LAMCORE subset is used for pulmonary LAM; normal lung uses
    # all normal-lung reference cells.
    lung_lam = h5ad_group_summary(ROOT / "data/processed/GSE135851_lam_states.h5ad", "LAM_lung_core_reanalysis", 'lamcore_candidate == True')
    normal_lung = h5ad_group_summary(ROOT / "data/processed/external/GSE122960_normal_lung.h5ad", "normal_lung_reference")
    # h5ad_group_summary reads the combined adaptation score table, so the
    # first argument is retained only for provenance and does not reload data.
    donor_scores = pd.concat([pd.DataFrame(rows), lung_lam.assign(group="pulmonary_LAM"), normal_lung.assign(group="normal_lung")], ignore_index=True, sort=False)
    donor_scores.to_csv(out / "four_group_donor_program_scores.csv", index=False)
    contrast_rows = []
    for program in PROGRAMS:
        col = f"{program}_median"
        group_values = donor_scores.groupby("group", observed=True)[col].median()
        pulmonary = group_values.get("pulmonary_LAM", np.nan)
        normal_lung_value = group_values.get("normal_lung", np.nan)
        lam_uterus = group_values.get("LAM_uterus", np.nan)
        normal_uterus = group_values.get("normal_uterus", np.nan)
        pulmonary_effect = pulmonary - normal_lung_value if np.isfinite(pulmonary) and np.isfinite(normal_lung_value) else np.nan
        uterine_effect = lam_uterus - normal_uterus if np.isfinite(lam_uterus) and np.isfinite(normal_uterus) else np.nan
        interaction = pulmonary_effect - uterine_effect if np.isfinite(pulmonary_effect) and np.isfinite(uterine_effect) else np.nan
        if np.isfinite(interaction) and np.isfinite(pulmonary_effect) and np.isfinite(uterine_effect):
            if abs(uterine_effect) < 0.1 and pulmonary_effect > 0.25:
                classification = "lung_acquired_candidate"
            elif pulmonary_effect > 0.25 and uterine_effect > 0.25:
                classification = "lineage_or_LAM_transformation_candidate"
            else:
                classification = "unresolved_or_reference_sensitive"
        else:
            classification = "insufficient"
        contrast_rows.append({
            "program": program,
            "pulmonary_LAM_minus_normal_lung": pulmonary_effect,
            "LAM_uterus_minus_normal_uterus": uterine_effect,
            "lung_acquired_interaction": interaction,
            "classification": classification,
            "interpretation": "population-level descriptive interaction; not patient-matched causal evidence",
        })
    contrasts = pd.DataFrame(contrast_rows)
    contrasts.to_csv(out / "four_group_lung_acquired_interactions.csv", index=False)
    manifest = {
        "groups": ["normal_uterus", "LAM_uterus", "normal_lung", "pulmonary_LAM"],
        "uterus_inputs": [{"group": x[0], "donor_id": x[1], "matrix": x[2], "genes": x[3]} for x in uterus_specs],
        "lung_inputs": ["results/adaptation/adaptation_cell_scores.csv"],
        "min_total_counts_for_uterus_matrix_market": args.min_total_counts,
        "contrasts": {
            "LAM_transformation": "LAM_uterus - normal_uterus",
            "pulmonary_disease_effect": "pulmonary_LAM - normal_lung",
            "lung_acquired_interaction": "(pulmonary_LAM - normal_lung) - (LAM_uterus - normal_uterus)",
        },
        "limitations": ["uterus references currently have one source specimen each", "groups are not patient matched", "mouse uterus remains a separate auxiliary lineage reference"],
        "interpretation": "Do not call a tissue difference lung-acquired unless the four-group interaction supports it.",
    }
    (out / "four_group_adaptation_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(json.dumps({"donor_rows": len(donor_scores), "program_rows": len(contrasts), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
