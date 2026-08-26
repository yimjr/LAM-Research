"""Validate the stable LINCS modules against an external RET-oriented profile.

GSE49414 compares DMSO with 24-hour RPI-1 treatment in the human TPC1 cell
line.  RPI-1 is a pharmacological RET-oriented reference, not a RET-selective
genetic perturbation, so this script reports module reproduction and gene-level
directional evidence without calling it single-target validation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests

from analyze_external import collapse_duplicate_genes, map_probe_symbols, read_series_matrix
from common import CANDIDATE_PROGRAMS, CANDIDATE_VALIDATION, ROOT, sha256_file, write_json


DATASET = "GSE49414"
MATRIX_DEFAULT = ROOT / "data/raw/GSE49414/GSE49414_series_matrix.txt.gz"
ANNOTATION_DEFAULT = ROOT / "data/raw/GPL17518/GPL17518_family.soft.gz"


def load_modules() -> dict[str, dict[str, object]]:
    common = set(
        pd.read_csv(
            CANDIDATE_PROGRAMS / "LINCS_cross_release_common_gene_comparison.csv.gz",
            usecols=["gene"],
        )["gene"].astype(str).str.upper()
    )
    matches = pd.read_csv(CANDIDATE_PROGRAMS / "LINCS_cross_release_module_matches.csv")
    matches = matches.loc[matches["same_module_first_pass"].eq(True)].reset_index(drop=True)
    modules: dict[str, dict[str, object]] = {}
    for idx, row in matches.iterrows():
        module_id = f"{row['scope']}__stable_{idx + 1}"
        genes = {
            gene.strip().upper()
            for gene in str(row["common_genes"]).split(";")
            if gene and gene != "nan"
        } & common
        modules[module_id] = {
            "module_id": module_id,
            "scope": str(row["scope"]),
            "genes": sorted(genes),
            "jaccard": float(row["jaccard"]),
        }
    return modules


def load_disease_panel() -> pd.DataFrame:
    panel = pd.read_csv(ROOT / "results/signatures/GSE179044_cmap_query_signatures.csv")
    panel = panel.loc[
        panel["contrast"].eq("tsc2_loss_plastic")
        & panel["default_cmap_query"].astype(bool),
        ["gene", "direction", "signed_score"],
    ].copy()
    panel["gene"] = panel["gene"].astype(str).str.upper()
    panel = panel.drop_duplicates("gene").set_index("gene")
    panel = panel.rename(columns={"signed_score": "disease_signed_score", "direction": "disease_direction"})
    return panel


def parse_external_profile(matrix_path: Path, annotation_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    table = read_series_matrix(matrix_path)
    if table.shape[1] != 6:
        raise ValueError(f"Expected six GSE49414 samples, found {table.shape[1]}")
    table = map_probe_symbols(table, annotation_path)
    table = collapse_duplicate_genes(table)
    table.columns = [
        "DMSO_24h_1",
        "DMSO_24h_2",
        "DMSO_24h_3",
        "RPI_24h_1",
        "RPI_24h_2",
        "RPI_24h_3",
    ]
    table = table.apply(pd.to_numeric, errors="coerce")
    effects = pd.DataFrame(
        {
            "rpi1_minus_dmso": table[["RPI_24h_1", "RPI_24h_2", "RPI_24h_3"]].mean(axis=1)
            - table[["DMSO_24h_1", "DMSO_24h_2", "DMSO_24h_3"]].mean(axis=1),
        }
    )
    p_values = []
    for gene in table.index:
        dmso = table.loc[gene, ["DMSO_24h_1", "DMSO_24h_2", "DMSO_24h_3"]].to_numpy(float)
        rpi = table.loc[gene, ["RPI_24h_1", "RPI_24h_2", "RPI_24h_3"]].to_numpy(float)
        p_values.append(ttest_ind(rpi, dmso, equal_var=False, nan_policy="omit").pvalue)
    effects["p_value_welch"] = p_values
    effects["q_value_bh"] = multipletests(effects["p_value_welch"].fillna(1.0), method="fdr_bh")[1]
    return table, effects


def score_modules(effects: pd.DataFrame, panel: pd.DataFrame, modules: dict[str, dict[str, object]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_rows = []
    module_rows = []
    for module_id, module in modules.items():
        genes = list(module["genes"])
        present = [gene for gene in genes if gene in effects.index and gene in panel.index]
        for gene in genes:
            disease_score = float(panel.loc[gene, "disease_signed_score"]) if gene in panel.index else np.nan
            drug_effect = float(effects.loc[gene, "rpi1_minus_dmso"]) if gene in effects.index else np.nan
            if not np.isfinite(disease_score) or not np.isfinite(drug_effect):
                direction = "not_available"
            elif disease_score * drug_effect < 0:
                direction = "reversal"
            elif disease_score * drug_effect > 0:
                direction = "mimic"
            else:
                direction = "neutral"
            expected = "reversal" if module["scope"] == "reversal_only" else "mimic"
            gene_rows.append(
                {
                    "dataset": DATASET,
                    "reference": "RPI-1",
                    "axis": "RET-oriented, non-selective proof",
                    "module_id": module_id,
                    "module_scope": module["scope"],
                    "gene": gene,
                    "disease_signed_score": disease_score,
                    "external_drug_effect": drug_effect,
                    "external_p_value": float(effects.loc[gene, "p_value_welch"]) if gene in effects.index else np.nan,
                    "external_q_value": float(effects.loc[gene, "q_value_bh"]) if gene in effects.index else np.nan,
                    "direction": direction,
                    "expected_module_direction": expected,
                    "matches_expected_module_direction": direction == expected,
                    "weighted_contribution_auxiliary": disease_score * drug_effect if direction != "not_available" else np.nan,
                }
            )
        module_effects = effects.reindex(present)
        module_panel = panel.reindex(present)
        products = module_effects["rpi1_minus_dmso"].to_numpy(float) * module_panel["disease_signed_score"].to_numpy(float)
        expected = products < 0 if module["scope"] == "reversal_only" else products > 0
        direction = np.where(products < 0, "reversal", np.where(products > 0, "mimic", "neutral"))
        module_rows.append(
            {
                "dataset": DATASET,
                "reference": "RPI-1",
                "axis": "RET-oriented, non-selective proof",
                "module_id": module_id,
                "module_scope": module["scope"],
                "module_jaccard_cross_release": float(module["jaccard"]),
                "n_module_genes": len(genes),
                "n_genes_analyzed": len(present),
                "n_genes_not_available": len(genes) - len(present),
                "external_reversal_fraction": float(np.mean(direction == "reversal")) if present else np.nan,
                "external_mimic_fraction": float(np.mean(direction == "mimic")) if present else np.nan,
                "external_neutral_fraction": float(np.mean(direction == "neutral")) if present else np.nan,
                "expected_module_direction_fraction": float(np.mean(expected)) if present else np.nan,
                "n_expected_direction": int(np.sum(expected)),
                "n_significant_genes_bh05": int(np.sum(module_effects["q_value_bh"].to_numpy(float) < 0.05)),
                "interpretation": "external pharmacological module reference; RPI-1 is not a RET-selective genetic perturbation",
            }
        )
    return pd.DataFrame(gene_rows), pd.DataFrame(module_rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=MATRIX_DEFAULT)
    parser.add_argument("--annotation", type=Path, default=ANNOTATION_DEFAULT)
    args = parser.parse_args()
    matrix_path = args.input if args.input.is_absolute() else ROOT / args.input
    annotation_path = args.annotation if args.annotation.is_absolute() else ROOT / args.annotation
    table, effects = parse_external_profile(matrix_path, annotation_path)
    panel = load_disease_panel()
    modules = load_modules()
    genes, modules_result = score_modules(effects, panel, modules)
    out = CANDIDATE_VALIDATION
    out.mkdir(parents=True, exist_ok=True)
    genes.to_csv(out / "GSE49414_RET_RPI1_module_gene_validation.csv.gz", index=False, compression="gzip")
    modules_result.to_csv(out / "GSE49414_RET_RPI1_module_validation.csv", index=False)
    write_json(ROOT / "manifests" / "GSE49414_RET_reference.json", {
        "dataset": DATASET,
        "role": "independent RET-oriented pharmacological reference",
        "treatment": "RPI-1, 24 hours",
        "control": "DMSO, 24 hours",
        "n_replicates_per_condition": 3,
        "platform": "GPL17518",
        "matrix_sha256": sha256_file(matrix_path),
        "annotation_sha256": sha256_file(annotation_path),
        "interpretation": "RPI-1 is an external pharmacological reference, not a clean RET-only validation; module-level directional evidence is reported without assigning causality to RET.",
        "n_genes_expression": int(table.shape[0]),
        "n_genes_panel": int(len(panel)),
        "n_gene_rows": int(len(genes)),
    })
    print(modules_result.to_string(index=False))


if __name__ == "__main__":
    main()
