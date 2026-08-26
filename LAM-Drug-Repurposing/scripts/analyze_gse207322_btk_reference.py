"""Validate stable modules against an external BTK/ibrutinib experiment.

GSE207322 profiles TMD8 cells with WT BTK or clinically observed C481
mutations, with DMSO or 10 nM ibrutinib.  The experiment is useful because
BTK-resistant alleles provide a pharmacological specificity contrast, but
ibrutinib remains a multi-target drug and this is not a LAM model.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

from common import CANDIDATE_PROGRAMS, CANDIDATE_VALIDATION, ROOT, sha256_file, write_json


DATASET = "GSE207322"
MATRIX_DEFAULT = ROOT / "data/raw/GSE207322/GSE207322_TPM.txt.gz"


def load_modules() -> dict[str, dict[str, object]]:
    common = set(
        pd.read_csv(
            CANDIDATE_PROGRAMS / "LINCS_cross_release_common_gene_comparison.csv.gz",
            usecols=["gene"],
        )["gene"].astype(str).str.upper()
    )
    matches = pd.read_csv(CANDIDATE_PROGRAMS / "LINCS_cross_release_module_matches.csv")
    matches = matches.loc[matches["same_module_first_pass"].eq(True)].reset_index(drop=True)
    modules = {}
    for idx, row in matches.iterrows():
        module_id = f"{row['scope']}__stable_{idx + 1}"
        modules[module_id] = {
            "module_id": module_id,
            "scope": str(row["scope"]),
            "genes": sorted(
                {
                    gene.strip().upper()
                    for gene in str(row["common_genes"]).split(";")
                    if gene and gene != "nan"
                }
                & common
            ),
            "jaccard": float(row["jaccard"]),
        }
    return modules


def load_disease_panel() -> pd.DataFrame:
    panel = pd.read_csv(ROOT / "results/signatures/GSE179044_cmap_query_signatures.csv")
    panel = panel.loc[
        panel["contrast"].eq("tsc2_loss_plastic")
        & panel["default_cmap_query"].astype(bool),
        ["gene", "signed_score"],
    ].copy()
    panel["gene"] = panel["gene"].astype(str).str.upper()
    return panel.drop_duplicates("gene").set_index("gene").rename(columns={"signed_score": "disease_signed_score"})


def read_tpm(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    raw = pd.read_csv(path, sep="\t", compression="infer")
    gene_column = raw.columns[0]
    sample_columns = list(raw.columns[1:])
    metadata = []
    for column in sample_columns:
        match = re.fullmatch(r"TMD8_(parental|C481[SFYR])_(DMSO|ibrutinib)_R([12]) TPM", str(column))
        if not match:
            raise ValueError(f"Unrecognized GSE207322 sample column: {column}")
        genotype, treatment, replicate = match.groups()
        metadata.append(
            {
                "sample_id": column,
                "genotype": "WT" if genotype == "parental" else genotype,
                "treatment": treatment,
                "replicate": int(replicate),
            }
        )
    metadata = pd.DataFrame(metadata)
    genes = raw[gene_column].astype(str).str.rsplit("_", n=1).str[-1].str.upper()
    table = raw.drop(columns=[gene_column]).apply(pd.to_numeric, errors="coerce")
    table.index = genes
    table = table.groupby(level=0, sort=False).mean()
    # TPM is not a log-scale measure; log1p makes between-gene effect sizes
    # more comparable while preserving the direction of the perturbation.
    table = np.log1p(table)
    return table, metadata


def score_modules(table: pd.DataFrame, metadata: pd.DataFrame, panel: pd.DataFrame, modules: dict[str, dict[str, object]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_rows = []
    module_rows = []
    for module_id, module in modules.items():
        module_genes = list(module["genes"])
        present = [gene for gene in module_genes if gene in table.index and gene in panel.index]
        for genotype in ["WT", "C481S", "C481F", "C481Y", "C481R"]:
            selected = metadata.loc[metadata["genotype"].eq(genotype)]
            dmso_columns = selected.loc[selected["treatment"].eq("DMSO"), "sample_id"].tolist()
            drug_columns = selected.loc[selected["treatment"].eq("ibrutinib"), "sample_id"].tolist()
            effects = table.loc[present, drug_columns].mean(axis=1) - table.loc[present, dmso_columns].mean(axis=1)
            disease = panel.loc[present, "disease_signed_score"]
            products = effects * disease
            directions = np.where(products < 0, "reversal", np.where(products > 0, "mimic", "neutral"))
            for gene in module_genes:
                effect = float(effects.loc[gene]) if gene in effects.index else np.nan
                disease_score = float(panel.loc[gene, "disease_signed_score"]) if gene in panel.index else np.nan
                if not np.isfinite(effect) or not np.isfinite(disease_score):
                    direction = "not_available"
                elif effect * disease_score < 0:
                    direction = "reversal"
                elif effect * disease_score > 0:
                    direction = "mimic"
                else:
                    direction = "neutral"
                gene_rows.append(
                    {
                        "dataset": DATASET,
                        "reference": "ibrutinib",
                        "axis": "BTK-oriented, C481 resistance contrast",
                        "module_id": module_id,
                        "module_scope": module["scope"],
                        "genotype": genotype,
                        "gene": gene,
                        "disease_signed_score": disease_score,
                        "ibrutinib_minus_dmso_log1p_tpm": effect,
                        "direction": direction,
                    }
                )
            expected = products < 0
            module_rows.append(
                {
                    "dataset": DATASET,
                    "reference": "ibrutinib",
                    "axis": "BTK-oriented, C481 resistance contrast",
                    "module_id": module_id,
                    "module_scope": module["scope"],
                    "genotype": genotype,
                    "module_jaccard_cross_release": float(module["jaccard"]),
                    "n_module_genes": len(module_genes),
                    "n_genes_analyzed": len(present),
                    "n_genes_not_available": len(module_genes) - len(present),
                    "reversal_fraction": float(np.mean(directions == "reversal")) if present else np.nan,
                    "mimic_fraction": float(np.mean(directions == "mimic")) if present else np.nan,
                    "neutral_fraction": float(np.mean(directions == "neutral")) if present else np.nan,
                    "btk_expected_reversal_fraction": float(np.mean(expected)) if present else np.nan,
                    "mean_effect_log1p_tpm": float(effects.mean()) if present else np.nan,
                    "interpretation": "external BTK-oriented pharmacological reference; C481 mutants are a resistance contrast, not a LAM validation",
                }
            )
    module_result = pd.DataFrame(module_rows)
    return pd.DataFrame(gene_rows), module_result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=MATRIX_DEFAULT)
    args = parser.parse_args()
    matrix_path = args.input if args.input.is_absolute() else ROOT / args.input
    table, metadata = read_tpm(matrix_path)
    panel = load_disease_panel()
    modules = load_modules()
    genes, modules_result = score_modules(table, metadata, panel, modules)
    out = CANDIDATE_VALIDATION
    out.mkdir(parents=True, exist_ok=True)
    genes.to_csv(out / "GSE207322_BTK_ibrutinib_module_gene_validation.csv.gz", index=False, compression="gzip")
    modules_result.to_csv(out / "GSE207322_BTK_ibrutinib_module_validation.csv", index=False)
    write_json(ROOT / "manifests" / "GSE207322_BTK_reference.json", {
        "dataset": DATASET,
        "role": "independent BTK-oriented pharmacological reference",
        "treatment": "ibrutinib, 10 nM, 24 hours",
        "genotypes": ["WT", "C481S", "C481F", "C481Y", "C481R"],
        "n_replicates_per_genotype_treatment": 2,
        "matrix_sha256": sha256_file(matrix_path),
        "interpretation": "C481 mutants provide a BTK-resistance contrast; ibrutinib is not a single-target drug and TMD8 is not a LAM model.",
        "n_genes_expression": int(table.shape[0]),
        "n_genes_panel": int(len(panel)),
        "n_gene_rows": int(len(genes)),
    })
    print(modules_result.to_string(index=False))


if __name__ == "__main__":
    main()
