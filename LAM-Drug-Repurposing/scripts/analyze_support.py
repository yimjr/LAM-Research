"""Analyze the non-factorial GSE16944 historical LAM-like support dataset.

GSE16944 does not contain TSC2-restored + rapamycin, so this script reports
descriptive comparisons only. It deliberately does not label any contrast as
formal residual or GxR evidence.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from analyze_external import map_probe_symbols, read_series_matrix
from common import ROOT, collapse_duplicate_genes, write_json


SAMPLES = pd.DataFrame([
    {"sample_id": "GSM424642", "genotype": "KO", "treatment": "rapamycin"},
    {"sample_id": "GSM426549", "genotype": "KO", "treatment": "rapamycin"},
    {"sample_id": "GSM426550", "genotype": "KO", "treatment": "rapamycin"},
    {"sample_id": "GSM426551", "genotype": "KO", "treatment": "vehicle"},
    {"sample_id": "GSM426552", "genotype": "KO", "treatment": "vehicle"},
    {"sample_id": "GSM426553", "genotype": "KO", "treatment": "vehicle"},
    {"sample_id": "GSM426554", "genotype": "WT_addback", "treatment": "vehicle"},
    {"sample_id": "GSM426555", "genotype": "WT_addback", "treatment": "vehicle"},
    {"sample_id": "GSM426556", "genotype": "WT_addback", "treatment": "vehicle"},
])


def group_mean(table: pd.DataFrame, genotype: str, treatment: str) -> pd.Series:
    ids = SAMPLES.loc[
        SAMPLES.genotype.eq(genotype) & SAMPLES.treatment.eq(treatment),
        "sample_id",
    ]
    return table.loc[:, list(ids)].mean(axis=1)


def main() -> None:
    input_path = ROOT / "data/raw/GSE16944/GSE16944_series_matrix.txt.gz"
    annotation_path = ROOT / "data/raw/GPL2895/GPL2895.annot.gz"
    table = map_probe_symbols(read_series_matrix(input_path), annotation_path)
    table = collapse_duplicate_genes(table)
    table = table.loc[:, SAMPLES.sample_id]

    ko_vehicle = group_mean(table, "KO", "vehicle")
    ko_rapa = group_mean(table, "KO", "rapamycin")
    wt_vehicle = group_mean(table, "WT_addback", "vehicle")
    result = pd.DataFrame({
        "KO_vehicle_minus_WT_vehicle": ko_vehicle - wt_vehicle,
        "KO_rapamycin_minus_KO_vehicle": ko_rapa - ko_vehicle,
        # Descriptive only: this is not D(rapamycin), because WT+rapamycin is absent.
        "KO_rapamycin_vs_WT_vehicle_descriptive": ko_rapa - wt_vehicle,
    })
    result.index.name = "gene"
    out_dir = ROOT / "results" / "tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_dir / "GSE16944_historical_support.csv")

    config = yaml.safe_load((ROOT / "config" / "analysis.yaml").read_text())
    module_rows = []
    for module, genes in config["module_sets"].items():
        available = [gene for gene in genes if gene in result.index]
        row = {"module": module, "n_genes": len(available)}
        for column in result.columns:
            row[column] = float(result.loc[available, column].mean()) if available else np.nan
        module_rows.append(row)
    pd.DataFrame(module_rows).to_csv(out_dir / "GSE16944_module_support.csv", index=False)

    key_genes = [gene for gene in ("MMP2", "MMP9", "MMP11", "COL1A1", "COL3A1", "FN1", "SPARC") if gene in result.index]
    key_values = result.loc[key_genes].to_dict(orient="index")
    write_json(ROOT / "manifests" / "GSE16944_analysis.json", {
        "dataset": "GSE16944",
        "role": "historical rapamycin-insensitive LAM-like support",
        "n_features": int(result.shape[0]),
        "n_samples": int(table.shape[1]),
        "formal_limit": "No TSC2-restored + rapamycin group; no formal D(rapamycin) or complete GxR interaction.",
        "key_genes": key_values,
    })
    print({"dataset": "GSE16944", "n_features": int(result.shape[0]), "key_genes": key_values})


if __name__ == "__main__":
    main()
