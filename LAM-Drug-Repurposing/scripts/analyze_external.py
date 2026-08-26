"""Analyze the external GSE27982 2x2 and create a comparable response table."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from common import ROOT, classify_ratio, collapse_duplicate_genes, signed_ratio, write_json


def read_series_matrix(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", comment="!", index_col=0, compression="infer")


def map_probe_symbols(table: pd.DataFrame, annotation_path: Path | None) -> pd.DataFrame:
    if annotation_path is None or not annotation_path.exists():
        return table
    # GPL339 annotation files contain a metadata preamble and a table marker;
    # start parsing at the actual tabular header, then ignore the end marker.
    opener = gzip.open if str(annotation_path).endswith(".gz") else open
    with opener(annotation_path, "rt", encoding="utf-8", errors="replace") as handle:
        lines = handle.readlines()
    start = next((i + 1 for i, line in enumerate(lines) if line.strip() == "!platform_table_begin"), None)
    if start is None:
        return table
    end = next(
        (i for i, line in enumerate(lines[start:], start=start) if line.strip() == "!platform_table_end"),
        len(lines),
    )
    annotation = pd.read_csv(
        annotation_path,
        sep="\t",
        skiprows=start,
        nrows=max(0, end - start - 1),
        comment="!",
        compression="infer",
        dtype=str,
    )
    id_col = next((column for column in annotation.columns if str(column).strip('"') in {"ID", "ID_REF"}), annotation.columns[0])
    symbol_col = next(
        (
            column
            for column in annotation.columns
            if "gene symbol" in str(column).lower()
            or str(column).lower() in {"gene_symbol", "symbol"}
        ),
        None,
    )
    if symbol_col is None:
        return table
    lookup = annotation.set_index(id_col)[symbol_col].astype(str).str.strip().replace({"nan": ""})
    mapped = table.copy()
    mapped.index = [lookup.get(str(probe), str(probe)) or str(probe) for probe in mapped.index]
    mapped = mapped.loc[~mapped.index.str.contains("///|//", regex=True)]
    return mapped


def analyze(table: pd.DataFrame, minimum_effect: float = 0.5) -> pd.DataFrame:
    samples = pd.DataFrame([
        {"sample_id": "GSM692432", "genotype": "WT", "treatment": "vehicle"},
        {"sample_id": "GSM692433", "genotype": "WT", "treatment": "vehicle"},
        {"sample_id": "GSM692434", "genotype": "WT", "treatment": "rapamycin"},
        {"sample_id": "GSM692435", "genotype": "WT", "treatment": "rapamycin"},
        {"sample_id": "GSM692436", "genotype": "KO", "treatment": "vehicle"},
        {"sample_id": "GSM692437", "genotype": "KO", "treatment": "vehicle"},
        {"sample_id": "GSM692438", "genotype": "KO", "treatment": "rapamycin"},
        {"sample_id": "GSM692439", "genotype": "KO", "treatment": "rapamycin"},
    ])
    if table.shape[1] != len(samples):
        raise ValueError(f"Expected 8 GSE27982 samples, found {table.shape[1]}")
    # GEO series matrices are already normalized/log-scaled expression values.
    # Do not apply a second log transform; only coalesce probes after mapping.
    table = table.apply(pd.to_numeric, errors="coerce")
    table = collapse_duplicate_genes(table)
    table.columns = samples.sample_id
    def mean(genotype: str, treatment: str) -> pd.Series:
        selected = samples.loc[
            samples.genotype.eq(genotype) & samples.treatment.eq(treatment),
            "sample_id",
        ].tolist()
        return table.loc[:, selected].mean(axis=1)
    d0 = mean("KO", "vehicle") - mean("WT", "vehicle")
    d1 = mean("KO", "rapamycin") - mean("WT", "rapamycin")
    ko_response = mean("KO", "rapamycin") - mean("KO", "vehicle")
    wt_response = mean("WT", "rapamycin") - mean("WT", "vehicle")
    result = pd.DataFrame({
        "tsc2_loss": d0,
        "residual_after_rapamycin": d1,
        "genotype_dependent_rapamycin_response": ko_response - wt_response,
    })
    result["signed_residual_ratio"] = signed_ratio(result["residual_after_rapamycin"], result["tsc2_loss"], minimum_effect)
    result["absolute_residual_ratio"] = result["signed_residual_ratio"].abs()
    result["residual_class"] = classify_ratio(result["signed_residual_ratio"])
    result.index.name = "probe_or_gene"
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/raw/GSE27982/GSE27982_series_matrix.txt.gz")
    parser.add_argument("--annotation", default="data/raw/GPL339/GPL339.annot.gz")
    args = parser.parse_args()
    path = ROOT / args.input
    if not path.exists():
        raise FileNotFoundError(path)
    table = read_series_matrix(path)
    table = map_probe_symbols(table, ROOT / args.annotation)
    analysis_config = yaml.safe_load((ROOT / "config" / "analysis.yaml").read_text())
    result = analyze(table, float(analysis_config["min_baseline_effect"]))
    out = ROOT / "results" / "tables" / "GSE27982_external_response.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(out)
    write_json(ROOT / "manifests" / "GSE27982_analysis.json", {
        "dataset": "GSE27982",
        "role": "external Tsc2 x rapamycin response validation",
        "interpretation": "Genotype-dependent rapamycin response is not automatically classified as escape under low serum.",
        "ratio_gate": float(analysis_config["min_baseline_effect"]),
        "n_features": int(result.shape[0]),
    })
    print({"dataset": "GSE27982", "n_features": int(result.shape[0])})


if __name__ == "__main__":
    main()
