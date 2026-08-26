#!/usr/bin/env python3
"""Score known LAM programs in public TSC2/rapamycin perturbation data."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def load_programs(path: Path) -> dict[str, list[str]]:
    payload = yaml.safe_load(path.read_text()) or {}
    programs = {}
    for item in payload.get("programs", []):
        name = str(item.get("program_name", ""))
        genes = [str(g).upper() for g in item.get("genes", [])]
        if name and genes:
            programs[name] = genes
    return programs


def gse179044_metadata() -> pd.DataFrame:
    rows = []
    labels = [
        ("S01", "WT", "hydrogel", "vehicle", 1), ("S02", "TSC2ko", "hydrogel", "vehicle", 1),
        ("S03", "WT", "hydrogel", "rapamycin", 1), ("S04", "TSC2ko", "hydrogel", "rapamycin", 1),
        ("S05", "WT", "plastic", "vehicle", 1), ("S06", "TSC2ko", "plastic", "vehicle", 1),
        ("S07", "WT", "plastic", "rapamycin", 1), ("S08", "TSC2ko", "plastic", "rapamycin", 1),
        ("S09", "WT", "hydrogel", "vehicle", 2), ("S10", "TSC2ko", "hydrogel", "vehicle", 2),
        ("S11", "WT", "hydrogel", "rapamycin", 2), ("S12", "TSC2ko", "hydrogel", "rapamycin", 2),
        ("S13", "WT", "plastic", "vehicle", 2), ("S14", "TSC2ko", "plastic", "vehicle", 2),
        ("S15", "WT", "plastic", "rapamycin", 2), ("S16", "TSC2ko", "plastic", "rapamycin", 2),
    ]
    for sample, genotype, environment, treatment, replicate in labels:
        rows.append({"sample_id": f"Pietrobon_{sample}", "genotype": genotype, "environment": environment, "treatment": treatment, "replicate": replicate})
    return pd.DataFrame(rows)


def score_gene_programs(matrix: pd.DataFrame, programs: dict[str, list[str]], sample_meta: pd.DataFrame) -> pd.DataFrame:
    matrix.index = matrix.index.astype(str).str.upper()
    rows = []
    for program, genes in programs.items():
        available = [gene for gene in genes if gene in matrix.index]
        if not available:
            continue
        for sample_id in matrix.columns:
            rows.append({
                "program": program,
                "sample_id": sample_id,
                "available_genes": len(available),
                "score": float(matrix.loc[available, sample_id].mean()),
            })
    result = pd.DataFrame(rows)
    return result.merge(sample_meta, on="sample_id", how="left")


def contrast_table(scores: pd.DataFrame, strata: list[str], treatment_a: str, treatment_b: str, label: str) -> pd.DataFrame:
    rows = []
    for keys, group in scores.groupby(["program", *strata], observed=True):
        if not isinstance(keys, tuple):
            keys = (keys,)
        values_a = group.loc[group["treatment"].eq(treatment_a), "score"]
        values_b = group.loc[group["treatment"].eq(treatment_b), "score"]
        if len(values_a) == 0 or len(values_b) == 0:
            continue
        row = {"contrast": label, **dict(zip(["program", *strata], keys)), "condition_a": treatment_a, "condition_b": treatment_b, "n_a": len(values_a), "n_b": len(values_b), "mean_a": float(values_a.mean()), "mean_b": float(values_b.mean())}
        row["effect_a_minus_b"] = row["mean_a"] - row["mean_b"]
        rows.append(row)
    return pd.DataFrame(rows)


def analyze_gse179044(programs: dict[str, list[str]], out: Path) -> None:
    path = ROOT / "data/raw/perturbation/GSE179044/GSE179044_rlog_counts_matrix.csv.gz"
    table = pd.read_csv(path)
    sample_cols = [c for c in table.columns if c.startswith("Pietrobon_S")]
    matrix = table.set_index("gene_name")[sample_cols].apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()
    meta = gse179044_metadata()
    scores = score_gene_programs(matrix, programs, meta)
    contrast_parts = [
        contrast_table(scores, ["genotype", "environment"], "rapamycin", "vehicle", "rapamycin_within_genotype_environment"),
    ]
    # Genotype and environment contrasts use the same score table but are
    # explicitly labeled so they cannot be confused with treatment effects.
    for program in sorted(scores["program"].unique()):
        subset = scores[scores["program"].eq(program)]
        for environment in sorted(subset["environment"].unique()):
            for treatment in sorted(subset["treatment"].unique()):
                g = subset[(subset["environment"] == environment) & (subset["treatment"] == treatment)]
                a = g.loc[g["genotype"].eq("TSC2ko"), "score"]
                b = g.loc[g["genotype"].eq("WT"), "score"]
                if len(a) and len(b):
                    contrast_parts.append(pd.DataFrame([{"contrast": "TSC2ko_minus_WT", "program": program, "environment": environment, "treatment": treatment, "n_a": len(a), "n_b": len(b), "mean_a": float(a.mean()), "mean_b": float(b.mean()), "effect_a_minus_b": float(a.mean() - b.mean())}]))
        for genotype in sorted(subset["genotype"].unique()):
            for treatment in sorted(subset["treatment"].unique()):
                g = subset[(subset["genotype"] == genotype) & (subset["treatment"] == treatment)]
                a = g.loc[g["environment"].eq("hydrogel"), "score"]
                b = g.loc[g["environment"].eq("plastic"), "score"]
                if len(a) and len(b):
                    contrast_parts.append(pd.DataFrame([{"contrast": "hydrogel_minus_plastic", "program": program, "genotype": genotype, "treatment": treatment, "n_a": len(a), "n_b": len(b), "mean_a": float(a.mean()), "mean_b": float(b.mean()), "effect_a_minus_b": float(a.mean() - b.mean())}]))
    contrasts = pd.concat([x for x in contrast_parts if not x.empty], ignore_index=True)
    meta.to_csv(out / "GSE179044_sample_metadata.csv", index=False)
    scores.to_csv(out / "GSE179044_program_scores.csv", index=False)
    contrasts.to_csv(out / "GSE179044_program_contrasts.csv", index=False)
    (out / "GSE179044_analysis_manifest.json").write_text(json.dumps({"source": "GSE179044", "matrix": str(path.relative_to(ROOT)), "matrix_type": "rlog", "design": "TSC2ko/WT x hydrogel/plastic x vehicle/rapamycin x 2 replicates", "interpretation": "Perturbation evidence for candidate persistence mechanisms; not patient-level treatment evidence."}, indent=2, ensure_ascii=False))


def analyze_gse84476(programs: dict[str, list[str]], out: Path) -> None:
    # GEO labels this as tabular text, but the public file is whitespace
    # separated; use a regex separator to preserve all sample columns.
    abundance = pd.read_csv(ROOT / "data/raw/perturbation/GSE84476/GSE84476_VKMJ-STAT3_kallisto_abundance.txt.gz", sep=r"\s+", engine="python")
    fold = pd.read_csv(ROOT / "data/raw/perturbation/GSE84476/GSE84476_VKMJ-STAT3_fold_changes.txt.gz", sep="\t")
    transcript_gene = fold.dropna(subset=["gene"]).drop_duplicates("target_id").set_index("target_id")["gene"].astype(str).str.upper().to_dict()
    tpm_cols = [c for c in abundance.columns if c.startswith("tpm.")]
    sample_names = [c.removeprefix("tpm.") for c in tpm_cols]
    matrix = abundance.set_index("target_id")[tpm_cols].apply(pd.to_numeric, errors="coerce")
    matrix.columns = sample_names
    matrix["gene"] = [transcript_gene.get(str(t), "") for t in matrix.index]
    matrix = matrix[matrix["gene"].ne("")].drop(columns=["gene"])
    matrix = matrix.groupby([transcript_gene.get(str(t), "") for t in matrix.index]).mean()
    metadata = []
    for sample in sample_names:
        cell_line = sample.split("_")[0]
        treatment = sample.split("_", 1)[1]
        metadata.append({"sample_id": sample, "cell_line": cell_line, "genotype": "TSC2-null" if cell_line == "102cell" else "TSC2-reexpressed", "treatment": treatment})
    meta = pd.DataFrame(metadata)
    scores = score_gene_programs(matrix, programs, meta)
    parts = [contrast_table(scores, ["genotype"], "Rap", "siCtrl", "rapamycin_within_genotype"), contrast_table(scores, ["genotype"], "siSTAT3", "siCtrl", "STAT3_perturbation_within_genotype")]
    for program in sorted(scores["program"].unique()):
        group = scores[(scores["program"] == program) & (scores["treatment"] == "siCtrl")]
        a = group.loc[group["genotype"].eq("TSC2-null"), "score"]
        b = group.loc[group["genotype"].eq("TSC2-reexpressed"), "score"]
        if len(a) and len(b):
            parts.append(pd.DataFrame([{"contrast": "TSC2_null_minus_reexpressed", "program": program, "n_a": len(a), "n_b": len(b), "mean_a": float(a.mean()), "mean_b": float(b.mean()), "effect_a_minus_b": float(a.mean() - b.mean())}]))
    contrasts = pd.concat([x for x in parts if not x.empty], ignore_index=True)
    meta.to_csv(out / "GSE84476_sample_metadata.csv", index=False)
    scores.to_csv(out / "GSE84476_program_scores.csv", index=False)
    contrasts.to_csv(out / "GSE84476_program_contrasts.csv", index=False)
    (out / "GSE84476_analysis_manifest.json").write_text(json.dumps({"source": "GSE84476", "matrix": "kallisto TPM", "transcript_to_gene": "mapped from GEO fold-change table", "genotype_mapping": "102cell=TSC2-null; 103cell=TSC2-reexpressed; confirmed from GEO family metadata", "interpretation": "Perturbation evidence for candidate TSC2-loss and rapamycin-associated mechanisms; one LAM-derived cell-line system, not patient-level evidence."}, indent=2, ensure_ascii=False))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/perturbation")
    parser.add_argument("--programs", default="config/known_lam_programs.yaml")
    args = parser.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    programs = load_programs(ROOT / args.programs)
    analyze_gse179044(programs, out)
    analyze_gse84476(programs, out)
    print(json.dumps({"datasets": ["GSE179044", "GSE84476"], "programs": len(programs), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
