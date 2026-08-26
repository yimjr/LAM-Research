#!/usr/bin/env python3
"""Classify TSC2-loss programs by effect-size retention after rapamycin.

The primary result is not a p-value test of "no decrease". For every gene we
compare the TSC2-loss effect before and after rapamycin and report a
suppression fraction plus replicate concordance.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from analyze_perturbation_programs import gse179044_metadata, load_programs

ROOT = Path(__file__).resolve().parents[1]


def read_matrix(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path)
    sample_cols = [c for c in table.columns if c.startswith("Pietrobon_S")]
    table["gene_name"] = table["gene_name"].astype(str).str.upper()
    return table.set_index("gene_name")[sample_cols].apply(pd.to_numeric, errors="coerce").groupby(level=0).mean()


def classify(suppression_fraction: float) -> str:
    if not np.isfinite(suppression_fraction):
        return "insufficient"
    if suppression_fraction < 0:
        return "enhanced_or_reprogrammed"
    if suppression_fraction < 0.25:
        return "fully_retained"
    if suppression_fraction < 0.75:
        return "partially_retained"
    return "suppressed"


def retention_table(matrix: pd.DataFrame, meta: pd.DataFrame, min_effect: float) -> pd.DataFrame:
    rows = []
    for gene in matrix.index.astype(str):
        for environment in ["hydrogel", "plastic"]:
            sub = meta[meta["environment"].eq(environment)]
            per_rep = []
            for replicate in sorted(sub["replicate"].unique()):
                rep = sub[sub["replicate"].eq(replicate)]
                def get(genotype: str, treatment: str) -> float:
                    sample = rep[(rep["genotype"].eq(genotype)) & (rep["treatment"].eq(treatment))]["sample_id"]
                    return float(matrix.loc[gene, sample.iloc[0]]) if len(sample) else np.nan
                wt_vehicle = get("WT", "vehicle")
                ko_vehicle = get("TSC2ko", "vehicle")
                wt_rap = get("WT", "rapamycin")
                ko_rap = get("TSC2ko", "rapamycin")
                e_pre = ko_vehicle - wt_vehicle
                e_post = ko_rap - wt_rap
                sf = 1.0 - (e_post / e_pre) if np.isfinite(e_pre) and abs(e_pre) >= min_effect else np.nan
                per_rep.append({"replicate": int(replicate), "e_pre": e_pre, "e_post": e_post, "suppression_fraction": sf})
            pre = np.asarray([x["e_pre"] for x in per_rep], dtype=float)
            post = np.asarray([x["e_post"] for x in per_rep], dtype=float)
            sf_values = np.asarray([x["suppression_fraction"] for x in per_rep], dtype=float)
            valid_pre = np.isfinite(pre) & (np.abs(pre) >= min_effect)
            pre_direction_concordant = bool(valid_pre.sum() == len(pre) and len(set(np.sign(pre))) == 1)
            post_direction_concordant = bool(np.isfinite(post).all() and (np.all(post == 0) or len(set(np.sign(post[post != 0]))) <= 1))
            mean_pre = float(np.nanmean(pre))
            mean_post = float(np.nanmean(post))
            mean_sf = 1.0 - (mean_post / mean_pre) if np.isfinite(mean_pre) and abs(mean_pre) >= min_effect else np.nan
            rows.append({
                "gene": gene,
                "environment": environment,
                "n_replicates": len(per_rep),
                "replicate_pre_effects": ";".join(f"{x:.6g}" for x in pre),
                "replicate_post_effects": ";".join(f"{x:.6g}" for x in post),
                "mean_pre_effect": mean_pre,
                "mean_post_effect": mean_post,
                "suppression_fraction": mean_sf,
                "category": classify(mean_sf),
                "baseline_effect_size_pass": bool(np.isfinite(mean_pre) and abs(mean_pre) >= min_effect),
                "baseline_direction_concordant": pre_direction_concordant,
                "post_direction_concordant": post_direction_concordant,
                "replicate_retention_eligible": bool(np.isfinite(mean_sf) and pre_direction_concordant),
            })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", default="data/raw/perturbation/GSE179044/GSE179044_rlog_counts_matrix.csv.gz")
    parser.add_argument("--programs", default="config/known_lam_programs.yaml")
    parser.add_argument("--output-dir", default="results/perturbation")
    parser.add_argument("--min-effect", type=float, default=0.5)
    args = parser.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    matrix = read_matrix(ROOT / args.matrix)
    meta = gse179044_metadata()
    genes = retention_table(matrix, meta, args.min_effect)
    genes.to_csv(out / "GSE179044_gene_retention_effects.csv", index=False)
    programs = load_programs(ROOT / args.programs)
    program_rows = []
    for program, program_genes in programs.items():
        subset = genes[genes["gene"].isin(set(program_genes))]
        for environment in ["hydrogel", "plastic"]:
            env = subset[subset["environment"].eq(environment)]
            eligible = env[env["replicate_retention_eligible"]]
            counts = env["category"].value_counts().to_dict()
            program_rows.append({
                "program": program,
                "environment": environment,
                "program_genes_in_matrix": int(len(env)),
                "eligible_genes": int(len(eligible)),
                "mean_suppression_fraction": float(eligible["suppression_fraction"].mean()) if len(eligible) else np.nan,
                "fully_retained_genes": int(counts.get("fully_retained", 0)),
                "partially_retained_genes": int(counts.get("partially_retained", 0)),
                "enhanced_or_reprogrammed_genes": int(counts.get("enhanced_or_reprogrammed", 0)),
                "suppressed_genes": int(counts.get("suppressed", 0)),
                "program_retention_candidate": bool(len(eligible) >= 2 and (eligible["category"].isin(["fully_retained", "partially_retained", "enhanced_or_reprogrammed"]).mean() >= 0.5)),
            })
    program_summary = pd.DataFrame(program_rows)
    program_summary.to_csv(out / "GSE179044_program_retention_effects.csv", index=False)
    cross = genes[genes["replicate_retention_eligible"]].groupby("gene", observed=True).agg(
        environments=("environment", "nunique"),
        mean_suppression_fraction=("suppression_fraction", "mean"),
        categories=("category", lambda x: ";".join(sorted(set(x)))),
    ).reset_index()
    cross["cross_environment_retention_candidate"] = (
        (cross["environments"] == 2)
        & cross["categories"].str.contains("fully_retained|partially_retained|enhanced_or_reprogrammed", regex=True)
    )
    cross.to_csv(out / "GSE179044_cross_environment_retained_genes.csv", index=False)
    manifest = {
        "source": "GSE179044",
        "matrix": args.matrix,
        "definition": "Retention is based on TSC2-loss effect before versus after rapamycin, not p>0.05.",
        "e_pre": "TSC2ko_vehicle - WT_vehicle within environment and replicate.",
        "e_post": "TSC2ko_rapamycin - WT_rapamycin within environment and replicate.",
        "suppression_fraction": "1 - e_post/e_pre when abs(e_pre) >= min_effect.",
        "categories": {"suppressed": ">=0.75", "partially_retained": "0.25 to <0.75", "fully_retained": "0 to <0.25", "enhanced_or_reprogrammed": "<0"},
        "min_effect": args.min_effect,
        "replicate_rule": "Both GSE179044 replicates must have concordant baseline TSC2-loss direction for a gene to be retention-eligible.",
        "no_patient_claim": "This is perturbation-model evidence, not patient-level sirolimus persistence.",
    }
    (out / "GSE179044_retention_analysis_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(json.dumps({"genes": len(genes), "retention_candidates": int(cross["cross_environment_retention_candidate"].sum()), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
