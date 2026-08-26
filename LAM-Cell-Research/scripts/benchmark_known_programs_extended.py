#!/usr/bin/env python3
"""Extended positive-control benchmark using expression and NMF loading links."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]


def load_programs(path: Path) -> dict[str, list[str]]:
    payload = yaml.safe_load(path.read_text()) or {}
    return {str(x["program_name"]): [str(g).upper() for g in x.get("genes", [])] for x in payload.get("programs", []) if x.get("program_name") and x.get("genes")}


def gene_symbols(adata: ad.AnnData) -> np.ndarray:
    for col in ["gene_symbol_upper", "gene_symbol"]:
        if col in adata.var:
            return adata.var[col].astype(str).str.upper().to_numpy()
    return adata.var_names.astype(str).str.upper().to_numpy()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/reproduction_core/GSE135851_core_reproduction.h5ad")
    parser.add_argument("--known-programs", default="config/known_lam_programs.yaml")
    parser.add_argument("--pooled-scores", default="results/program_discovery/pooled_program_scores.csv")
    parser.add_argument("--pooled-genes", default="results/program_discovery/pooled_program_genes.csv")
    parser.add_argument("--output-dir", default="results/program_discovery/benchmark")
    args = parser.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    programs = load_programs(ROOT / args.known_programs)
    adata = ad.read_h5ad(ROOT / args.input, backed="r")
    symbols = gene_symbols(adata)
    known_union = sorted(set(sum(programs.values(), [])))
    indices = np.flatnonzero(np.isin(symbols, known_union))
    raw = adata[:, indices].layers["counts"] if "counts" in adata.layers else adata[:, indices].X
    raw = raw.toarray() if hasattr(raw, "toarray") else np.asarray(raw)
    raw = np.asarray(raw, dtype=float)
    totals = pd.to_numeric(adata.obs.get("total_counts", pd.Series(raw.sum(axis=1))), errors="coerce").to_numpy(float)
    norm = np.log1p(raw / np.maximum(totals[:, None], 1.0) * 1e4)
    selected_symbols = symbols[indices]
    obs = adata.obs.copy()
    obs["_cell_id"] = obs["cell_id"].astype(str) if "cell_id" in obs else obs.index.astype(str)
    obs["_donor"] = obs["donor_id"].astype(str)
    obs["_assay"] = obs["assay"].astype(str)
    obs["_candidate"] = obs.get("lamcore_candidate_author_style", False).astype(bool)
    score_table = pd.DataFrame({"cell_id": obs["_cell_id"].to_numpy(), "donor_id": obs["_donor"].to_numpy(), "assay": obs["_assay"].to_numpy(), "candidate": obs["_candidate"].to_numpy()})
    for program, genes in programs.items():
        rows = np.flatnonzero(np.isin(selected_symbols, genes))
        score_table[program] = np.nanmean(norm[:, rows], axis=1) if len(rows) else np.nan

    donor_rows = []
    for (donor_id, assay, program), group in score_table.melt(id_vars=["cell_id", "donor_id", "assay", "candidate"], var_name="known_program", value_name="score").groupby(["donor_id", "assay", "known_program"], observed=True):
        cand = group.loc[group["candidate"], "score"].dropna()
        bg = group.loc[~group["candidate"], "score"].dropna()
        donor_rows.append({
            "donor_id": donor_id,
            "assay": assay,
            "known_program": program,
            "candidate_cells": int(len(cand)),
            "background_cells": int(len(bg)),
            "candidate_mean": float(cand.mean()) if len(cand) else np.nan,
            "background_mean": float(bg.mean()) if len(bg) else np.nan,
            "candidate_minus_background": float(cand.mean() - bg.mean()) if len(cand) and len(bg) else np.nan,
            "expression_direction": "positive" if len(cand) and len(bg) and cand.mean() > bg.mean() else "non_positive_or_insufficient",
        })
    donor_expression = pd.DataFrame(donor_rows)
    donor_expression.to_csv(out / "known_program_donor_expression_scores.csv", index=False)

    pooled_scores = pd.read_csv(ROOT / args.pooled_scores)
    known_by_cell = score_table.set_index("cell_id")
    loading_rows = []
    for (pool, candidate_program), group in pooled_scores.groupby(["pool", "candidate_program"], observed=True):
        common = group["cell_id"].astype(str).isin(known_by_cell.index)
        joined = group.loc[common, ["cell_id", "score", "donor_id"]].copy()
        if joined.empty:
            continue
        known = known_by_cell.loc[joined["cell_id"].astype(str)]
        for known_program in programs:
            x = joined["score"].to_numpy(float)
            y = known[known_program].to_numpy(float)
            finite = np.isfinite(x) & np.isfinite(y)
            corr = float(pd.Series(x[finite]).corr(pd.Series(y[finite]), method="spearman")) if finite.sum() >= 10 else np.nan
            loading_rows.append({
                "pool": pool,
                "candidate_program": candidate_program,
                "known_program": known_program,
                "n_common_cells": int(finite.sum()),
                "spearman_nmf_score_vs_known_expression": corr,
                "loading_similarity_status": "positive_association" if np.isfinite(corr) and corr >= 0.3 else "weak_or_unavailable",
            })
    pd.DataFrame(loading_rows).to_csv(out / "known_program_nmf_loading_similarity.csv", index=False)

    # Program-level sufficiency is deliberately conservative. The benchmark
    # is considered biologically usable only when expression recovery and at
    # least one representation-level link are present; otherwise unknown
    # cross-PatientID matching remains gated.
    summary_rows = []
    for program in programs:
        expr = donor_expression[(donor_expression["known_program"].eq(program)) & donor_expression["donor_id"].isin(["LAM1", "LAM2", "LAM3", "LAM4"])]
        expr = expr[expr["candidate_cells"] >= 3]
        positive_fraction = float((expr["expression_direction"].eq("positive")).mean()) if len(expr) else np.nan
        load = pd.DataFrame(loading_rows)
        load = load[load["known_program"].eq(program)] if not load.empty else load
        loading_positive = bool((load["loading_similarity_status"].eq("positive_association")).any()) if not load.empty else False
        expression_pass = bool(np.isfinite(positive_fraction) and positive_fraction >= 0.75)
        summary_rows.append({
            "known_program": program,
            "lam_donors_tested": int(len(expr)),
            "expression_positive_fraction": positive_fraction,
            "expression_benchmark_pass": expression_pass,
            "loading_similarity_positive": loading_positive,
            "regulon_status": "not_run_current_runtime",
            "pathway_status": "not_run_current_runtime",
            "benchmark_sufficient_for_unknown_interpretation": bool(expression_pass and loading_positive),
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(out / "known_program_extended_benchmark_summary.csv", index=False)
    core_required = [x for x in ["CORE1", "CORE2", "ECM_remodeling", "SLS_stem_like", "IS_inflammatory"] if x in summary["known_program"].tolist()]
    manifest = {
        "input": args.input,
        "methods": ["donor_level_expression_enrichment", "pooled_NMF_score_vs_known_expression_spearman", "existing_top_gene_jaccard_benchmark"],
        "required_programs_for_unknown_gate": core_required,
        "unknown_gate": "Unknown cross-PatientID matches remain gated unless the relevant known positive controls pass expression and loading benchmarks.",
        "regulon_pathway_status": "not_run_current_runtime; requires curated regulon/pathway resource before final novelty claims.",
        "CORE3_rule": "CORE3 remains a structured identity+depth-adjusted-low-activity+translation model, not an ordinary gene set.",
    }
    (out / "extended_benchmark_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(json.dumps({"programs": len(programs), "donor_expression_rows": len(donor_expression), "loading_rows": len(loading_rows), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
