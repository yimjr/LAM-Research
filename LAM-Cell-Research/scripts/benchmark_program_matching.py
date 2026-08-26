#!/usr/bin/env python3
"""Benchmark donor-wise program matching before interpreting unknown programs.

The benchmark asks whether the current discovery/matching representation can
recover programs that are already expected to be present. It is a method
calibration layer, not evidence that any tested program is a new LAM state.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]

DEFAULT_PROGRAMS = {
    "CORE1": ["ACTA2", "MYH11", "TAGLN", "CNN1", "TPM2", "DES", "RGS5", "COL4A1"],
    "CORE2": ["PMEL", "MLANA", "MITF", "TYR", "DCT", "GPNMB", "CTSK", "FIGF"],
    "contractile": ["ACTA2", "TAGLN", "MYH11", "CNN1", "TPM2", "DES", "CALD1", "MYL9"],
    "ECM_remodeling": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "FN1", "VCAN", "MMP2", "MMP9", "CTSK", "TIMP1"],
    "SLS_stem_like": ["MDK", "SOX2", "PROM1", "NES", "LGR5", "EPCAM", "KLF4", "NANOG", "DPPA4", "HMGA2"],
    "IS_inflammatory": ["IL6", "IL1B", "CXCL8", "CCL2", "CCL7", "NFKBIA", "TNFAIP3", "STAT1", "ISG15", "IFIT1"],
    "MDK_dormancy_persistence": ["MDK", "LGALS3", "BCL2", "BCL2L1", "JUN", "FOS", "DDIT4", "TXNIP", "KLF6", "NR4A1"],
}


def load_known_programs(path: Path) -> dict[str, list[str]]:
    if not path.exists():
        return DEFAULT_PROGRAMS
    payload = yaml.safe_load(path.read_text()) or {}
    programs = {}
    for item in payload.get("programs", []):
        name = str(item.get("program_name", ""))
        genes = [str(g).upper() for g in item.get("genes", [])]
        if name and genes:
            programs[name] = genes
    for name, genes in DEFAULT_PROGRAMS.items():
        programs.setdefault(name, genes)
    return programs


def jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def random_null(left_size: int, right_size: int, universe: list[str], n: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    if not left_size or not right_size or len(universe) < max(left_size, right_size):
        return np.asarray([], dtype=float)
    values = []
    for _ in range(n):
        a = set(rng.choice(universe, size=left_size, replace=False))
        b = set(rng.choice(universe, size=right_size, replace=False))
        values.append(jaccard(a, b))
    null = np.asarray(values, dtype=float)
    return null


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/reproduction_core/GSE135851_core_reproduction.h5ad")
    parser.add_argument("--known-programs", default="config/known_lam_programs.yaml")
    parser.add_argument("--pooled-programs", default="results/program_discovery/pooled_program_genes.csv")
    parser.add_argument("--donor-programs", default="results/program_discovery/donor_wise_program_genes.csv")
    parser.add_argument("--output-dir", default="results/program_discovery/benchmark")
    parser.add_argument("--top-k", type=int, nargs="+", default=[25, 50, 100])
    parser.add_argument("--null-iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260823)
    args = parser.parse_args()

    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    programs = load_known_programs(ROOT / args.known_programs)
    pooled = pd.read_csv(ROOT / args.pooled_programs)
    donor = pd.read_csv(ROOT / args.donor_programs)
    pooled = pooled[pooled["pool"].isin(["high_confidence", "broad_lam_like", "unrestricted_lam"])].copy()
    donor = donor[donor["pool"].isin(["high_confidence", "broad_lam_like", "unrestricted_lam"])].copy()

    # This benchmark operates on the already-generated program tables. Avoid
    # loading the multi-gigabyte AnnData matrix just to obtain a gene universe.
    # The source AnnData path remains in the manifest for provenance.
    gene_universe = sorted(
        set(pooled["gene"].astype(str).str.upper())
        | set(donor["gene"].astype(str).str.upper())
        | {g for genes in programs.values() for g in genes}
    )
    null_cache: dict[tuple[int, int], np.ndarray] = {}
    rows = []
    for program_name, genes in programs.items():
        for pool, group in pooled.groupby("pool", observed=True):
            for candidate, pgroup in group.groupby("candidate_program", observed=True):
                ranked = pgroup.sort_values("rank_position")["gene"].astype(str).str.upper().tolist()
                for top_k in args.top_k:
                    known_set = set(str(g).upper() for g in genes)
                    candidate_set = set(ranked[:top_k])
                    overlap = len(known_set & candidate_set)
                    observed = overlap / len(known_set | candidate_set) if known_set | candidate_set else 0.0
                    cache_key = (len(known_set), len(candidate_set))
                    if cache_key not in null_cache:
                        null_cache[cache_key] = random_null(*cache_key, gene_universe, args.null_iterations, args.seed + top_k)
                    null = null_cache[cache_key]
                    null_q95 = float(np.quantile(null, 0.95)) if len(null) else np.nan
                    null_p = float((null >= observed).mean()) if len(null) else np.nan
                    rows.append({
                        "comparison": "pooled_to_known",
                        "pool": pool,
                        "donor_id": "pooled",
                        "candidate_program": candidate,
                        "known_program": program_name,
                        "top_k": top_k,
                        "known_gene_count": len(known_set),
                        "candidate_gene_count": len(candidate_set),
                        "overlap_genes": overlap,
                        "jaccard": observed,
                        "null_q95": null_q95,
                        "null_empirical_p": null_p,
                        "benchmark_detectable": bool(np.isfinite(null_q95) and observed > null_q95),
                    })
        for (pool, donor_id, candidate), group in donor.groupby(["pool", "donor_id", "candidate_program"], observed=True):
            ranked = group.sort_values("rank_position")["gene"].astype(str).str.upper().tolist()
            for top_k in args.top_k:
                known_set = set(str(g).upper() for g in genes)
                candidate_set = set(ranked[:top_k])
                overlap = len(known_set & candidate_set)
                observed = overlap / len(known_set | candidate_set) if known_set | candidate_set else 0.0
                cache_key = (len(known_set), len(candidate_set))
                if cache_key not in null_cache:
                    null_cache[cache_key] = random_null(*cache_key, gene_universe, args.null_iterations, args.seed + top_k)
                null = null_cache[cache_key]
                null_q95 = float(np.quantile(null, 0.95)) if len(null) else np.nan
                null_p = float((null >= observed).mean()) if len(null) else np.nan
                rows.append({
                    "comparison": "donor_to_known",
                    "pool": pool,
                    "donor_id": donor_id,
                    "candidate_program": candidate,
                    "known_program": program_name,
                    "top_k": top_k,
                    "known_gene_count": len(known_set),
                    "candidate_gene_count": len(candidate_set),
                    "overlap_genes": overlap,
                    "jaccard": observed,
                    "null_q95": null_q95,
                    "null_empirical_p": null_p,
                    "benchmark_detectable": bool(np.isfinite(null_q95) and observed > null_q95),
                })

    benchmark = pd.DataFrame(rows)
    benchmark.to_csv(out / "known_program_matching_benchmark.csv", index=False)
    summary = benchmark.groupby(["comparison", "known_program", "top_k"], dropna=False).agg(
        rows=("jaccard", "size"),
        median_jaccard=("jaccard", "median"),
        max_jaccard=("jaccard", "max"),
        detectable_fraction=("benchmark_detectable", "mean"),
    ).reset_index()
    summary.to_csv(out / "known_program_matching_summary.csv", index=False)
    manifest = {
        "input": args.input,
        "known_program_source": args.known_programs,
        "pooled_program_source": args.pooled_programs,
        "donor_program_source": args.donor_programs,
        "top_k": args.top_k,
        "null_iterations": args.null_iterations,
        "seed": args.seed,
        "interpretation": "Positive-control benchmark for method sensitivity; does not establish a biological claim.",
        "unknown_program_interpretation_gate": "Do not interpret weak cross-PatientID matching as patient-specific biology until known-program benchmark is adequate.",
    }
    (out / "benchmark_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(json.dumps({"benchmark_rows": len(benchmark), "known_programs": len(programs), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
