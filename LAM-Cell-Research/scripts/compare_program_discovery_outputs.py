#!/usr/bin/env python3
"""Match candidate programs across independent discovery runs.

This compares top-gene sets only; it does not convert overlap into biological
confirmation. PatientID overlap is explicitly recorded so same-donor RNA and
orthogonal assay evidence cannot be counted as independent replication.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import pandas as pd
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]

DATASETS = {
    "GSE135851_core": {"dir": "results/program_discovery", "patient_ids": ["LAM1", "LAM2", "LAM3", "LAM4"]},
    "GSE190260": {"dir": "results/program_discovery/external_GSE190260", "patient_ids": ["LAM1110", "LAM1158", "LAM1163", "LAM1164"]},
    "GSE217108": {"dir": "results/program_discovery/external_GSE217108", "patient_ids": ["LAM32", "LAM44"]},
    "GSE302356": {"dir": "results/program_discovery/external_GSE302356", "patient_ids": ["LAM32", "LAM18", "LAM3", "LAM50"]},
}


def jaccard(a: set[str], b: set[str]) -> float:
    return len(a & b) / len(a | b) if a | b else 0.0


def load_programs(dataset: str, pool: str) -> dict[str, set[str]]:
    path = ROOT / DATASETS[dataset]["dir"] / "pooled_program_genes.csv"
    if not path.exists():
        return {}
    frame = pd.read_csv(path)
    frame = frame[frame["pool"] == pool]
    return {name: set(group.sort_values("rank_position").head(50)["gene"].astype(str).str.upper()) for name, group in frame.groupby("candidate_program")}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/program_discovery/cross_dataset")
    parser.add_argument("--threshold", type=float, default=0.15)
    args = parser.parse_args()
    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    pair_rows = []
    meta_rows = []
    for pool in ["high_confidence", "broad_lam_like", "unrestricted_lam"]:
        loaded = {dataset: load_programs(dataset, pool) for dataset in DATASETS}
        for left, right in combinations(DATASETS, 2):
            lp, rp = loaded[left], loaded[right]
            if not lp or not rp:
                continue
            left_names, right_names = list(lp), list(rp)
            scores = [[jaccard(lp[a], rp[b]) for b in right_names] for a in left_names]
            rows, cols = linear_sum_assignment([[-x for x in row] for row in scores])
            for i, j in zip(rows, cols):
                shared = sorted(lp[left_names[i]] & rp[right_names[j]])
                patient_overlap = sorted(set(DATASETS[left]["patient_ids"]) & set(DATASETS[right]["patient_ids"]))
                pair_rows.append({"pool": pool, "dataset_left": left, "program_left": left_names[i], "dataset_right": right, "program_right": right_names[j], "top_gene_jaccard": float(scores[i][j]), "shared_genes": ",".join(shared), "patient_id_overlap": ",".join(patient_overlap), "independence_note": "same_patient_overlap_present" if patient_overlap else "different_patient_sets"})
        for dataset, programs in loaded.items():
            for program, genes in programs.items():
                matches = [row for row in pair_rows if row["pool"] == pool and ((row["dataset_left"] == dataset and row["program_left"] == program) or (row["dataset_right"] == dataset and row["program_right"] == program))]
                qualifying = [row for row in matches if row["top_gene_jaccard"] >= args.threshold]
                independent = [row for row in qualifying if row["independence_note"] == "different_patient_sets"]
                meta_rows.append({"pool": pool, "dataset": dataset, "program": program, "n_pairwise_matches_ge_threshold": len(qualifying), "n_matches_with_different_patient_sets": len(independent), "matched_datasets": ",".join(sorted({row["dataset_right"] if row["dataset_left"] == dataset else row["dataset_left"] for row in qualifying})), "max_top_gene_jaccard": max([row["top_gene_jaccard"] for row in matches], default=0.0), "classification": "candidate_cross_dataset_program" if len(independent) >= 1 else "same_cohort_or_same_patient_overlap_only"})
    pair = pd.DataFrame(pair_rows)
    meta = pd.DataFrame(meta_rows)
    pair.to_csv(out / "cross_dataset_program_matches.csv", index=False)
    meta.to_csv(out / "cross_dataset_meta_programs.csv", index=False)
    manifest = {"datasets": DATASETS, "top_gene_threshold": args.threshold, "interpretation": "Gene-set overlap is a candidate matching signal; independent PatientID and orthogonal assay evidence must be assessed separately.", "outputs": ["cross_dataset_program_matches.csv", "cross_dataset_meta_programs.csv"]}
    (out / "cross_dataset_comparison_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(json.dumps({"pair_matches": len(pair), "meta_rows": len(meta), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
