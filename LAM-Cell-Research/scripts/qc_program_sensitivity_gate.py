#!/usr/bin/env python3
"""Apply the fixed baseline-versus-strict-QC gate to candidate programs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def read_program_genes(path: Path, source: str) -> pd.DataFrame:
    frame = pd.read_csv(path)
    if source == "pooled":
        required = {"pool", "candidate_program", "gene", "rank_position"}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        frame["donor_id"] = "pooled"
    else:
        required = {"pool", "donor_id", "candidate_program", "gene", "rank_position"}
        missing = required - set(frame.columns)
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
    frame["gene"] = frame["gene"].astype(str).str.upper()
    return frame


def load_candidate_labels(path: Path, strict: bool) -> pd.DataFrame:
    frame = pd.read_csv(path)
    if strict:
        if "qc_pass" not in frame.columns:
            raise ValueError("strict QC requires qc_pass in the candidate table")
        frame = frame[frame["qc_pass"].astype(bool)].copy()
    return frame


def jaccard(a: set[str], b: set[str]) -> float:
    union = a | b
    return len(a & b) / len(union) if union else 0.0


def robust_direction(value: float, tolerance: float = 0.05) -> str:
    """Classify a signed effect without pretending it is inferential evidence."""
    if not np.isfinite(value) or abs(value) < tolerance:
        return "near_zero_or_unavailable"
    return "positive" if value > 0 else "negative"


def summarize_score(scores: pd.Series, labels: pd.Series, donor: pd.Series, mask: pd.Series,
                    program: str, source: str, pool: str = "metadata") -> list[dict]:
    """Return candidate-vs-background summaries for one score under one QC mask."""
    rows: list[dict] = []
    selected = mask & scores.notna()
    for donor_id, idx in donor[selected].groupby(donor[selected], observed=True).groups.items():
        idx = pd.Index(idx)
        candidate = labels.loc[idx].astype(bool)
        value = pd.to_numeric(scores.loc[idx], errors="coerce")
        cand_values = value[candidate]
        bg_values = value[~candidate]
        cand_mean = float(cand_values.mean()) if len(cand_values) else np.nan
        bg_mean = float(bg_values.mean()) if len(bg_values) else np.nan
        effect = cand_mean - bg_mean if np.isfinite(cand_mean) and np.isfinite(bg_mean) else np.nan
        rows.append({
            "source": source,
            "pool": pool,
            "program": program,
            "donor_id": str(donor_id),
            "n_cells_scored": int(value.notna().sum()),
            "candidate_cells_scored": int(len(cand_values)),
            "background_cells_scored": int(len(bg_values)),
            "candidate_mean": cand_mean,
            "background_mean": bg_mean,
            "candidate_minus_background": effect,
            "direction": robust_direction(effect),
        })
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/reproduction_core/GSE135851_core_reproduction.h5ad")
    parser.add_argument("--programs", default="results/program_discovery/pooled_program_genes.csv")
    parser.add_argument("--donor-programs", default="results/program_discovery/donor_wise_program_genes.csv")
    parser.add_argument("--output-dir", default="results/program_discovery/qc_sensitivity")
    parser.add_argument("--top-k", type=int, default=50)
    parser.add_argument("--pooled-scores", default="results/program_discovery/pooled_program_scores.csv")
    parser.add_argument("--core3-scores", default="results/program_discovery/core3_structured_scores.csv")
    args = parser.parse_args()

    out = ROOT / args.output_dir
    out.mkdir(parents=True, exist_ok=True)
    # Only obs metadata is required here; backed mode avoids loading the full
    # expression matrix while preserving the source AnnData provenance.
    adata = ad.read_h5ad(ROOT / args.input, backed="r")
    if "qc_pass" not in adata.obs:
        raise ValueError("Input AnnData has no qc_pass field")
    obs = adata.obs.copy()
    baseline_mask_raw = obs["qc_pass"].astype(bool)
    limits = obs["mt_pct_limit"].astype(float)
    strict_mask_raw = (
        (obs["n_genes_by_counts"] >= 500)
        & (obs["total_counts"] >= 1000)
        & (obs["pct_counts_mt"] <= np.maximum(limits - 5.0, 5.0))
    )
    baseline_labels = obs[baseline_mask_raw].copy()
    strict_labels = obs[strict_mask_raw].copy()
    pooled = read_program_genes(ROOT / args.programs, "pooled")
    donor = read_program_genes(ROOT / args.donor_programs, "donor")

    rows = []
    for source, frame in [("pooled", pooled), ("donor_wise", donor)]:
        grouping = ["pool", "candidate_program"] if source == "pooled" else ["pool", "donor_id", "candidate_program"]
        for keys, group in frame.groupby(grouping, observed=True):
            if not isinstance(keys, tuple):
                keys = (keys,)
            label = dict(zip(grouping, keys))
            genes = set(group.sort_values("rank_position").head(args.top_k)["gene"])
            rows.append({"source": source, **label, "top_k": args.top_k, "genes": ",".join(sorted(genes))})
    programs = pd.DataFrame(rows)

    candidate_col = "lamcore_candidate_author_style"
    if candidate_col not in baseline_labels:
        raise ValueError(f"Input AnnData has no {candidate_col} field")
    candidate_mask = baseline_labels[candidate_col].astype(bool).to_numpy()
    strict_candidate_mask = strict_labels[candidate_col].astype(bool).to_numpy()
    summary_rows = []
    for donor_id in sorted(baseline_labels["donor_id"].astype(str).unique()):
        base_group = baseline_labels[baseline_labels["donor_id"].astype(str).eq(donor_id)]
        strict_group = strict_labels[strict_labels["donor_id"].astype(str).eq(donor_id)]
        summary_rows.append({
            "donor_id": donor_id,
            "baseline_cells": len(base_group),
            "strict_qc_cells": len(strict_group),
            "baseline_candidate_cells": int(base_group[candidate_col].astype(bool).sum()),
            "strict_qc_candidate_cells": int(strict_group[candidate_col].astype(bool).sum()),
            "baseline_candidate_fraction": float(base_group[candidate_col].astype(bool).mean()),
            "strict_qc_candidate_fraction": float(strict_group[candidate_col].astype(bool).mean()) if len(strict_group) else np.nan,
        })

    # Compare pooled and donor-wise program gene lists under the fixed gate.
    program_rows = []
    for _, row in programs.iterrows():
        genes = set(str(row["genes"]).split(",")) if row["genes"] else set()
        related = donor[(donor["pool"] == row["pool"]) & (donor["candidate_program"] == row["candidate_program"])]
        donor_sets = {str(d): set(g.sort_values("rank_position").head(args.top_k)["gene"]) for d, g in related.groupby("donor_id", observed=True)}
        pairwise = [jaccard(genes, donor_genes) for donor_genes in donor_sets.values()]
        program_rows.append({
            "source": row["source"],
            "pool": row["pool"],
            "donor_id": row.get("donor_id", "pooled"),
            "candidate_program": row["candidate_program"],
            "top_k": args.top_k,
            "baseline_candidate_cells": int(candidate_mask.sum()),
            "strict_qc_candidate_cells": int(strict_candidate_mask.sum()),
            "donorwise_to_pooled_median_jaccard": float(np.median(pairwise)) if pairwise else np.nan,
            "donorwise_to_pooled_max_jaccard": float(np.max(pairwise)) if pairwise else np.nan,
            "qc_gate_status": "gene_list_comparison_only",
        })

    # Score-level sensitivity. The core reproduction object already stores the
    # state scores used by the baseline analysis. Reusing these columns avoids
    # silently changing normalization while still testing the fixed QC masks.
    state_score_columns = {
        "contractile": "state_contractile",
        "LAM_myogenic_contractile": "state_contractile",
        "lineage_uterine_smooth_muscle": "state_contractile",
        "CORE1": "state_contractile",
        "ECM_remodeling": "state_ecm_remodeling",
        "TGFbeta_fibroblast": "microenv_fibroblast",
        "normal_lung_interstitial": "microenv_fibroblast",
        "SLS_stem_like": "state_metabolic",
        "IS_inflammatory": "state_inflammatory",
        "MDK_dormancy_persistence": "state_stress_hypoxia",
        "hypoxia_stress": "state_stress_hypoxia",
        "cell_cycle": "state_proliferative",
        "hormone_related": "state_hormone_related",
        "mTOR_translation": "state_mTOR_related",
        "LAF_niche": "microenv_fibroblast",
        "macrophage_TREM2_TYROBP": "microenv_immune",
        "IL6_AT2_repair": "microenv_immune",
        "protease_ECM_niche": "state_ecm_remodeling",
        "biomarker_VEGFD_PMEL_CCL14_MMP8": "known_marker_score",
        "CORE2": "known_marker_score",
    }
    score_rows: list[dict] = []
    score_columns_used: dict[str, str] = {}
    base_mask = pd.Series(baseline_mask_raw.to_numpy(), index=obs.index)
    strict_mask = pd.Series(strict_mask_raw.to_numpy(), index=obs.index)
    labels = obs[candidate_col].astype(bool)
    donors = obs["donor_id"].astype(str)
    for program, column in state_score_columns.items():
        if column not in obs:
            continue
        score_columns_used[program] = column
        scores = pd.to_numeric(obs[column], errors="coerce")
        score_rows.extend(summarize_score(scores, labels, donors, base_mask, program, "adata_obs", "all_cells"))
        score_rows.extend(summarize_score(scores, labels, donors, strict_mask, program, "adata_obs", "all_cells"))

    # The NMF scores are only available for the candidate pools that were used
    # to discover them. They are therefore reported as a complementary score
    # shift, not as a fresh QC-aware discovery run.
    pooled_score_path = ROOT / args.pooled_scores
    if pooled_score_path.exists():
        nmf_scores = pd.read_csv(pooled_score_path)
        nmf_scores["cell_id"] = nmf_scores["cell_id"].astype(str)
        nmf_scores = nmf_scores.set_index("cell_id")
        common = obs.index.astype(str).intersection(nmf_scores.index)
        if len(common):
            local_obs = obs.copy()
            local_obs.index = local_obs.index.astype(str)
            for (pool, program), group in nmf_scores.groupby(["pool", "candidate_program"], observed=True):
                # Scores are per cell in a candidate pool. Summarize score
                # distributions and candidate fraction under both QC masks.
                for label, mask in [("baseline", base_mask), ("strict_qc", strict_mask)]:
                    use_ids = common[mask.reindex(common).fillna(False).to_numpy()]
                    values = group.reindex(use_ids)["score"].dropna()
                    score_rows.append({
                        "source": "pooled_nmf_score",
                        "pool": str(pool),
                        "program": str(program),
                        "donor_id": "pooled",
                        "qc_view": label,
                        "n_cells_scored": int(len(values)),
                        "candidate_cells_scored": np.nan,
                        "background_cells_scored": np.nan,
                        "candidate_mean": float(values.mean()) if len(values) else np.nan,
                        "background_mean": np.nan,
                        "candidate_minus_background": np.nan,
                        "direction": "score_distribution_only",
                        "score_median": float(values.median()) if len(values) else np.nan,
                        "score_q25": float(values.quantile(0.25)) if len(values) else np.nan,
                        "score_q75": float(values.quantile(0.75)) if len(values) else np.nan,
                    })

    score_frame = pd.DataFrame(score_rows)
    if not score_frame.empty:
        # The adata_obs summaries do not carry this label yet because they are
        # generated by a common helper; infer it from the cell count mask by
        # writing an explicit row-level view here.
        if "qc_view" not in score_frame.columns:
            score_frame["qc_view"] = np.nan
        score_frame.loc[score_frame["source"].eq("adata_obs") & score_frame["qc_view"].isna(), "qc_view"] = "unspecified"
        score_frame.to_csv(out / "program_qc_gate_scores.csv", index=False)

        # Direction comparison is only made where both views have an effect.
        direction_rows = []
        adata_scores = score_frame[score_frame["source"].eq("adata_obs")].copy()
        # Helper output has two rows per donor/program; identify baseline vs
        # strict by recomputing the expected row order from cell counts.
        for program, column in score_columns_used.items():
            p = adata_scores[adata_scores["program"].eq(program)].copy()
            # Re-run the compact summary with an explicit QC label.
            for qc_view, mask in [("baseline", base_mask), ("strict_qc", strict_mask)]:
                compact = summarize_score(pd.to_numeric(obs[column], errors="coerce"), labels, donors, mask, program, "adata_obs", "all_cells")
                for row in compact:
                    row["qc_view"] = qc_view
                    direction_rows.append(row)
        direction = pd.DataFrame(direction_rows)
        if not direction.empty:
            direction.to_csv(out / "program_qc_gate_direction_by_donor.csv", index=False)
            comparison_rows = []
            for (program, donor_id), group in direction.groupby(["program", "donor_id"], observed=True):
                views = group.set_index("qc_view")
                b = views.loc["baseline"] if "baseline" in views.index else None
                s = views.loc["strict_qc"] if "strict_qc" in views.index else None
                bdir = b["direction"] if b is not None else "missing"
                sdir = s["direction"] if s is not None else "missing"
                comparison_rows.append({
                    "program": program,
                    "donor_id": donor_id,
                    "baseline_effect": b["candidate_minus_background"] if b is not None else np.nan,
                    "strict_qc_effect": s["candidate_minus_background"] if s is not None else np.nan,
                    "baseline_direction": bdir,
                    "strict_qc_direction": sdir,
                    "direction_consistent": bool(bdir == sdir and bdir in {"positive", "negative"}),
                    "qc_direction_status": "qc_stable" if bdir == sdir and bdir in {"positive", "negative"} else "qc_sensitive_or_insufficient",
                })
            pd.DataFrame(comparison_rows).to_csv(out / "program_qc_gate_direction_comparison.csv", index=False)

    pd.DataFrame(summary_rows).to_csv(out / "baseline_strict_qc_cell_summary.csv", index=False)
    pd.DataFrame(program_rows).to_csv(out / "program_qc_gate_summary.csv", index=False)
    manifest = {
        "input": args.input,
        "programs": args.programs,
        "donor_programs": args.donor_programs,
        "top_k": args.top_k,
        "baseline_definition": "all cells in the core reproduction AnnData",
        "strict_qc_definition": "n_genes_by_counts >= 500 and total_counts >= 1000 and pct_counts_mt <= max(mt_pct_limit - 5, 5), matching scripts/run_targeted_robustness.py",
        "fixed_gate": "Every Hypothesis Card must compare baseline and strict-QC candidate count, donor distribution, program direction and loading before upgrade.",
        "score_columns_used": score_columns_used,
        "score_outputs": ["program_qc_gate_scores.csv", "program_qc_gate_direction_by_donor.csv", "program_qc_gate_direction_comparison.csv"],
        "note": "State scores stored in the baseline AnnData are summarized under both masks; pooled NMF scores are reported as score-distribution sensitivity only, not a fresh QC-aware NMF run.",
    }
    (out / "qc_sensitivity_gate_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False))
    print(json.dumps({"baseline_cells": len(baseline_labels), "strict_qc_cells": len(strict_labels), "baseline_candidates": int(candidate_mask.sum()), "strict_qc_candidates": int(strict_candidate_mask.sum()), "output_dir": str(out)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
