"""Run targeted Phase 2 checks against the Phase 1 core reproduction."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from scipy.stats import rankdata
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


ROOT = Path(__file__).resolve().parents[1]


def formal_genes(adata: ad.AnnData) -> list[str]:
    table = pd.read_csv(ROOT / "data/raw/reference/LAM_core_signature_genes.csv")
    table.columns = table.columns.astype(str).str.strip()
    requested = table["Gene"].dropna().astype(str).tolist()
    upper = {}
    for actual, symbol in zip(adata.raw.var_names.astype(str), adata.raw.var["gene_symbol_upper"].astype(str)):
        upper.setdefault(symbol.upper(), actual)
    return list(dict.fromkeys(upper[g.upper()] for g in requested if g.upper() in upper))


def filter_summary(obs: pd.DataFrame, retained: pd.Series, label: str) -> pd.DataFrame:
    group_cols = ["condition", "donor_id", "sample_id", "assay"]
    rows = []
    for keys, group in obs.groupby(group_cols, observed=True):
        group_retained = group.loc[retained.loc[group.index]]
        candidate = group_retained["lamcore_candidate_author_style"].astype(bool)
        rows.append({
            "filter": label,
            **dict(zip(group_cols, keys)),
            "retained_cells": int(len(group_retained)),
            "candidate_cells_retained": int(candidate.sum()),
            "candidate_fraction_retained": float(candidate.mean()) if len(candidate) else np.nan,
            "predicted_doublets_retained": int(group_retained["doublet_predicted"].sum()),
            "mean_marker_score": float(group_retained["known_marker_score"].mean()) if len(group_retained) else np.nan,
        })
    return pd.DataFrame(rows)


def donor_program_rows(obs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for donor_id, group in obs[obs["condition"].astype(str).eq("LAM")].groupby("donor_id", observed=True):
        candidate = group["lamcore_candidate_author_style"].astype(bool)
        row = {
            "donor_id": donor_id,
            "assay": str(group["assay"].iloc[0]),
            "cells": int(len(group)),
            "candidate_cells": int(candidate.sum()),
            "candidate_fraction": float(candidate.mean()),
            "candidate_marker_score": float(group.loc[candidate, "known_marker_score"].mean()) if candidate.any() else np.nan,
            "other_marker_score": float(group.loc[~candidate, "known_marker_score"].mean()) if (~candidate).any() else np.nan,
            "candidate_777_score": float(group.loc[candidate, "lamcore_score_777_consistency"].mean()) if candidate.any() else np.nan,
            "other_777_score": float(group.loc[~candidate, "lamcore_score_777_consistency"].mean()) if (~candidate).any() else np.nan,
        }
        row["marker_score_difference"] = row["candidate_marker_score"] - row["other_marker_score"]
        row["score_777_difference"] = row["candidate_777_score"] - row["other_777_score"]
        rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    config = yaml.safe_load((ROOT / "config/reproduction_core.yaml").read_text())
    input_path = ROOT / "data/processed/reproduction_core/GSE135851_core_reproduction.h5ad"
    result_dir = ROOT / "results/robustness"
    table_dir = result_dir / "tables"
    result_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_path)
    obs = adata.obs.copy()
    obs["lamcore_candidate_author_style"] = obs["lamcore_candidate_author_style"].astype(bool)
    base = obs["qc_pass"].astype(bool)
    no_doublet_removal = base
    remove_doublets = base & ~obs["doublet_predicted"].astype(bool)
    limits = obs["mt_pct_limit"].astype(float)
    loose = (obs["n_genes_by_counts"] >= 100) & (obs["total_counts"] >= 300) & (obs["pct_counts_mt"] <= limits + 10)
    strict = (obs["n_genes_by_counts"] >= 500) & (obs["total_counts"] >= 1000) & (obs["pct_counts_mt"] <= np.maximum(limits - 5, 5))
    filters = {
        "phase1_no_doublet_removal": no_doublet_removal,
        "phase2_remove_predicted_doublets": remove_doublets,
        "loose_qc_no_doublet_removal": loose,
        "strict_qc_no_doublet_removal": strict,
    }
    qc_table = pd.concat([filter_summary(obs, retained, label) for label, retained in filters.items()], ignore_index=True)
    qc_table.to_csv(table_dir / "targeted_qc_doublet_sensitivity.csv", index=False)

    # Compare the same fixed author-style candidate labels under alternative scores.
    score_rows = []
    lam = obs[obs["condition"].astype(str).eq("LAM")]
    for label, group in lam.groupby("donor_id", observed=True):
        candidate = group["lamcore_candidate_author_style"].astype(bool)
        score_rows.append({
            "donor_id": label,
            "assay": group["assay"].iloc[0],
            "candidate_cells": int(candidate.sum()),
            "candidate_known_marker_score": float(group.loc[candidate, "known_marker_score"].mean()) if candidate.any() else np.nan,
            "other_known_marker_score": float(group.loc[~candidate, "known_marker_score"].mean()) if (~candidate).any() else np.nan,
            "candidate_777_module_score": float(group.loc[candidate, "lamcore_score_777_consistency"].mean()) if candidate.any() else np.nan,
            "other_777_module_score": float(group.loc[~candidate, "lamcore_score_777_consistency"].mean()) if (~candidate).any() else np.nan,
        })
    score_table = pd.DataFrame(score_rows)
    score_table["known_marker_difference"] = score_table["candidate_known_marker_score"] - score_table["other_known_marker_score"]
    score_table["777_difference_consistency_only"] = score_table["candidate_777_module_score"] - score_table["other_777_module_score"]
    score_table.to_csv(table_dir / "targeted_score_comparison.csv", index=False)

    genes = formal_genes(adata)
    matrix = adata.raw[:, genes].X
    dense = matrix.toarray().astype(np.float32) if sp.issparse(matrix) else np.asarray(matrix, dtype=np.float32)
    ranks = rankdata(dense, axis=1, method="average")
    rank_score = ranks.mean(axis=1) / max(1, dense.shape[1])
    obs["lamcore_rank_score_777"] = rank_score
    rank_rows = []
    for donor_id, group in obs[obs["condition"].astype(str).eq("LAM")].groupby("donor_id", observed=True):
        candidate = group["lamcore_candidate_author_style"].astype(bool)
        rank_rows.append({
            "donor_id": donor_id,
            "candidate_rank_score_777": float(group.loc[candidate, "lamcore_rank_score_777"].mean()) if candidate.any() else np.nan,
            "other_rank_score_777": float(group.loc[~candidate, "lamcore_rank_score_777"].mean()) if (~candidate).any() else np.nan,
            "candidate_cells": int(candidate.sum()),
        })
    pd.DataFrame(rank_rows).to_csv(table_dir / "targeted_rank_based_score_comparison.csv", index=False)

    # Same neighbor graph, different seeds/resolutions: a targeted clustering sensitivity.
    stability_rows = []
    reference = obs["leiden_r0_8"].astype(str).to_numpy()
    for resolution in (0.4, 0.8, 1.2):
        for seed in (0, 1, 2):
            key = f"targeted_leiden_{str(resolution).replace('.', '_')}_{seed}"
            sc.tl.leiden(adata, resolution=resolution, random_state=seed, key_added=key)
            labels = adata.obs[key].astype(str).to_numpy()
            stability_rows.append({
                "resolution": resolution,
                "seed": seed,
                "clusters": int(pd.Series(labels).nunique()),
                "ari_vs_baseline_leiden_r0_8": float(adjusted_rand_score(reference, labels)),
                "nmi_vs_baseline_leiden_r0_8": float(normalized_mutual_info_score(reference, labels)),
            })
    pd.DataFrame(stability_rows).to_csv(table_dir / "targeted_clustering_seed_resolution.csv", index=False)

    # Donor-level and leave-one-donor-out evidence for the fixed candidate set.
    donor_table = donor_program_rows(obs)
    donor_table.to_csv(table_dir / "targeted_donor_level_summary.csv", index=False)
    loo_rows = []
    for excluded in sorted(lam["donor_id"].astype(str).unique()):
        remaining = lam[lam["donor_id"].astype(str) != excluded]
        candidate = remaining["lamcore_candidate_author_style"].astype(bool)
        loo_rows.append({
            "excluded_donor": excluded,
            "remaining_donors": int(remaining["donor_id"].nunique()),
            "remaining_candidate_cells": int(candidate.sum()),
            "remaining_candidate_fraction": float(candidate.mean()),
            "candidate_marker_score": float(remaining.loc[candidate, "known_marker_score"].mean()) if candidate.any() else np.nan,
            "other_marker_score": float(remaining.loc[~candidate, "known_marker_score"].mean()) if (~candidate).any() else np.nan,
            "marker_difference": float(remaining.loc[candidate, "known_marker_score"].mean() - remaining.loc[~candidate, "known_marker_score"].mean()) if candidate.any() and (~candidate).any() else np.nan,
        })
    pd.DataFrame(loo_rows).to_csv(table_dir / "targeted_leave_one_donor_out.csv", index=False)

    assay_rows = (
        lam.groupby("assay", observed=True)
        .agg(cells=("lamcore_candidate_author_style", "size"), candidate_cells=("lamcore_candidate_author_style", "sum"), mean_marker_score=("known_marker_score", "mean"))
        .reset_index()
    )
    assay_rows["candidate_fraction"] = assay_rows["candidate_cells"] / assay_rows["cells"]
    assay_rows.to_csv(table_dir / "targeted_assay_sensitivity.csv", index=False)

    summary = {
        "phase1_cells": int(no_doublet_removal.sum()),
        "phase2_cells_after_doublet_removal": int(remove_doublets.sum()),
        "predicted_doublets": int(obs["doublet_predicted"].sum()),
        "phase1_candidate_cells": int((obs["lamcore_candidate_author_style"] & no_doublet_removal).sum()),
        "phase2_candidate_cells_after_doublet_removal": int((obs["lamcore_candidate_author_style"] & remove_doublets).sum()),
        "loose_candidate_cells": int((obs["lamcore_candidate_author_style"] & loose).sum()),
        "strict_candidate_cells": int((obs["lamcore_candidate_author_style"] & strict).sum()),
        "candidate_by_donor": {str(k): int(v) for k, v in lam.groupby("donor_id", observed=True)["lamcore_candidate_author_style"].sum().items()},
        "candidate_by_assay": {str(k): int(v) for k, v in lam.groupby("assay", observed=True)["lamcore_candidate_author_style"].sum().items()},
        "777_rank_score_status": "completed",
        "ssGSEA_status": "not_run_in_phase2_baseline",
        "interpretation": "targeted sensitivity filter; no mechanical donor pass/fail gate",
    }
    (result_dir / "targeted_robustness_summary.json").write_text(json.dumps(summary, indent=2))
    zh = [
        "# LAMCORE 针对性稳健性验证",
        "",
        "## 目的",
        "",
        "第二阶段针对具体候选结果进行质量过滤，而不是把所有方法比较做成项目终点。Phase 1 候选标签固定后，分别比较 doublet 是否去除、QC 宽严、聚类种子/分辨率、777 module score、rank-based score、assay 分层和 leave-one-donor-out。",
        "",
        f"- Phase 1 不去除 doublet：{summary['phase1_cells']} 个细胞，候选 {summary['phase1_candidate_cells']} 个。",
        f"- 去除预测 doublet：{summary['phase2_cells_after_doublet_removal']} 个细胞，候选 {summary['phase2_candidate_cells_after_doublet_removal']} 个。",
        f"- 预测 doublet 总数：{summary['predicted_doublets']}。",
        f"- 宽松/严格 QC 候选数：{summary['loose_candidate_cells']} / {summary['strict_candidate_cells']}。",
        "",
        "这些结果回答的是候选群对技术处理是否敏感，不等同于独立 donor 验证，也不设置机械的 3/4 donor 通过门槛。",
    ]
    en = [
        "# Targeted Robustness Validation for LAMCORE",
        "",
        "## Purpose",
        "",
        "Phase 2 is a quality filter for concrete candidate findings rather than the endpoint of a methods-comparison project. With the Phase 1 candidate labels fixed, we compare doublet removal, loose/strict QC, clustering seeds/resolutions, the 777-gene module score, a rank-based score, assay strata and leave-one-donor-out analyses.",
        "",
        f"- Phase 1 without doublet removal: {summary['phase1_cells']} cells and {summary['phase1_candidate_cells']} candidates.",
        f"- After removing predicted doublets: {summary['phase2_cells_after_doublet_removal']} cells and {summary['phase2_candidate_cells_after_doublet_removal']} candidates.",
        f"- Predicted doublets: {summary['predicted_doublets']}.",
        f"- Loose/strict QC candidate counts: {summary['loose_candidate_cells']} / {summary['strict_candidate_cells']}.",
        "",
        "These checks assess technical sensitivity; they are not independent donor validation and do not use a mechanical 3/4-donor pass rule.",
    ]
    (result_dir / "LAM_robustness_report_zh.md").write_text("\n".join(zh) + "\n")
    (result_dir / "LAM_robustness_report_en.md").write_text("\n".join(en) + "\n")

    manifest_path = ROOT / "manifests/run_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("steps", [])
    manifest["steps"] = [step for step in manifest["steps"] if step.get("name") != "phase2_targeted_robustness"]
    manifest["steps"].append({"name": "phase2_targeted_robustness", "completed_at": pd.Timestamp.now(tz="UTC").isoformat(), "input": str(input_path.relative_to(ROOT)), "outputs": [str(path.relative_to(ROOT)) for path in sorted(table_dir.glob("targeted_*.csv"))] + [str((result_dir / "targeted_robustness_summary.json").relative_to(ROOT))], "doublet_comparison": "included"})
    manifest["status"] = "phase1_core_reproduction_ready_phase2_completed_phase3_in_progress"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
