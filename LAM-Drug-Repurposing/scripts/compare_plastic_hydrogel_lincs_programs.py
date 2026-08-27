"""Compare independent plastic-panel and hydrogel-panel LINCS analyses."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
PROGRAMS = ROOT / "data/processed/candidate_analysis/programs"
REPORTS = ROOT / "reports/06_candidate_analysis"
SIGNATURES = ROOT / "results/signatures/GSE179044_cmap_query_signatures.csv"
PREFIX = "tsc2_loss_plastic_or_hydrogel_replicated_concordant"
PLASTIC_SUMMARY = PROGRAMS / f"{PREFIX}_LINCS_drug_gene_response_summary.csv.gz"
HYDROGEL_SUMMARY = PROGRAMS / f"{PREFIX}_hydrogel_panel_LINCS_drug_gene_response_summary.csv.gz"


def top_panel(signatures: pd.DataFrame, contrast: str) -> dict[str, set[str]]:
    group = signatures.loc[signatures["contrast"].eq(contrast)].copy()
    if group.empty:
        raise ValueError(f"missing signature contrast: {contrast}")
    result: dict[str, set[str]] = {}
    for direction, ascending in (("up", False), ("down", True)):
        selected = (
            group.loc[group["direction"].eq(direction)]
            .sort_values(["moderated_q", "signed_score"], ascending=[True, ascending])
            .head(150)
        )
        result[direction] = set(selected["gene"].astype(str).str.upper().str.strip())
    return result


def safe_jaccard(left: set[str], right: set[str]) -> float:
    union = left | right
    return len(left & right) / len(union) if union else np.nan


def safe_overlap_coefficient(left: set[str], right: set[str]) -> float:
    denominator = min(len(left), len(right))
    return len(left & right) / denominator if denominator else np.nan


def compare_genes(plastic: pd.DataFrame, hydrogel: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    p = plastic.loc[plastic["entity_type"].eq("compound")].copy()
    h = hydrogel.loc[hydrogel["entity_type"].eq("compound")].copy()
    keys = ["entity_id", "dataset", "scope", "gene"]
    merged = p.merge(h, on=keys, how="outer", suffixes=("_plastic", "_hydrogel"), validate="one_to_one", indicator=True)
    if merged.empty:
        raise ValueError("no common compound × dataset × scope × gene rows")
    merged["pert_iname"] = merged["pert_iname_plastic"].fillna(merged["pert_iname_hydrogel"])

    status_p = merged["direction_status_plastic"].fillna("not_available").astype(str)
    status_h = merged["direction_status_hydrogel"].fillna("not_available").astype(str)
    merged["direction_pair"] = status_p + "__" + status_h
    merged["panel_presence"] = merged["_merge"].map({"both": "both_panels", "left_only": "plastic_only_panel", "right_only": "hydrogel_only_panel"})
    merged["status_same"] = status_p.eq(status_h)
    both_panels = merged["_merge"].eq("both")
    nonneutral = both_panels & status_p.isin(["reversal", "mimic"]) & status_h.isin(["reversal", "mimic"])
    merged["nonneutral_direction_concordant"] = nonneutral & status_p.eq(status_h)
    merged["direction_switch"] = (status_p.eq("reversal") & status_h.eq("mimic")) | (status_p.eq("mimic") & status_h.eq("reversal"))
    p_effect = pd.to_numeric(merged["drug_effect_median_plastic"], errors="coerce")
    h_effect = pd.to_numeric(merged["drug_effect_median_hydrogel"], errors="coerce")
    merged["effect_delta_hydrogel_minus_plastic"] = h_effect - p_effect
    merged["effect_sign_same"] = np.sign(p_effect).eq(np.sign(h_effect)) & p_effect.notna() & h_effect.notna()

    rows = []
    for (entity_id, dataset, scope), group in merged.groupby(["entity_id", "dataset", "scope"], sort=True):
        p_status = group["direction_status_plastic"].fillna("not_available").astype(str)
        h_status = group["direction_status_hydrogel"].fillna("not_available").astype(str)
        both_panel_rows = group["panel_presence"].eq("both_panels")
        p_reversal = set(group.loc[p_status.eq("reversal"), "gene"])
        h_reversal = set(group.loc[h_status.eq("reversal"), "gene"])
        p_mimic = set(group.loc[p_status.eq("mimic"), "gene"])
        h_mimic = set(group.loc[h_status.eq("mimic"), "gene"])
        p_effect = pd.to_numeric(group["drug_effect_median_plastic"], errors="coerce").to_numpy(float)
        h_effect = pd.to_numeric(group["drug_effect_median_hydrogel"], errors="coerce").to_numpy(float)
        valid = np.isfinite(p_effect) & np.isfinite(h_effect)
        rho = float(spearmanr(p_effect[valid], h_effect[valid]).statistic) if valid.sum() >= 3 else np.nan
        both_non_neutral = both_panel_rows & p_status.isin(["reversal", "mimic"]) & h_status.isin(["reversal", "mimic"])
        rows.append(
            {
                "entity_id": entity_id,
                "pert_iname": group["pert_iname"].iloc[0],
                "dataset": dataset,
                "scope": scope,
                "n_genes_union": len(group),
                "n_genes_in_both_panels": int(both_panel_rows.sum()),
                "n_genes_plastic_only_panel": int((group["panel_presence"] == "plastic_only_panel").sum()),
                "n_genes_hydrogel_only_panel": int((group["panel_presence"] == "hydrogel_only_panel").sum()),
                "n_status_same": int(group["status_same"].sum()),
                "status_agreement": float(group["status_same"].mean()),
                "n_status_same_both_panels": int(group.loc[both_panel_rows, "status_same"].sum()),
                "status_agreement_both_panels": float(group.loc[both_panel_rows, "status_same"].mean()) if both_panel_rows.any() else np.nan,
                "n_non_neutral_both": int(both_non_neutral.sum()),
                "n_non_neutral_direction_concordant": int(group["nonneutral_direction_concordant"].sum()),
                "non_neutral_direction_concordance": float(group.loc[both_non_neutral, "nonneutral_direction_concordant"].mean()) if both_non_neutral.any() else np.nan,
                "n_direction_switch": int(group["direction_switch"].sum()),
                "n_reversal_plastic": len(p_reversal),
                "n_reversal_hydrogel": len(h_reversal),
                "n_reversal_both": len(p_reversal & h_reversal),
                "reversal_jaccard": safe_jaccard(p_reversal, h_reversal),
                "reversal_overlap_coefficient": safe_overlap_coefficient(p_reversal, h_reversal),
                "n_mimic_plastic": len(p_mimic),
                "n_mimic_hydrogel": len(h_mimic),
                "n_mimic_both": len(p_mimic & h_mimic),
                "mimic_jaccard": safe_jaccard(p_mimic, h_mimic),
                "n_plastic_only_reversal": len(p_reversal - h_reversal),
                "n_hydrogel_only_reversal": len(h_reversal - p_reversal),
                "n_effects_compared": int(valid.sum()),
                "effect_spearman": rho,
                "effect_sign_concordance": float(group.loc[valid, "effect_sign_same"].mean()) if valid.any() else np.nan,
                "mean_effect_delta_hydrogel_minus_plastic": float(np.nanmean(h_effect - p_effect)) if valid.any() else np.nan,
            }
        )
    summary = pd.DataFrame(rows)
    return merged, summary


def reversal_set_comparison(gene_comparison: pd.DataFrame) -> pd.DataFrame:
    """Return the actual reversal sets for each drug × LINCS release."""
    rows = []
    for (entity_id, dataset, scope), group in gene_comparison.groupby(["entity_id", "dataset", "scope"], sort=True):
        p_status = group["direction_status_plastic"].fillna("not_available").astype(str)
        h_status = group["direction_status_hydrogel"].fillna("not_available").astype(str)
        p_reversal = set(group.loc[p_status.eq("reversal"), "gene"].astype(str))
        h_reversal = set(group.loc[h_status.eq("reversal"), "gene"].astype(str))
        shared = p_reversal & h_reversal
        rows.append(
            {
                "entity_id": entity_id,
                "pert_iname": group["pert_iname"].iloc[0],
                "dataset": dataset,
                "scope": scope,
                "n_reversal_plastic": len(p_reversal),
                "n_reversal_hydrogel": len(h_reversal),
                "n_reversal_both": len(shared),
                "n_plastic_only_reversal": len(p_reversal - h_reversal),
                "n_hydrogel_only_reversal": len(h_reversal - p_reversal),
                "reversal_jaccard": safe_jaccard(p_reversal, h_reversal),
                "reversal_overlap_coefficient": safe_overlap_coefficient(p_reversal, h_reversal),
                "n_direction_switch": int(group["direction_switch"].sum()),
                "plastic_reversal_genes": ";".join(sorted(p_reversal)),
                "hydrogel_reversal_genes": ";".join(sorted(h_reversal)),
                "shared_reversal_genes": ";".join(sorted(shared)),
            }
        )
    return pd.DataFrame(rows)


def aggregate_reversal_sets(reversal_sets: pd.DataFrame) -> pd.DataFrame:
    """Aggregate set-level comparisons across the two LINCS releases."""
    rows = []
    for (entity_id, scope), group in reversal_sets.groupby(["entity_id", "scope"], sort=True):
        jaccard = pd.to_numeric(group["reversal_jaccard"], errors="coerce")
        overlap = pd.to_numeric(group["reversal_overlap_coefficient"], errors="coerce")
        rows.append(
            {
                "entity_id": entity_id,
                "pert_iname": group["pert_iname"].iloc[0],
                "scope": scope,
                "n_release_rows": int(len(group)),
                "n_releases_with_jaccard_ge_0_5": int((jaccard >= 0.5).sum()),
                "median_n_reversal_plastic": float(group["n_reversal_plastic"].median()),
                "median_n_reversal_hydrogel": float(group["n_reversal_hydrogel"].median()),
                "median_n_reversal_both": float(group["n_reversal_both"].median()),
                "median_reversal_jaccard": float(jaccard.median()) if jaccard.notna().any() else np.nan,
                "median_reversal_overlap_coefficient": float(overlap.median()) if overlap.notna().any() else np.nan,
                "median_n_plastic_only_reversal": float(group["n_plastic_only_reversal"].median()),
                "median_n_hydrogel_only_reversal": float(group["n_hydrogel_only_reversal"].median()),
                "total_direction_switches": int(group["n_direction_switch"].sum()),
            }
        )
    return pd.DataFrame(rows)


def reversal_gene_frequency(gene_comparison: pd.DataFrame) -> pd.DataFrame:
    """Count recurrent reversal calls separately for each disease panel."""
    # compare_genes() has already restricted this table to compound perturbations.
    source = gene_comparison.copy()
    rows = []
    for panel, status_column in (("plastic", "direction_status_plastic"), ("hydrogel", "direction_status_hydrogel")):
        status = source[status_column].fillna("not_available").astype(str)
        available = status.ne("not_available")
        panel_source = source.loc[available, ["entity_id", "dataset", "scope", "gene"]].copy()
        panel_source["is_reversal"] = status.loc[available].eq("reversal").to_numpy()
        if panel_source.empty:
            continue
        frequency = (
            panel_source.groupby("gene", sort=True)["is_reversal"]
            .agg(n_drug_release_reversal="sum", n_drug_release_observed="size")
            .reset_index()
        )
        frequency["n_unique_drugs_reversal"] = (
            panel_source.loc[panel_source["is_reversal"]].groupby("gene")["entity_id"].nunique()
        ).reindex(frequency["gene"]).fillna(0).astype(int).to_numpy()
        frequency["n_unique_releases_reversal"] = (
            panel_source.loc[panel_source["is_reversal"]].groupby("gene")["dataset"].nunique()
        ).reindex(frequency["gene"]).fillna(0).astype(int).to_numpy()
        frequency["fraction_reversal_among_observed"] = (
            frequency["n_drug_release_reversal"] / frequency["n_drug_release_observed"]
        )
        frequency.insert(0, "panel", panel)
        rows.append(frequency)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def main() -> None:
    plastic = pd.read_csv(PLASTIC_SUMMARY, low_memory=False)
    hydrogel = pd.read_csv(HYDROGEL_SUMMARY, low_memory=False)
    gene_comparison, drug_comparison = compare_genes(plastic, hydrogel)
    reversal_sets = reversal_set_comparison(gene_comparison)
    drug_reversal_summary = aggregate_reversal_sets(reversal_sets)
    reversal_frequency = reversal_gene_frequency(gene_comparison)

    signatures = pd.read_csv(SIGNATURES)
    p_panel = top_panel(signatures, "tsc2_loss_plastic")
    h_panel = top_panel(signatures, "tsc2_loss_hydrogel")
    panel_rows = []
    for direction in ("up", "down"):
        p_genes, h_genes = p_panel[direction], h_panel[direction]
        panel_rows.append(
            {
                "direction": direction,
                "plastic_genes": len(p_genes),
                "hydrogel_genes": len(h_genes),
                "same_direction_overlap": len(p_genes & h_genes),
                "opposite_direction_overlap": len(p_genes & h_panel["down" if direction == "up" else "up"]),
                "jaccard_same_direction": safe_jaccard(p_genes, h_genes),
            }
        )
    p_all, h_all = p_panel["up"] | p_panel["down"], h_panel["up"] | h_panel["down"]
    panel_rows.append(
        {
            "direction": "all",
            "plastic_genes": len(p_all),
            "hydrogel_genes": len(h_all),
            "same_direction_overlap": len((p_panel["up"] & h_panel["up"]) | (p_panel["down"] & h_panel["down"])),
            "opposite_direction_overlap": len((p_panel["up"] & h_panel["down"]) | (p_panel["down"] & h_panel["up"])),
            "jaccard_same_direction": safe_jaccard((p_panel["up"] & h_panel["up"]) | (p_panel["down"] & h_panel["down"]), p_all | h_all),
        }
    )
    panel_comparison = pd.DataFrame(panel_rows)

    PROGRAMS.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)
    gene_path = PROGRAMS / f"{PREFIX}_plastic_vs_hydrogel_gene_comparison.csv.gz"
    drug_path = PROGRAMS / f"{PREFIX}_plastic_vs_hydrogel_drug_program_comparison.csv"
    reversal_set_path = PROGRAMS / f"{PREFIX}_plastic_vs_hydrogel_drug_reversal_set_comparison.csv"
    reversal_summary_path = PROGRAMS / f"{PREFIX}_plastic_vs_hydrogel_drug_reversal_summary.csv"
    reversal_frequency_path = PROGRAMS / f"{PREFIX}_plastic_vs_hydrogel_panel_reversal_gene_frequency.csv"
    panel_path = PROGRAMS / f"{PREFIX}_plastic_vs_hydrogel_signature_panel_comparison.csv"
    gene_comparison.to_csv(gene_path, index=False, compression="gzip")
    drug_comparison.to_csv(drug_path, index=False)
    reversal_sets.to_csv(reversal_set_path, index=False)
    drug_reversal_summary.to_csv(reversal_summary_path, index=False)
    reversal_frequency.to_csv(reversal_frequency_path, index=False)
    panel_comparison.to_csv(panel_path, index=False)

    report_path = REPORTS / f"{PREFIX}_plastic_vs_hydrogel_comparison.md"
    overall = {
        "n_drug_comparison_rows": int(len(drug_comparison)),
        "n_unique_drugs": int(drug_comparison["entity_id"].nunique()),
        "n_gene_comparison_rows": int(len(gene_comparison)),
        "n_datasets": int(drug_comparison["dataset"].nunique()),
        "median_status_agreement": float(drug_comparison["status_agreement"].median()),
        "median_status_agreement_both_panels": float(drug_comparison["status_agreement_both_panels"].median()),
        "median_reversal_jaccard_release_level": float(reversal_sets["reversal_jaccard"].median()),
        "median_reversal_overlap_coefficient_release_level": float(reversal_sets["reversal_overlap_coefficient"].median()),
        "median_n_reversal_both_release_level": float(reversal_sets["n_reversal_both"].median()),
        "median_effect_spearman_secondary": float(drug_comparison["effect_spearman"].median()),
        "n_direction_switch_rows": int(drug_comparison["n_direction_switch"].sum()),
    }
    top = drug_reversal_summary.sort_values(
        ["median_n_reversal_both", "median_reversal_jaccard", "median_reversal_overlap_coefficient"],
        ascending=[False, False, False],
    ).head(20)
    top_columns = [
        "entity_id",
        "pert_iname",
        "n_release_rows",
        "n_releases_with_jaccard_ge_0_5",
        "median_n_reversal_plastic",
        "median_n_reversal_hydrogel",
        "median_n_reversal_both",
        "median_reversal_jaccard",
        "median_reversal_overlap_coefficient",
        "median_n_plastic_only_reversal",
        "median_n_hydrogel_only_reversal",
    ]
    top = top[top_columns]
    release_top = reversal_sets.sort_values(
        ["n_reversal_both", "reversal_jaccard", "reversal_overlap_coefficient"],
        ascending=[False, False, False],
    ).head(12)
    release_top = release_top[
        [
            "entity_id",
            "pert_iname",
            "dataset",
            "n_reversal_plastic",
            "n_reversal_hydrogel",
            "n_reversal_both",
            "n_plastic_only_reversal",
            "n_hydrogel_only_reversal",
            "reversal_jaccard",
            "reversal_overlap_coefficient",
            "shared_reversal_genes",
        ]
    ]
    common_effects_by_dataset = (
        drug_comparison.groupby("dataset")["n_effects_compared"].median().round(1).astype(float).to_dict()
    )
    frequency_columns = [
        "panel",
        "gene",
        "n_drug_release_reversal",
        "n_drug_release_observed",
        "n_unique_drugs_reversal",
        "n_unique_releases_reversal",
        "fraction_reversal_among_observed",
    ]
    top_frequency = reversal_frequency.sort_values(
        ["n_unique_drugs_reversal", "fraction_reversal_among_observed", "n_drug_release_reversal"],
        ascending=[False, False, False],
    ).head(30)[frequency_columns]
    plastic_frequency = reversal_frequency.loc[reversal_frequency["panel"].eq("plastic")].copy()
    hydrogel_frequency = reversal_frequency.loc[reversal_frequency["panel"].eq("hydrogel")].copy()
    shared_frequency = plastic_frequency.merge(hydrogel_frequency, on="gene", suffixes=("_plastic", "_hydrogel"))
    shared_frequency["min_fraction_reversal"] = shared_frequency[
        ["fraction_reversal_among_observed_plastic", "fraction_reversal_among_observed_hydrogel"]
    ].min(axis=1)
    top_shared_frequency = shared_frequency.sort_values(
        ["min_fraction_reversal", "n_unique_drugs_reversal_plastic", "n_unique_drugs_reversal_hydrogel"],
        ascending=[False, False, False],
    ).head(30)
    top_shared_frequency = top_shared_frequency[
        [
            "gene",
            "n_unique_drugs_reversal_plastic",
            "n_unique_drugs_reversal_hydrogel",
            "fraction_reversal_among_observed_plastic",
            "fraction_reversal_among_observed_hydrogel",
            "min_fraction_reversal",
        ]
    ]
    lines = [
        "# Plastic vs hydrogel LINCS signature comparison",
        "",
        "同一批 66 个候选药物分别使用 `tsc2_loss_plastic` 与 `tsc2_loss_hydrogel` top150 up + top150 down 面板；两个 LINCS release 先独立汇总，再比较共同基因。",
        "这里的主要问题不是重新比较同一套 LINCS perturbation effect，而是比较两个 disease panel 各自判定出的 reversal 基因集合：哪些基因在两个 panel 都被 reversal、哪些只在 plastic 或 hydrogel panel 中出现。",
        "两边使用同一批 LINCS drug signatures，因此共同可比较基因上的 drug effect 数值本来就相同；effect correlation 只作为结构性可比性诊断，不作为生物学证据。",
        "",
        f"- plastic panel: {len(p_all)} genes; hydrogel panel: {len(h_all)} genes",
        f"- same-direction panel overlap: {len((p_panel['up'] & h_panel['up']) | (p_panel['down'] & h_panel['down']))}; panel union: {len(p_all | h_all)}",
        f"- median common drug-effect genes compared per release: {common_effects_by_dataset}",
        "",
        "## Primary: reversal-set comparison",
        "",
        "对每个 drug × release 定义 plastic reversal set 与 hydrogel reversal set，并报告集合大小、交集、Jaccard、overlap coefficient 及两侧特异基因。跨 release 的主表使用这些集合指标的中位数；不把同一 LINCS effect 的相关性当成复现。",
        "",
        pd.DataFrame([overall]).to_string(index=False),
        "",
        "### Drugs with the largest shared reversal sets",
        "",
        top.to_string(index=False),
        "",
        "### Release-level shared reversal genes",
        "",
        release_top.to_string(index=False),
        "",
        "这些共享基因和各自特异基因的完整列表见 reversal-set comparison 输出；它们比共同基因上的 effect correlation 更直接回答两个 disease panel 是否把同一药物映射到相同的 reversal program。",
        "",
        "### Recurrent reversal genes by panel",
        "",
        "下面按 drug × release 统计每个 disease panel 中反复被判为 reversal 的基因；这是 panel-level reversal 频率，不是同一 LINCS effect 的重复测量。完整频率表同时保留 plastic、hydrogel、观察次数和药物数。",
        "",
        top_frequency.to_string(index=False),
        "",
        "#### Genes recurrent in both panels",
        "",
        top_shared_frequency.to_string(index=False),
        "",
        "## Secondary diagnostics",
        "",
        "status agreement、共同基因上的 effect Spearman 和 effect sign concordance 仅用于确认数据拼接/可比性；由于 drug effect 来自同一套 LINCS perturbation，不能把它们解释为独立生物学复现。",
        "",
        "解释边界：reversal 集合重叠仍不等于 hydrogel-specific causal mechanism。候选是否真正依赖 3D 环境，仍需结合 G×E / G×R×E、人体 niche 和选择性扰动验证。",
    ]
    report_path.write_text("\n".join(lines), encoding="utf-8")

    manifest = {
        "plastic_summary": str(PLASTIC_SUMMARY.relative_to(ROOT)),
        "hydrogel_summary": str(HYDROGEL_SUMMARY.relative_to(ROOT)),
        "signature_file": str(SIGNATURES.relative_to(ROOT)),
        "n_unique_drugs": int(drug_comparison["entity_id"].nunique()),
        "n_gene_comparison_rows": int(len(gene_comparison)),
        "panel_comparison": panel_comparison.to_dict("records"),
        "comparison_definitions": {
            "status_agreement": "same direction_status across the two panels",
            "non_neutral_direction_concordance": "same reversal/mimic direction among genes non-neutral in both panels",
            "reversal_jaccard": "Jaccard overlap of reversal gene sets",
            "reversal_overlap_coefficient": "intersection divided by the size of the smaller reversal gene set",
            "reversal_gene_frequency": "per-panel recurrence of reversal calls across drug × LINCS release rows",
            "effect_spearman": "Spearman correlation of drug_effect_median on common genes",
        },
        "outputs": {
            "gene_comparison": str(gene_path.relative_to(ROOT)),
            "drug_comparison": str(drug_path.relative_to(ROOT)),
            "reversal_set_comparison": str(reversal_set_path.relative_to(ROOT)),
            "drug_reversal_summary": str(reversal_summary_path.relative_to(ROOT)),
            "reversal_gene_frequency": str(reversal_frequency_path.relative_to(ROOT)),
            "panel_comparison": str(panel_path.relative_to(ROOT)),
            "report": str(report_path.relative_to(ROOT)),
        },
    }
    (ROOT / "manifests" / f"{PREFIX}_plastic_vs_hydrogel_comparison.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(overall, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
