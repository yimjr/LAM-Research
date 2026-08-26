"""Run Phase 3 exploratory state, microenvironment and hypothesis-card analyses."""

from __future__ import annotations

import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from scipy.stats import spearmanr
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score


ROOT = Path(__file__).resolve().parents[1]

STATE_PROGRAMS = {
    "contractile": ["ACTA2", "TAGLN", "MYL9", "CNN1", "MYH11", "TPM2"],
    "ecm_remodeling": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "SPARC"],
    "proliferative": ["MKI67", "TOP2A", "PCNA", "UBE2C", "CENPF"],
    "stress_hypoxia": ["HIF1A", "VEGFA", "DDIT3", "ATF4", "HSPA1A", "HSPA1B"],
    "inflammatory": ["IL6", "CXCL8", "CCL2", "TNF", "NFKBIA", "IFITM3"],
    "hormone_related": ["ESR1", "PGR", "GATA3", "FOXA1", "GREB1"],
    "metabolic": ["PPARG", "CPT1A", "ACACA", "SREBF1", "LDHA", "PDK1"],
    "mTOR_related": ["TSC2", "MTOR", "RPTOR", "RICTOR", "EIF4EBP1", "RPS6KB1", "AKT1", "MAPK1"],
}

MICROENVIRONMENT_PROGRAMS = {
    "lymphatic_endothelial": ["PDPN", "LYVE1", "CCL21", "FLT4", "PROX1", "CAV1"],
    "vascular_endothelial": ["PECAM1", "VWF", "KDR", "EMCN", "ESAM", "ENG"],
    "fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM", "COL3A1", "PDGFRA"],
    "immune": ["PTPRC", "LST1", "CD74", "CD3D", "NKG7", "MS4A1"],
}


def resolve_symbols(adata: ad.AnnData, requested: list[str]) -> tuple[list[str], list[str]]:
    raw_names = [str(name) for name in adata.raw.var_names]
    upper_lookup: dict[str, str] = {}
    if "gene_symbol_upper" in adata.raw.var:
        for actual, upper in zip(raw_names, adata.raw.var["gene_symbol_upper"].astype(str)):
            upper_lookup.setdefault(upper.upper(), actual)
    available, missing = [], []
    for gene in requested:
        actual = gene if gene in raw_names else upper_lookup.get(gene.upper())
        if actual is None:
            missing.append(gene)
        elif actual not in available:
            available.append(actual)
    return available, missing


def score_programs(adata: ad.AnnData, programs: dict[str, list[str]], prefix: str = "program") -> dict[str, dict]:
    status = {}
    for name, genes in programs.items():
        available, missing = resolve_symbols(adata, genes)
        column = f"{prefix}_{name}"
        if len(available) >= 2:
            try:
                sc.tl.score_genes(adata, gene_list=available, score_name=column, use_raw=True, ctrl_size=50, random_state=0)
                state = "completed"
            except Exception as exc:
                adata.obs[column] = np.nan
                state = f"failed: {type(exc).__name__}: {exc}"
        else:
            adata.obs[column] = np.nan
            state = "insufficient_genes"
        status[name] = {"requested": genes, "available": available, "missing": missing, "status": state}
    return status


def external_score(path: Path, dataset: str, programs: dict[str, list[str]], output: Path) -> dict:
    if not path.exists():
        return {"dataset": dataset, "status": "not_available"}
    ext = ad.read_h5ad(path)
    if "counts" not in ext.layers:
        ext.layers["counts"] = ext.X.copy()
    ext.var["gene_symbol_upper"] = ext.var_names.astype(str).str.upper().values
    ext.X = ext.X.astype(np.float32)
    sc.pp.normalize_total(ext, target_sum=10000)
    sc.pp.log1p(ext)
    ext.raw = ext.copy()
    statuses = score_programs(ext, {**STATE_PROGRAMS, **MICROENVIRONMENT_PROGRAMS}, prefix="external")
    rows = []
    for donor_id, group in ext.obs.groupby("donor_id", observed=True):
        row = {"dataset": dataset, "donor_id": donor_id, "cells": int(len(group))}
        for name in [*STATE_PROGRAMS, *MICROENVIRONMENT_PROGRAMS]:
            column = f"external_{name}"
            row[f"mean_{name}"] = float(group[column].mean()) if column in group else np.nan
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    return {"dataset": dataset, "status": "completed", "cells": int(ext.n_obs), "donors": int(ext.obs["donor_id"].nunique()), "program_status": statuses}


def pseudobulk(adata: ad.AnnData, output: Path) -> None:
    counts = adata.layers["counts"]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    rows, labels = [], []
    for donor_id, idx in adata.obs.groupby("donor_id", observed=True).groups.items():
        donor_mask = adata.obs.index.isin(idx)
        for label, candidate_mask in [("candidate", adata.obs["lamcore_candidate_author_style"].astype(bool).to_numpy()), ("other", ~adata.obs["lamcore_candidate_author_style"].astype(bool).to_numpy())]:
            mask = donor_mask & candidate_mask
            if mask.any():
                rows.append(np.asarray(counts[mask].sum(axis=0)).ravel())
                labels.append(f"{donor_id}__{label}")
    if rows:
        pd.DataFrame(rows, index=labels, columns=adata.var_names).to_csv(output)


def write_card(card_id: str, title_zh: str, title_en: str, body_zh: list[str], body_en: list[str]) -> None:
    card_dir = ROOT / "results/hypothesis_cards"
    card_dir.mkdir(parents=True, exist_ok=True)
    (card_dir / f"{card_id}_zh.md").write_text("\n".join([f"# {title_zh}", "", *body_zh]) + "\n")
    (card_dir / f"{card_id}_en.md").write_text("\n".join([f"# {title_en}", "", *body_en]) + "\n")


def main() -> None:
    input_path = ROOT / "data/processed/reproduction_core/GSE135851_core_reproduction.h5ad"
    result_dir = ROOT / "results/discovery"
    table_dir = result_dir / "tables"
    result_dir.mkdir(parents=True, exist_ok=True)
    table_dir.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_path)
    statuses = score_programs(adata, STATE_PROGRAMS, prefix="state")
    micro_statuses = score_programs(adata, MICROENVIRONMENT_PROGRAMS, prefix="microenv")

    state_columns = [f"state_{name}" for name in STATE_PROGRAMS]
    candidate = adata.obs["lamcore_candidate_author_style"].astype(bool)
    state_rows = []
    for donor_id, group in adata.obs[adata.obs["condition"].astype(str).eq("LAM")].groupby("donor_id", observed=True):
        cand = group["lamcore_candidate_author_style"].astype(bool)
        for name, column in zip(STATE_PROGRAMS, state_columns):
            state_rows.append({
                "donor_id": donor_id,
                "assay": group["assay"].iloc[0],
                "state": name,
                "candidate_cells": int(cand.sum()),
                "other_cells": int((~cand).sum()),
                "candidate_mean": float(group.loc[cand, column].mean()) if cand.any() else np.nan,
                "other_mean": float(group.loc[~cand, column].mean()) if (~cand).any() else np.nan,
                "difference": float(group.loc[cand, column].mean() - group.loc[~cand, column].mean()) if cand.any() and (~cand).any() else np.nan,
            })
    state_table = pd.DataFrame(state_rows)
    state_table["candidate_higher"] = state_table["difference"] > 0
    state_table.to_csv(table_dir / "lamcore_state_programs_by_donor.csv", index=False)

    cell_columns = ["sample_id", "donor_id", "condition", "assay", "lamcore_candidate_author_style", *state_columns, *[f"microenv_{x}" for x in MICROENVIRONMENT_PROGRAMS]]
    adata.obs[cell_columns].to_csv(table_dir / "lamcore_candidate_state_scores.csv")
    pseudobulk(adata[adata.obs["condition"].astype(str).eq("LAM")].copy(), table_dir / "donor_pseudobulk_counts_author_style.csv")

    candidate_cells = adata.obs.loc[candidate, state_columns].dropna()
    heterogeneity = []
    if len(candidate_cells) >= 10:
        matrix = candidate_cells.to_numpy()
        for k in (2, 3, 4):
            if len(candidate_cells) <= k:
                continue
            labels = KMeans(n_clusters=k, random_state=20260822, n_init=20).fit_predict(matrix)
            heterogeneity.append({"k": k, "cells": int(len(candidate_cells)), "silhouette": float(silhouette_score(matrix, labels)), "donors": int(adata.obs.loc[candidate_cells.index, "donor_id"].nunique())})
    pd.DataFrame(heterogeneity).to_csv(table_dir / "lamcore_state_heterogeneity_exploration.csv", index=False)

    association_rows = []
    lam_obs = adata.obs[adata.obs["condition"].astype(str).eq("LAM")]
    for state in STATE_PROGRAMS:
        state_col = f"state_{state}"
        for env in MICROENVIRONMENT_PROGRAMS:
            env_col = f"microenv_{env}"
            usable = lam_obs[[state_col, env_col]].dropna()
            rho, pvalue = spearmanr(usable[state_col], usable[env_col]) if len(usable) >= 10 else (np.nan, np.nan)
            association_rows.append({"state": state, "microenvironment_program": env, "cells": int(len(usable)), "spearman_rho": float(rho), "p_value": float(pvalue), "interpretation": "candidate association only; not direct cell communication"})
    pd.DataFrame(association_rows).to_csv(table_dir / "candidate_microenvironment_associations.csv", index=False)

    donor_state = state_table.groupby("state", observed=True).agg(
        donors_tested=("donor_id", "nunique"),
        donors_candidate_higher=("candidate_higher", "sum"),
        median_difference=("difference", "median"),
    ).reset_index()
    donor_state["fraction_candidate_higher"] = donor_state["donors_candidate_higher"] / donor_state["donors_tested"]
    donor_state.to_csv(table_dir / "lamcore_state_stability_summary.csv", index=False)

    external = []
    external.append(external_score(ROOT / "data/processed/external/GSE122960_normal_lung.h5ad", "GSE122960_normal_lung", STATE_PROGRAMS, table_dir / "external_GSE122960_state_scores.csv"))
    external.append(external_score(ROOT / "data/processed/external/GSE118180_wildtype_uterus.h5ad", "GSE118180_wildtype_uterus", STATE_PROGRAMS, table_dir / "external_GSE118180_state_scores.csv"))
    (result_dir / "external_validation_status.json").write_text(json.dumps(external, indent=2))

    # These are research leads, not claims of new LAM subtypes.
    stable = donor_state.sort_values("fraction_candidate_higher", ascending=False)
    top_states = stable.loc[stable["fraction_candidate_higher"] >= 0.5, "state"].astype(str).tolist()
    card1_zh = [
        "## 分类：高可信候选线索（同一数据集，尚无独立 LAM donor 验证）",
        "",
        f"## 观察到了什么\n在作者风格 marker 候选细胞中，连续状态程序在多个 LAM donor 中呈现方向一致的候选信号：{', '.join(top_states) if top_states else '没有达到多数 donor 方向一致'}。具体差值见 `lamcore_state_programs_by_donor.csv`。",
        "",
        "## donor、细胞和 pathway\n单位是 LAM1–LAM4 donor；涉及 contractile、ECM remodeling、stress/hypoxia、inflammatory、hormone-related、metabolic 和 mTOR-related 程序。",
        "",
        "## 稳健性\n候选群本身来自已知 marker + 作者风格图聚类；Phase 2 已比较 doublet、QC、聚类种子/分辨率、777 module score、rank-based score、assay 分层和 leave-one-donor-out。由于同一批数据用于发现与评估，仍不能称独立验证。",
        "",
        "## 替代解释\n候选定义与部分 contractile/ECM marker 有重叠；组织处理、assay、细胞周期、应激和供体差异都可能贡献信号。",
        "",
        "## 下一步验证\n在独立 LAM donor、空间转录组、蛋白或 snATAC 数据中，先测试这些连续程序是否仍出现在相同 marker-defined LAMCORE-like 细胞中，再做 TSC2/mTOR、ECM 和淋巴管相关实验。",
        "",
        "## 新颖性/可信度/优先级\n新颖性：中；可信度：中；推荐优先级：高。该卡片是可供研究者继续验证的线索，不是新亚型结论。",
    ]
    card1_en = [
        "## Class: high-confidence candidate clue (same dataset; no independent LAM-donor validation yet)",
        "",
        f"## Observation\nContinuous state programs in the author-style marker candidates showed directional signals across multiple LAM donors: {', '.join(top_states) if top_states else 'no majority-donor directional signal'}. See `lamcore_state_programs_by_donor.csv` for the donor-level differences.",
        "",
        "## Donors, cells and pathways\nThe unit is donor across LAM1–LAM4; programs include contractile, ECM remodeling, stress/hypoxia, inflammatory, hormone-related, metabolic and mTOR-related states.",
        "",
        "## Robustness\nCandidates were defined by known markers plus the author-style graph. Phase 2 compared doublets, QC, clustering seeds/resolutions, the 777-gene module score, a rank-based score, assay strata and leave-one-donor-out analyses. Because discovery and assessment use the same cohort, this is not independent validation.",
        "",
        "## Alternatives\nThe candidate definition overlaps some contractile/ECM markers; tissue processing, assay, cell cycle, stress and donor biology may contribute.",
        "",
        "## Next validation\nTest the programs in independent LAM donors, spatial transcriptomics, protein or snATAC data within the same marker-defined LAMCORE-like cells, then assess TSC2/mTOR, ECM and lymphatic mechanisms experimentally.",
        "",
        "## Novelty / confidence / priority\nNovelty: medium; confidence: medium; priority: high. This is a research lead, not a new-subtype claim.",
    ]
    write_card("card_01_continuous_states", "LAMCORE 可能包含连续表达状态", "LAMCORE May Contain Continuous Expression States", card1_zh, card1_en)

    lam2 = state_table[state_table["donor_id"].astype(str).eq("LAM2")]
    lam2_diff = ", ".join(lam2.loc[lam2["difference"].notna()].sort_values("difference").head(3)["state"].astype(str).tolist()) if not lam2.empty else ""
    card2_zh = [
        "## 分类：探索性假说",
        "",
        f"## 观察到了什么\nLAM2 在本操作性规则下有 {int((adata.obs.loc[adata.obs['donor_id'].astype(str).eq('LAM2'), 'lamcore_candidate_author_style']).sum())} 个候选，明显少于 LAM1/3/4。其相对差异较突出的程序包括：{lam2_diff or '当前表中没有足够细胞'}。",
        "",
        "## 候选解释\n可能涉及 donor-specific biology、治疗状态、组织处理或较弱的 LAMCORE-like 状态；不能把差异简单归因于 snRNA，因为 LAM2 不是 snRNA。",
        "",
        "## 如何验证\n需要更多明确临床信息的 LAM donor，尤其是 sirolimus 用药和疾病阶段，并在相同 assay 条件下比较；同时检查 doublet、cell cycle、批次和 marker detection。",
        "",
        "## 新颖性/可信度/优先级\n新颖性：中；可信度：低到中；推荐优先级：中。",
    ]
    card2_en = [
        "## Class: exploratory hypothesis",
        "",
        f"## Observation\nLAM2 has {int((adata.obs.loc[adata.obs['donor_id'].astype(str).eq('LAM2'), 'lamcore_candidate_author_style']).sum())} candidates under the operational rule, fewer than LAM1/3/4. Its most distinct programs in this comparison include: {lam2_diff or 'not enough cells for a stable comparison'}.",
        "",
        "## Candidate explanations\nDonor-specific biology, treatment state, tissue handling or a weaker LAMCORE-like state may contribute. The difference cannot simply be attributed to snRNA because LAM2 is not snRNA.",
        "",
        "## Validation\nUse more LAM donors with clinical metadata, especially sirolimus exposure and disease stage, compare matched assays, and check doublets, cell cycle, batch and marker detection.",
        "",
        "## Novelty / confidence / priority\nNovelty: medium; confidence: low-to-medium; priority: medium.",
    ]
    write_card("card_02_lam2_difference", "LAM2 的 LAMCORE-like 信号较弱或不同", "LAM2 Has a Weaker or Different LAMCORE-like Signal", card2_zh, card2_en)

    card3_zh = [
        "## 分类：探索性假说",
        "",
        "## 观察到了什么\nLAMCORE-like 候选与淋巴管内皮、血管内皮、成纤维细胞和免疫程序之间存在可计算的表达关联。关联表只描述共变，不是 ligand–receptor 真实通信。",
        "",
        "## 如何验证\n用空间转录组或成像确认细胞邻近关系，再用蛋白/功能实验验证具体 ligand–receptor、ECM 或 cytokine/growth-factor 关系。",
        "",
        "## 新颖性/可信度/优先级\n新颖性：中；可信度：低到中；推荐优先级：中。",
    ]
    card3_en = [
        "## Class: exploratory hypothesis",
        "",
        "## Observation\nLAMCORE-like candidates show computable expression associations with lymphatic endothelial, vascular endothelial, fibroblast and immune programs. These are co-variation signals, not evidence of physical ligand–receptor communication.",
        "",
        "## Validation\nUse spatial transcriptomics or imaging to establish proximity, then protein or functional assays for specific ligand–receptor, ECM or cytokine/growth-factor relationships.",
        "",
        "## Novelty / confidence / priority\nNovelty: medium; confidence: low-to-medium; priority: medium.",
    ]
    write_card("card_03_microenvironment", "LAMCORE 与微环境程序的候选关联", "Candidate LAMCORE–Microenvironment Associations", card3_zh, card3_en)

    summary = {
        "state_program_status": statuses,
        "microenvironment_program_status": micro_statuses,
        "candidate_cells": int(candidate.sum()),
        "candidate_donors": sorted(adata.obs.loc[candidate, "donor_id"].astype(str).unique().tolist()),
        "heterogeneity_models": heterogeneity,
        "external_validation": external,
        "cards": ["card_01_continuous_states", "card_02_lam2_difference", "card_03_microenvironment"],
        "communication_interpretation": "candidate association only; not direct cell communication",
    }
    (result_dir / "discovery_summary.json").write_text(json.dumps(summary, indent=2))
    zh = [
        "# LAMCORE 新生物学探索（Phase 3）",
        "",
        "## 结论定位",
        "",
        "本阶段从已经建立的作者风格 marker 候选出发，分析连续表达程序、LAM2 差异和微环境关联。结果用于生成研究线索，不把算法结构直接称为新亚型或真实细胞通信。",
        "",
        f"候选细胞数：{int(candidate.sum())}；候选供体：{', '.join(sorted(adata.obs.loc[candidate, 'donor_id'].astype(str).unique()))}。",
        "",
        "## 已生成的主要结果",
        "",
        "- `lamcore_state_programs_by_donor.csv`：每个 donor 的状态程序差值；",
        "- `lamcore_state_heterogeneity_exploration.csv`：连续状态空间的探索性聚类；",
        "- `candidate_microenvironment_associations.csv`：候选表达关联，不等同真实通信；",
        "- `external_validation_status.json`：GSE122960 正常肺和 GSE118180 小鼠子宫的处理状态；",
        "- `results/hypothesis_cards/`：三张中英文研究线索卡。",
        "",
        "当前最值得继续追踪的是 LAMCORE 内部连续状态，其次是 LAM2 的弱/异质信号以及与淋巴管和 ECM 的候选关联。独立 LAM donor 验证仍未完成。",
    ]
    en = [
        "# New Biological Exploration of LAMCORE (Phase 3)",
        "",
        "## Interpretation",
        "",
        "This phase starts from the author-style marker candidates and explores continuous expression programs, the LAM2 difference and microenvironment associations. Results are research leads; algorithmic structures are not called new subtypes or physical cell communication.",
        "",
        f"Candidate cells: {int(candidate.sum())}; candidate donors: {', '.join(sorted(adata.obs.loc[candidate, 'donor_id'].astype(str).unique()))}.",
        "",
        "## Main outputs",
        "",
        "- `lamcore_state_programs_by_donor.csv`: donor-level state differences;",
        "- `lamcore_state_heterogeneity_exploration.csv`: exploratory clustering in continuous state space;",
        "- `candidate_microenvironment_associations.csv`: candidate expression associations, not direct communication;",
        "- `external_validation_status.json`: processing status for normal-lung GSE122960 and mouse-uterus GSE118180;",
        "- `results/hypothesis_cards/`: three bilingual research-lead cards.",
        "",
        "The highest-priority lead is continuous LAMCORE state variation, followed by the weaker/heterogeneous LAM2 signal and candidate lymphatic/ECM associations. Independent LAM-donor validation remains outstanding.",
    ]
    (result_dir / "LAM_discovery_report_zh.md").write_text("\n".join(zh) + "\n")
    (result_dir / "LAM_discovery_report_en.md").write_text("\n".join(en) + "\n")

    manifest_path = ROOT / "manifests/run_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("steps", [])
    manifest["steps"] = [step for step in manifest["steps"] if step.get("name") != "phase3_biological_discovery"]
    manifest["steps"].append({"name": "phase3_biological_discovery", "completed_at": pd.Timestamp.now(tz="UTC").isoformat(), "input": str(input_path.relative_to(ROOT)), "outputs": ["results/discovery", "results/hypothesis_cards"], "external_validation": "normal-lung specificity and auxiliary mouse-uterus processing; independent LAM donor validation pending"})
    manifest["status"] = "phase1_core_reproduction_phase2_phase3_completed_pending_external_lam_validation"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    adata.write_h5ad(input_path, compression="gzip")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
