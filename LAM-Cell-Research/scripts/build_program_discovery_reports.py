#!/usr/bin/env python3
"""Build bilingual reports and research-line cards from program-discovery outputs."""

from __future__ import annotations

import json
import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "program_discovery"
CARDS = ROOT / "results" / "hypothesis_cards" / "unknown_programs"


def read_csv(name: str) -> pd.DataFrame:
    path = OUT / name
    return pd.read_csv(path) if path.exists() and path.stat().st_size else pd.DataFrame()


def top_genes(genes: pd.DataFrame, pool: str, program: str, n: int = 12) -> str:
    frame = genes[(genes["pool"] == pool) & (genes["candidate_program"] == program)].sort_values("rank_position").head(n)
    return ", ".join(frame["gene"].astype(str)) if not frame.empty else "无可用基因 / unavailable"


def strongest_known(known: pd.DataFrame, pool: str, program: str) -> str:
    frame = known[(known["pool"] == pool) & (known["candidate_program"] == program)].copy()
    if frame.empty:
        return "尚未匹配到已知程序 / no mapped known program"
    frame["abs_corr"] = frame["score_correlation"].abs().fillna(0)
    frame = frame.sort_values(["jaccard_top_genes", "abs_corr"], ascending=False).iloc[0]
    return f"{frame['known_program']} ({frame['evidence_scope']}; overlap={int(frame['overlap_n'])})"


def card(pool: str, program: str, genes: pd.DataFrame, known: pd.DataFrame, matches: pd.DataFrame, external: pd.DataFrame, language: str) -> str:
    g = top_genes(genes, pool, program)
    k = strongest_known(known, pool, program)
    m = matches[(matches["pool"] == pool) & (matches["best_pooled_program"] == program)] if not matches.empty else pd.DataFrame()
    donors = ", ".join(sorted(m.loc[m["independently_discovered"], "donor_id"].astype(str).unique())) if not m.empty else "none"
    ext_ready = ", ".join(external.loc[external["available"], "accession"].astype(str)) if not external.empty and external["available"].any() else "none"
    if language == "zh":
        return f"""# LAM Research Hypothesis Card：{pool} / {program}

## 当前分级

探索性假说（仅来自 GSE135851 同一队列；尚未达到独立验证标准）。

## 观察到了什么

在 `{pool}` 候选池的 pooled NMF 中出现了 `{program}`。当前排名靠前的基因为：{g}。

## donor 重复性

逐 donor 独立程序发现后，达到当前基因重叠标准的 donor：{donors}。这不是“在 pooled 模型中打分”，而是独立 donor 发现的初步匹配；若为空，说明目前没有足够证据称其跨 donor 重复。

## 与已知框架的关系

最强的已知程序匹配：{k}。部分重叠不会自动淘汰该候选；还需要判断它是否是已知状态在 LAMCORE 中的新实现、是否有 LAM-specific 基因或新的 TF/regulon。

## 外部证据

当前可用的外部 AnnData：{ext_ready}。ATAC、空间或蛋白证据尚未用于本卡片的结论。

## 替代解释

- donor、assay 或批次特异信号；
- cell cycle、doublet、测序深度或低质量造成的程序；
- 现有 CORE/SLS/IS/ECM 等状态的部分投影；
- 宽松候选集中仍未确认的细胞身份。

## 下一步验证

在 GSE190260、GSE217108 和 GSE302356 中按 PatientID 重复发现；按 donor 独立提取程序；比较已知程序解释比例；再检查 ATAC、空间或蛋白支持。未经这些步骤，不将其命名为新 LAMCORE 亚型或机制。

## 新颖性 / 可信度 / 优先级

- 新颖性：未评估；
- 当前可信度：低到中等，仅为同队列候选；
- 推荐优先级：中，取决于外部 donor 是否独立重现。
"""
    return f"""# LAM Research Hypothesis Card: {pool} / {program}

## Current tier

Exploratory hypothesis only. It comes from the GSE135851 cohort and has not met the independent-validation standard.

## Observation

`{program}` emerged in pooled NMF of the `{pool}` candidate pool. Leading genes are: {g}.

## Donor reproducibility

Donors meeting the current top-gene matching rule after donor-wise discovery: {donors}. This is independent donor discovery evidence, not passive scoring of a pooled model; an empty value means that cross-donor support is currently insufficient.

## Relation to known frameworks

Strongest known-program match: {k}. Partial overlap does not automatically invalidate a candidate; test whether it is a new LAMCORE implementation of a known state, has a LAM-specific component, or adds a new TF/regulon.

## External evidence

External AnnData currently available: {ext_ready}. ATAC, spatial, and protein evidence have not yet been used for the claim.

## Alternative explanations

- donor, assay, or batch-specific signal;
- cell cycle, doublet, sequencing depth, or low-quality effects;
- a projection of an existing CORE/SLS/IS/ECM state;
- unresolved identity in the broad candidate pool.

## Next validation

Re-discover the program by PatientID in GSE190260, GSE217108, and GSE302356; compare known-program explained variance; then test ATAC, spatial, or protein support. Do not name a new LAMCORE subtype or mechanism before these checks.

## Novelty / confidence / priority

- Novelty: not assessed;
- Current confidence: low to moderate, same-cohort candidate only;
- Priority: medium, conditional on independent donor re-discovery.
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="results/program_discovery")
    parser.add_argument("--cards-dir", default="results/hypothesis_cards/unknown_programs")
    args = parser.parse_args()
    global OUT, CARDS
    OUT = ROOT / args.output_dir
    CARDS = ROOT / args.cards_dir
    manifest = json.loads((OUT / "program_discovery_run_manifest.json").read_text())
    dataset_label = Path(manifest.get("input", "unknown_dataset")).stem
    pools = read_csv("candidate_pool_labels.csv")
    summary = read_csv("pooled_nmf_summary.csv")
    genes = read_csv("pooled_program_genes.csv")
    known = read_csv("known_program_comparisons.csv")
    matches = read_csv("donor_meta_program_matches.csv")
    external = read_csv("external_data_status.csv")
    counts = {c: int(pools[c].sum()) for c in ["pool_high_confidence", "pool_broad_lam_like", "pool_unrestricted_lam"] if c in pools}
    selected = summary.sort_values("selection_score", ascending=False).drop_duplicates("pool") if not summary.empty else pd.DataFrame()
    external_ready = not external.empty and external["available"].all()
    fast_mode = manifest.get("run_mode") == "fast_smoke"
    external_note_zh = (
        ("外部 GSE190260、GSE217108 和 GSE302356 已转换为可分析 AnnData；本报告是快速 smoke run，不能替代完整参数和独立性分析。" if fast_mode else "外部 GSE190260、GSE217108 和 GSE302356 已转换为可分析 AnnData，并已完成完整参数运行；结果仍需正交模态验证。")
        if external_ready
        else "外部 GSE190260、GSE217108 和 GSE302356 尚未全部转换为可分析 AnnData。"
    )
    external_note_en = (
        ("External GSE190260, GSE217108, and GSE302356 are now available as analysis-ready AnnData; this report is a fast smoke run and does not replace the full parameterized analysis." if fast_mode else "External GSE190260, GSE217108, and GSE302356 are now available as analysis-ready AnnData and have completed the full configured run; orthogonal modality validation remains." )
        if external_ready
        else "External GSE190260, GSE217108, and GSE302356 are not all available as analysis-ready AnnData."
    )
    lines_zh = [
        f"# LAMCORE 未知状态程序发现报告：{dataset_label}",
        "",
        f"> 本报告记录的是 {dataset_label} 中的程序发现候选，不是独立验证完成后的新生物学结论。",
        "",
        f"运行模式：`{manifest.get('run_mode', 'unknown')}`。" + ("当前参数用于流程验收，不作为最终模块数和稳定性结论。" if fast_mode else "当前使用配置中的完整 rank、seed 和 donor-wise 分析参数。"),
        "",
        "## 当前完成内容",
        "",
        f"- 输入：{manifest['n_cells']} 个细胞、{manifest['n_genes']} 个基因；",
        f"- 候选池：高置信 {counts.get('pool_high_confidence', 0)}，宽松 {counts.get('pool_broad_lam_like', 0)}，unrestricted guardrail {counts.get('pool_unrestricted_lam', 0)}；",
        "- 同时运行 pooled NMF 与 donor-wise 独立 NMF，并用核心基因重叠建立初步 meta-program；",
        "- 已知程序只做事后比较，没有在主分析前回归掉；",
        "- CORE3 使用 identity、深度校正后的低活性和 translation enrichment 三部分评分。",
        "",
        "## 当前结果如何解释",
        "",
        f"目前所有程序仍是候选。{external_note_zh} " + ("即使外部数据已可运行，也不能仅凭 smoke run 计算最终独立验证等级。" if fast_mode else "即使完整 RNA 分析已完成，也仍需正交模态和 PatientID 独立性检查。"),
        "",
        "## 运行摘要",
        "",
    ]
    if not selected.empty:
        for _, row in selected.iterrows():
            lines_zh.append(f"- `{row['pool']}`：选择 rank={int(row['rank'])}、seed={int(row['seed'])}，factor stability={row['factor_stability']:.3f}。")
    lines_zh += [
        "",
        "## 下一步",
        "",
        "1. 下载并检查外部公开处理后矩阵；",
        "2. 按 PatientID 进行 donor-wise 独立发现和 meta-program matching；",
        "3. 对跨 donor 候选做 doublet、assay、深度、已知程序和 leave-one-donor-out 验证；",
        "4. 用 GSE217108 的 ATAC、GSE302356 的 ATAC/空间和蛋白数据提高证据等级；",
        "5. 只有在身份、独立 donor 和正交证据同时支持时，才升级为高可信研究线索。",
    ]
    lines_en = [
        f"# LAMCORE Unknown-State Program Discovery Report: {dataset_label}",
        "",
        f"> This report documents program-discovery candidates in {dataset_label}. It is not a completed independent-validation claim.",
        "",
        f"Run mode: `{manifest.get('run_mode', 'unknown')}`. " + ("Parameters are for pipeline acceptance and not final rank or stability claims." if fast_mode else "The configured full rank, seed, and donor-wise analysis parameters were used."),
        "",
        "## Completed in this stage",
        "",
        f"- Input: {manifest['n_cells']} cells and {manifest['n_genes']} genes;",
        f"- Pools: high-confidence {counts.get('pool_high_confidence', 0)}, broad {counts.get('pool_broad_lam_like', 0)}, unrestricted guardrail {counts.get('pool_unrestricted_lam', 0)};",
        "- Pooled NMF and donor-wise independent NMF were both run, followed by preliminary meta-program matching;",
        "- Known programs were compared post hoc and were not regressed out before primary discovery;",
        "- CORE3 was modeled as identity, depth-adjusted low activity, and translation enrichment.",
        "",
        "## Interpretation",
        "",
        f"All programs remain candidates. {external_note_en} " + ("Even with external data available, a smoke run alone cannot establish the final independent-validation tier." if fast_mode else "Even after the full RNA analysis, orthogonal modalities and PatientID independence checks remain necessary."),
        "",
        "## Run summary",
        "",
    ]
    if not selected.empty:
        for _, row in selected.iterrows():
            lines_en.append(f"- `{row['pool']}`: selected rank={int(row['rank'])}, seed={int(row['seed'])}, factor stability={row['factor_stability']:.3f}.")
    lines_en += [
        "",
        "## Next steps",
        "",
        "1. Download and inspect the public processed matrices;",
        "2. Perform PatientID-aware donor-wise discovery and meta-program matching;",
        "3. Test cross-donor candidates against doublet, assay, depth, known-program, and leave-one-donor-out sensitivities;",
        "4. Use GSE217108 ATAC, GSE302356 ATAC/spatial, and protein data for orthogonal evidence;",
        "5. Upgrade a candidate only when identity, independent donors, and orthogonal evidence support it together.",
    ]
    (OUT / "LAM_unknown_program_report_zh.md").write_text("\n".join(lines_zh) + "\n")
    (OUT / "LAM_unknown_program_report_en.md").write_text("\n".join(lines_en) + "\n")

    CARDS.mkdir(parents=True, exist_ok=True)
    if not genes.empty:
        for pool in sorted(genes["pool"].unique()):
            for program in sorted(genes.loc[genes["pool"] == pool, "candidate_program"].unique()):
                stem = f"{pool}_{program}"
                (CARDS / f"{stem}_zh.md").write_text(card(pool, program, genes, known, matches, external, "zh"))
                (CARDS / f"{stem}_en.md").write_text(card(pool, program, genes, known, matches, external, "en"))
    print(json.dumps({"reports": [str(OUT / "LAM_unknown_program_report_zh.md"), str(OUT / "LAM_unknown_program_report_en.md")], "cards_dir": str(CARDS)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
