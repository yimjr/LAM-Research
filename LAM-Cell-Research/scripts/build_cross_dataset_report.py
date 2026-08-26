#!/usr/bin/env python3
"""Write a cautious bilingual summary of cross-dataset program matching."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/program_discovery/cross_dataset"


def main() -> None:
    matches = pd.read_csv(OUT / "cross_dataset_program_matches.csv")
    meta = pd.read_csv(OUT / "cross_dataset_meta_programs.csv")
    threshold = 0.15
    strong = matches[matches["top_gene_jaccard"] >= threshold]
    different_patient = strong[strong["independence_note"] == "different_patient_sets"]
    same_patient = strong[strong["independence_note"] == "same_patient_overlap_present"]
    top = matches.sort_values("top_gene_jaccard", ascending=False).head(10)
    zh = [
        "# 跨数据集程序比较报告",
        "",
        "> 本报告比较不同运行中候选程序的 top-50 gene overlap，不把 overlap 直接当作生物学复现。",
        "",
        "## 当前结果",
        "",
        f"- 共比较 {len(matches)} 对程序；",
        f"- Jaccard ≥ {threshold:.2f} 的匹配：{len(strong)} 对；",
        f"- 来自不同 PatientID 集合的强匹配：{len(different_patient)} 对；",
        f"- 含同一 PatientID 的匹配：{len(same_patient)} 对，主要涉及 GSE217108 与 GSE302356 的 LAM32。",
        "",
        "当前没有发现达到该阈值、且来自不同 PatientID 集合的稳定 meta-program。这个结果只能说明当前候选定义、特征选择和 top-gene 匹配下没有强信号，不能证明不存在新的跨患者程序。",
        "",
        "## 最高重叠匹配",
        "",
        "| pool | 数据集 | 程序 | 数据集 | 程序 | Jaccard | PatientID 关系 |",
        "|---|---|---|---|---|---:|---|",
    ]
    for _, row in top.iterrows():
        zh.append(f"| {row['pool']} | {row['dataset_left']} | {row['program_left']} | {row['dataset_right']} | {row['program_right']} | {row['top_gene_jaccard']:.3f} | {row['independence_note']} |")
    zh += [
        "",
        "## 下一步",
        "",
        "1. 在 donor 内独立发现结果上进行更稳健的 meta-program matching，而不是只比较 pooled NMF；",
        "2. 对候选程序进行 rank-based score、已知状态解释比例和 leave-one-donor-out 验证；",
        "3. 使用 GSE217108 ATAC、GSE302356 ATAC/空间数据检查正交支持；",
        "4. 将 LAM32 的同患者跨 assay 重复与真正不同患者的复现分开报告。",
    ]
    en = [
        "# Cross-Dataset Program Comparison Report",
        "",
        "> This report compares top-50 gene overlap across discovery runs; overlap is not treated as biological replication by itself.",
        "",
        "## Current result",
        "",
        f"- Program pairs compared: {len(matches)};",
        f"- Matches with Jaccard ≥ {threshold:.2f}: {len(strong)};",
        f"- Strong matches from different PatientID sets: {len(different_patient)};",
        f"- Matches involving the same PatientID: {len(same_patient)}, mainly LAM32 across GSE217108 and GSE302356.",
        "",
        "No stable meta-program from different PatientID sets reached this threshold. This means that the current candidate definition, feature selection, and top-gene matching did not produce a strong signal; it does not prove that no cross-patient program exists.",
        "",
        "## Highest-overlap matches",
        "",
        "| pool | dataset | program | dataset | program | Jaccard | PatientID relation |",
        "|---|---|---|---|---|---:|---|",
    ]
    for _, row in top.iterrows():
        en.append(f"| {row['pool']} | {row['dataset_left']} | {row['program_left']} | {row['dataset_right']} | {row['program_right']} | {row['top_gene_jaccard']:.3f} | {row['independence_note']} |")
    en += [
        "",
        "## Next steps",
        "",
        "1. Match independently discovered donor programs, not only pooled NMF programs;",
        "2. Test rank-based scores, known-state explained variance, and leave-one-donor-out stability;",
        "3. Use GSE217108 ATAC and GSE302356 ATAC/spatial data for orthogonal support;",
        "4. Report same-patient cross-assay repetition separately from true different-patient replication.",
    ]
    (OUT / "cross_dataset_program_report_zh.md").write_text("\n".join(zh) + "\n")
    (OUT / "cross_dataset_program_report_en.md").write_text("\n".join(en) + "\n")


if __name__ == "__main__":
    main()
