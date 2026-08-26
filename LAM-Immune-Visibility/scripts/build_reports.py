"""Build concise bilingual reports and hypothesis card from generated tables."""

from __future__ import annotations

import os
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(_PROJECT_ROOT / ".cache" / "matplotlib"))

import matplotlib.pyplot as plt
import pandas as pd

from common import PROJECT_ROOT, ensure_output_path, project_relative, write_json


def read_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path) if path.exists() else pd.DataFrame()


def build_figures(ranking: pd.DataFrame, patient_summary: pd.DataFrame) -> None:
    figure_dir = PROJECT_ROOT / "results" / "figures"
    if not ranking.empty:
        plot = ranking.sort_values("patient_consistency_any_detected", ascending=False).head(10).iloc[::-1]
        fig, ax = plt.subplots(figsize=(8, 4.5))
        ax.barh(plot["gene"], plot["patient_consistency_any_detected"], color="#5b8ff9")
        ax.set_xlabel("Patients with any detected signal / assayed patients")
        ax.set_title("Candidate antigen-associated expression consistency")
        fig.tight_layout()
        fig.savefig(ensure_output_path(figure_dir / "candidate_antigen_consistency.png"), dpi=160)
        plt.close(fig)
    if not patient_summary.empty:
        subset = patient_summary[patient_summary["identity_pool"].eq("high_confidence")]
        pivot = subset.pivot_table(index="patient_id", columns="module", values="mean_score", aggfunc="mean")
        columns = [column for column in ["antigen_associated", "presentation_machinery", "immune_evasion"] if column in pivot]
        if columns:
            fig, ax = plt.subplots(figsize=(9, 4.5))
            pivot[columns].plot.bar(ax=ax, color=["#61dDAA", "#65789B", "#F6BD16"][: len(columns)])
            ax.set_ylabel("Mean module score")
            ax.set_title("LAMCORE immune-visibility modules by patient")
            ax.legend(title="Module", fontsize=8)
            fig.tight_layout()
            fig.savefig(ensure_output_path(figure_dir / "patient_visibility_modules.png"), dpi=160)
            plt.close(fig)


def main() -> None:
    ranking = read_csv(PROJECT_ROOT / "results" / "candidate_antigens" / "candidate_antigen_ranking.csv")
    patients = read_csv(PROJECT_ROOT / "results" / "patient_summaries" / "patient_module_summary.csv")
    associations = read_csv(PROJECT_ROOT / "results" / "immune_context" / "patient_state_visibility_associations.csv")
    immune = read_csv(PROJECT_ROOT / "results" / "immune_context" / "patient_lamcore_immune_associations.csv")
    retention = read_csv(PROJECT_ROOT / "results" / "sirolimus_bridge" / "visibility_cross_environment_retention.csv")
    build_figures(ranking, patients)

    top = ranking.sort_values(["patient_consistency_any_detected", "lam_minus_normal_lung_detected_fraction"], ascending=False).head(10) if not ranking.empty else ranking
    top_lines = [
        f"- {row.gene}: 患者任一检出比例={row.patient_consistency_any_detected:.3f}, LAM-正常肺检出比例差={row.lam_minus_normal_lung_detected_fraction:.3f}" 
        for row in top.itertuples()
        if pd.notna(row.patient_consistency_any_detected)
    ]
    zh = "\n".join([
        "# LAMCORE 免疫可见性首轮计算报告",
        "",
        "## 结论定位",
        "",
        "这是基于既有公开处理后矩阵的计算研究，不是临床结论，也不证明真实 HLA 肽呈递或细胞间通信。",
        "",
        "## 缺失值与表达状态规则",
        "",
        "非零原始 count 即为已检出。`detected_low` 表示已经检出但处于预先定义的低表达区间；低/高只描述已检出信号的表达区间。基因不在矩阵/panel 中标记为 `not_assayed`，不能当作阴性。",
        "",
        f"## 当前结果规模",
        "",
        f"已生成候选抗原排序 {len(ranking)} 个、患者级模块汇总 {len(patients)} 行、状态关联 {len(associations)} 行、免疫上下文关联 {len(immune)} 行。",
        "",
        "## 候选表达一致性",
        "",
        *top_lines,
        "",
        "## 解释边界",
        "",
        "候选抗原只代表抗原相关表达；presentation module 只代表呈递机器表达。当前数据没有 immunopeptidomics，因此不能据此确认具体 HLA 肽。单细胞零值只表示当前 assay 中未检出，不能解释为真实不表达。",
        "空间和免疫关联均为候选关联；同一患者的多模态不增加独立患者数。rapamycin retention 只提供扰动支持的持久性线索。",
        "",
        "## 下一步",
        "",
        "优先在具有 HLA 型别、免疫肽组或功能 T-cell recognition 数据的样本中验证 PMEL/gp100 阳性对照及排名靠前候选；若要补充正常组织或肽呈递数据，应在本目录中新增数据，不修改源项目。",
    ])
    en = "\n".join([
        "# First-pass LAMCORE Immune Visibility Report",
        "",
        "## Positioning",
        "",
        "This is a computational analysis of existing public processed matrices. It does not establish clinical efficacy, peptide presentation, or direct cellular communication.",
        "",
        "## Expression-state rule",
        "",
        "Any nonzero raw count is detected. `detected_low` means detected and located in a predefined low-expression interval; low/high describe the expression interval after detection. Genes absent from a matrix or panel are `not_assayed`, not negative.",
        "",
        f"## Current outputs",
        "",
        f"The run produced {len(ranking)} antigen-associated candidates, {len(patients)} patient-level module rows, {len(associations)} state-association rows, and {len(immune)} immune-context association rows.",
        "",
        "## Interpretation",
        "",
        "Antigen modules represent antigen-associated expression, while presentation modules represent machinery expression. Without immunopeptidomics, no specific HLA peptide is confirmed. A single-cell zero is reported as not detected in the current assay, not as true biological absence.",
        "Spatial and immune associations are candidate associations only; same-patient modalities do not increase independent donor counts. Rapamycin retention is perturbation-supported persistence evidence, not patient-level treatment causality.",
    ])
    card_zh = "\n".join([
        "# Hypothesis Card：LAMCORE 免疫可见性",
        "",
        "## 分类",
        "探索性计算假说，尚需独立患者和正交实验验证。",
        "",
        "## 观察对象",
        "抗原相关表达、呈递机器和免疫逃逸状态在 LAMCORE-like 细胞中的患者级差异。",
        "",
        "## 严格表达定义",
        "非零 count 是已检出；低表达是已检出后位于低表达区间；未检出和未测量必须分别报告。",
        "",
        "## 实验验证方向",
        "HLA 分型/免疫肽组、PMEL/gp100 特异 T 细胞识别、以及 presentation/evasion 状态的蛋白验证。",
    ])
    reports = PROJECT_ROOT / "reports"
    ensure_output_path(reports / "immune_visibility_report_zh.md").write_text(zh)
    ensure_output_path(reports / "immune_visibility_report_en.md").write_text(en)
    ensure_output_path(reports / "hypothesis_cards" / "immune_visibility_state_zh.md").write_text(card_zh)
    write_json(PROJECT_ROOT / "manifests" / "report_manifest.json", {
        "reports": [project_relative(reports / "immune_visibility_report_zh.md"), project_relative(reports / "immune_visibility_report_en.md")],
        "strict_expression_language": "nonzero count is detected; detected_low is detected and in a low-expression interval",
        "retention_rows": len(retention),
    })
    print(f"reports written under {reports}")


if __name__ == "__main__":
    main()
