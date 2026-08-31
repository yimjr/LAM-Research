#!/usr/bin/env python3
"""Build the Stage 24 read-only project material package.

This script inventories and summarizes existing artifacts.  It deliberately
does not train models, recluster cells, alter candidate assignments, or write
to any artifact outside the requested Stage 24 output directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import date
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = PROJECT_ROOT / "results" / "stage24_final"
DATASETS = ["GSE135851", "GSE190260", "GSE217108", "GSE302356"]


def read_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}
    return value if isinstance(value, dict) else {}


def read_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except (OSError, pd.errors.ParserError, UnicodeDecodeError):
        return pd.DataFrame()


def scalar(frame: pd.DataFrame, row: int, column: str, default: Any = "") -> Any:
    if frame.empty or column not in frame.columns or row >= len(frame):
        return default
    value = frame.iloc[row][column]
    return default if pd.isna(value) else value


def norm_state(value: Any) -> str:
    try:
        return str(int(float(value)))
    except (TypeError, ValueError):
        return str(value).replace("State_", "").replace("State ", "").strip()


def json_text(value: Any) -> str:
    if isinstance(value, (dict, list)):
        return json.dumps(value, ensure_ascii=False, sort_keys=True)
    return str(value)


def rel(path: Path) -> str:
    try:
        return path.resolve().relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return str(path)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text.rstrip() + "\n", encoding="utf-8")


def write_frame(path: Path, frame: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8")


def markdown_table(frame: pd.DataFrame) -> str:
    """Render a small dataframe without adding a tabulate dependency."""
    if frame.empty:
        return "(no rows)"
    columns = [str(column) for column in frame.columns]

    def clean(value: Any) -> str:
        if pd.isna(value):
            return ""
        return str(value).replace("|", "\\|").replace("\n", " ")

    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join("---" for _ in columns) + " |",
    ]
    for values in frame.itertuples(index=False, name=None):
        lines.append("| " + " | ".join(clean(value) for value in values) + " |")
    return "\n".join(lines)


def sha256_small(path: Path, max_bytes: int = 50_000_000) -> str:
    if not path.is_file() or path.stat().st_size > max_bytes:
        return "not_computed_large_or_directory"
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_shape(path: Path) -> tuple[str, str]:
    suffix = path.suffix.lower()
    if suffix == ".csv":
        try:
            with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
                reader = csv.reader(handle)
                header = next(reader, [])
                rows = sum(1 for _ in reader)
            return str(rows), str(len(header))
        except OSError:
            return "", ""
    if suffix == ".json":
        return "json", ""
    if suffix == ".md":
        return "markdown", ""
    if suffix == ".npz":
        return "npz", ""
    if suffix in {".h5ad", ".pt"}:
        return "binary_artifact", ""
    return "", ""


STAGE_SPECS: list[dict[str, Any]] = [
    {
        "stage": 1,
        "title": "输入盘点与上游继承清单",
        "question": "四套 AnnData、四套 candidate/program 结果和共享配置是否可追溯、可用？",
        "inputs": "LAM-Cell-Research / data-temp / data/upstream；config 与 manifests",
        "method": "只读检查路径、字段、counts layer、三个 pool 层级，并递归发现同类 upstream 结果。",
        "outputs": "results/stage1_6/input_inventory.{csv,json}; upstream_annotation_manifest.csv",
        "checkpoint": "四套 AnnData 和四套 candidate_pool_labels 均 ready；normal 为可选输入。",
        "impact": "确定直接继承转换后的 AnnData 和 upstream annotations，不重做 GEO 转换。",
        "later_revision": "无；后续阶段继续使用这份继承边界。",
    },
    {
        "stage": 2,
        "title": "继承与准备",
        "question": "能否把 upstream annotation、registry 和 raw counts 安全映射到统一分析表？",
        "inputs": "四套上游 AnnData、candidate_pool_labels、core3/program 结果、donor_registry.yaml",
        "method": "按 cell ID/sample key 合并，canonicalize FIGF→VEGFD，固定 high-confidence candidate 和 boundary，验证 counts。",
        "outputs": "data/interim/inherited；results/stage1_6/upstream_inheritance.json",
        "checkpoint": "lam_candidate 只等于 pool_high_confidence；unrestricted 仅审计。",
        "impact": "为 QC/harmonization 建立统一 prepared 输入，同时把旧标签留在 upstream 命名空间。",
        "later_revision": "Stage 16 的新 gate 只作诊断，不替换这里冻结的主线 candidate pool。",
    },
    {
        "stage": 3,
        "title": "QC 与 gene harmonization",
        "question": "在不重复核心 QC 的前提下，能否为外部数据补齐本项目 QC 并建立共同 gene universe？",
        "inputs": "继承后的四套 AnnData；core qc_pass；外部 raw counts",
        "method": "GSE135851 继承原 qc_pass；外部使用 min_genes=200、min_counts=500、scRNA mt<20%、snRNA mt<10%；保留 doublet。",
        "outputs": "state_model_prepared.h5ad；qc_summary.csv",
        "checkpoint": "外部 QC 后：GSE190260 38,701、GSE217108 6,941、GSE302356 21,771 cells pass。",
        "impact": "后续 NMF/scVI 共享 alias-corrected gene universe；GSE135851 不被二次过滤。",
        "later_revision": "Stage 15–17 发现的是 upstream candidate gate 问题，不是本阶段 QC 合并错误。",
    },
    {
        "stage": 4,
        "title": "PCA/NMF baseline",
        "question": "在 State Modeling 自己的输入矩阵上，能否建立不依赖旧 state 的 baseline？",
        "inputs": "state_model_prepared.h5ad；counts layer",
        "method": "counts→library-size normalization→log1p→HVG(4000)→PCA/UMAP/Leiden/NMF；NMF 5 components、top 2000 features、max_iter 400、最多 12000 cells。",
        "outputs": "state_model_baseline.h5ad；baseline_cluster_summary.csv；nmf_cell_scores.csv；nmf_top_genes.csv",
        "checkpoint": "获得独立 baseline；旧 program/state 只保留作 post-hoc 对照。",
        "impact": "明确 NMF 不使用 scaled/PCA/scVI latent，scVI 也不使用 NMF 矩阵。",
        "later_revision": "没有用 baseline cluster 数直接定义最终 20 states。",
    },
    {
        "stage": 5,
        "title": "scVI latent model",
        "question": "只用 dataset 做 batch correction 时，是否能构建可复用 latent space？",
        "inputs": "state_model_prepared.h5ad 的 layers[counts]",
        "method": "scVI layer=counts、batch_key=dataset、无 categorical covariate、n_latent=20、n_layers=2、n_hidden=128、max_epochs=200、early stopping；CUDA 优先。",
        "outputs": "state_model_scvi.h5ad；scvi_model/model.pt；scvi_training_manifest.json",
        "checkpoint": "55,238 cells 的 200-epoch 模型完成；assay 留在 metadata，不进入 covariate。",
        "impact": "Stage 6–23 均复用同一 X_scVI，不再重训。",
        "later_revision": "无；后续所有结构分析都明确排除重新训练。",
    },
    {
        "stage": 6,
        "title": "LAM-only latent structure checkpoint",
        "question": "high-confidence LAM candidate 内是否存在跨患者共享、且不由 patient/dataset/assay 主导的结构？",
        "inputs": "X_scVI；5,378 pool_high_confidence cells；boundary/normal 仅辅助",
        "method": "LAM-only kNN→Leiden；独立网格 n_neighbors 15/30/50×resolution 0.2/0.4/0.6；保留全体 33-cluster 标签作审计，不用于 LAM 内部判定。",
        "outputs": "stage6_parameter_grid.csv；stage6_checkpoint.json；grid/cluster/ARI 表",
        "checkpoint": "GO；12 个 LAM-only reference clusters、12 个跨至少 2 patients 的 qualified clusters；patient ARI 0.152554、dataset ARI 0.086786、assay ARI 0.009094。",
        "impact": "允许进入 consensus；也暴露出单一参数和固定 preprocess 参数不应混用。",
        "later_revision": "Stage 7 将 9 个 grid configuration 等权并形成 20 个 consensus states；12 与 20 是不同对象，不是覆盖。",
    },
    {
        "stage": 7,
        "title": "Consensus stability",
        "question": "哪些 LAM-only 结构跨 grid/seed 稳定？",
        "inputs": "同一 X_scVI 的 5,378 candidate cells；9 configurations、21 raw partitions",
        "method": "每个 configuration 内先平均 seed co-assignment，再让 9 configurations 等权；最终生成 5,378×5,378 float32 co-assignment 并做完整距离 average-linkage。",
        "outputs": "state_consensus_assignments.csv；state_stability_summary.csv；cluster_matching_across_grid.csv；coassignment_matrix.npz；consensus_dendrogram.npz",
        "checkpoint": "形成 20 个当前 consensus states；seed stability 单独记录。",
        "impact": "把 Stage 6 的单一参考 partition 转为可追踪的 20-state consensus。",
        "later_revision": "20-state labels 在 Stage 8–23 冻结；后续 hierarchy 只作解释。",
    },
    {
        "stage": 8,
        "title": "Leave-one-out robustness",
        "question": "去掉一个 patient/dataset 后，consensus state 是否仍可恢复？",
        "inputs": "Stage 7 consensus；X_scVI；5,378 candidate cells",
        "method": "full-data 30/0.4 reference；LOO 前将 reference、consensus 和 LOO cluster 都限制在 retained cells；允许 split/merge matching。",
        "outputs": "loo_runs.csv；loo_cluster_matches.csv；loo_state_summary.csv；full_reference_consensus_matches.csv",
        "checkpoint": "产生连续的 baseline Jaccard、LOO recovery 和 additional loss。",
        "impact": "为每个 state 的结构稳定性提供 patient/dataset 支持，而不是只看 full partition。",
        "later_revision": "Stage 11 将这些连续指标与患者级生物学复现合并。",
    },
    {
        "stage": 9,
        "title": "Hierarchy and continuum",
        "question": "20 states 是平级离散状态，还是存在 parent/substate/连续连接？",
        "inputs": "consensus labels；X_scVI；boundary transitions",
        "method": "state latent distance、connectivity/PAGA 等价汇总、split/merge tree、boundary transitions；不重新聚类。",
        "outputs": "state_distance_matrix.csv；state_connectivity.csv；state_split_merge_tree.csv；boundary_state_transitions.csv",
        "checkpoint": "为 atlas 提供描述性 parent/substate 结构。",
        "impact": "避免把 resolution 改变产生的层级现象硬写成 20 个互不相干状态。",
        "later_revision": "Stage 22 进一步把全局问题收窄到 State15 周围局部分支。",
    },
    {
        "stage": 10,
        "title": "Per-state biology and DE",
        "question": "每个 consensus state 的表达/程序差异是否有患者级证据？",
        "inputs": "state_model_prepared.h5ad full gene universe/raw counts；consensus state",
        "method": "每个 State_k 单独对同患者 Rest_of_LAM 做 patient×group pseudobulk；设计 ~ patient_id + group；至少 3 patients 才正式 DE；每 state 独立 FDR。",
        "outputs": "state_de_results.csv（159,059 rows）；state_markers.csv；state_pseudobulk_counts.csv；state_program_scores.csv；pathway/regulon placeholders",
        "checkpoint": "20 个独立 state-vs-rest 分析路径，不存在 all-state 多分类模型。",
        "impact": "提供 state-level markers/program evidence，同时承认 pathway/regulon 当前 unavailable。",
        "later_revision": "Stage 24 不重新做 DE，只整理其支持范围。",
    },
    {
        "stage": 11,
        "title": "Patient-level reproducibility",
        "question": "细胞层面的 state 是否能转化为跨患者 evidence？",
        "inputs": "Stage 8 LOO；Stage 10 DE/pseudobulk；patient metadata",
        "method": "patient×state matrix；cells/fraction/signature/log2FC/direction/dataset coverage；保留 structural stability 与 biological reproducibility 两维。",
        "outputs": "patient_state_matrix.csv；state_reproducibility_summary.csv",
        "checkpoint": "State15 structural 0.854111、biological 0.374758；State16 0.719633/0.323088；State18 0.940033/0.408960；State20 0.912545/0.420718。",
        "impact": "把“cluster 很稳定”与“跨患者 biology 可复现”分开。",
        "later_revision": "Stage 19 对 State15 的患者组成和最大患者敏感性作专门审计。",
    },
    {
        "stage": 12,
        "title": "Boundary and normal auxiliary validation",
        "question": "候选边界和正常参考能否帮助解释 state 的外部邻域？",
        "inputs": "X_scVI；boundary 16,883；normal 32,977",
        "method": "boundary/normal 邻域与距离；不加入 state 数量，不改变 candidate/state 定义。",
        "outputs": "boundary_validation.csv；normal_validation.csv；state_auxiliary_summary.csv",
        "checkpoint": "normal 与 boundary 仅为辅助；State15 normal mean distance 3.981922，State18 3.023263，State20 3.240648。",
        "impact": "为后续 State15 anchor 与 manifold 分析提供参照。",
        "later_revision": "Stage 20–22 保持 normal remote、boundary projection 的辅助地位。",
    },
    {
        "stage": 13,
        "title": "State atlas and hypotheses",
        "question": "如何用连续证据汇总 20 states，而不把单一分数当 confidence？",
        "inputs": "Stages 7–12 的 state、DE、program、LOO、auxiliary 结果",
        "method": "汇总 structural/biological/coverage/normal/boundary/upstream correspondence；不设硬性 high/medium/low 门槛。",
        "outputs": "state_atlas.csv/json；state_hypothesis_candidates.csv；state_atlas.h5ad",
        "checkpoint": "形成第一版 20-state atlas 与 10 个 hypothesis candidates。",
        "impact": "为 Stage15–22 提供固定状态编号和证据背景。",
        "later_revision": "Stage 24 重新解释这些 labels 的生物学含义，但不改写 atlas artifact。",
    },
    {
        "stage": 14,
        "title": "Consensus/upstream merge",
        "question": "如何把 frozen consensus 与 upstream annotation 放在同一逐细胞表？",
        "inputs": "Stage 7 consensus；upstream candidate/state/program annotations",
        "method": "按 cell ID 一对一合并，保留 upstream 命名空间；不让旧 labels 进入 scVI/NMF。",
        "outputs": "state_consensus_with_upstream_annotations.csv；merge manifest",
        "checkpoint": "5,378 cells 的合并字段可追溯；Stage 15 audit later confirmed merge inconsistencies=0。",
        "impact": "支持 Stage 15 candidate identity audit 和后续 post-hoc comparison。",
        "later_revision": "外部 candidate/state 结果的完整继承在 Stage 1–6 manifest 中保留。",
    },
    {
        "stage": 15,
        "title": "Candidate identity audit",
        "question": "5,378 candidate 是否主要由过宽的 marker-combo gate 产生？",
        "inputs": "原 candidate_pool_labels；merged consensus；原始 marker_expr/counts",
        "method": "ID/字段一对一核对；重算原规则；FIGF/VEGFD duplicate audit；组合、UMI 和 state-level diagnostics；只读。",
        "outputs": "annotation_merge_audit.csv；rule_recalculation_audit.csv；marker_patterns.csv；root_cause_evidence.csv",
        "checkpoint": "A merge error=0；FIGF/VEGFD duplicate-pass=0；5,238 marker-combo、140 author/formal；1,443 个 marker-combo candidate 仅有 1-UMI 支持。",
        "impact": "根因定位到 C：任意两个 marker>0 gate 特异性不足，而非 Stage Modeling 搬错。",
        "later_revision": "Stage 16 提出连续 identity gate，但没有写回主线 candidate pool。",
    },
    {
        "stage": 16,
        "title": "Identity gate reconstruction",
        "question": "能否用 identity anchors+support+competing lineage 重建更可解释的 gate？",
        "inputs": "所有 condition=LAM cells；连续 module scores；独立正/负参照",
        "method": "PMEL/MLANA/MITF+LAMCORE/CORE evidence 为 identity；ACTA2/ESR1/VEGFD/CTSK 为 support；竞争谱系为 penalty；不使用旧 state 调参。",
        "outputs": "cell_identity_evidence.csv；reference_calibration.csv；LODO；new_candidate_assignment.csv",
        "checkpoint": "LAM_core 208、boundary 65,930、non_LAM_like 24,503；这是独立诊断 gate，不替换 frozen candidate。",
        "impact": "使 Stage 17 能够专门追踪外部 positive reference 的漏检。",
        "later_revision": "Stage 17 发现 GSE190260 的 dataset calibration/dropout/penalty 问题；Stage 24 不继续调 gate。",
    },
    {
        "stage": 17,
        "title": "Cross-dataset identity calibration audit",
        "question": "GSE190260 为什么漏掉 upstream CORE3-positive？",
        "inputs": "Stage 16 artifacts；formal 777-gene signature；raw counts",
        "method": "正参考逐细胞 score decomposition、dropout/depth、competing penalty、raw/z/percentile counterfactual、LODO；不生成新 assignment。",
        "outputs": "positive_reference_failures.csv；component_scores_by_dataset.csv；counterfactual_calibration.csv；root_cause_by_dataset.csv",
        "checkpoint": "GSE190260 2,117 positives：core recovery 0、core+boundary 0.708077；median final score -1.006577、score shift -6.610707、median penalty 4.950320；primary category competing_penalty_only，且 identity/support dropout 均高。",
        "impact": "证明跨数据集 score scale/marker dropout/penalty 联合作用，不能把外部 0 author-style 当 negative。",
        "later_revision": "Stage 18 formal LAMCORE anchor validation独立使用 777 genes；不改 Stage16。",
    },
    {
        "stage": 18,
        "title": "State15 anchor validation",
        "question": "冻结的 State15 能否作为 LAM-core reference anchor？",
        "inputs": "State15 200 cells；777-gene formal LAMCORE；author labels；normal/boundary/comparators",
        "method": "LAMCORE score、author enrichment/Fisher、marker/program、patient pseudobulk、comparator 和 latent-neighborhood validation。",
        "outputs": "state15_anchor_summary.json；state15_lamcore_summary.csv；state15_vs_comparators.csv；anchor report",
        "checkpoint": "State15 LAMCORE median 0.5125；normal 0.0513；overall author enrichment 13.9408、Fisher p=2.229e-36；decision provisional_reference_candidate_not_formally_upgraded。",
        "impact": "State15 被冻结为后续 reference anchor candidate，而不是正式 classifier。",
        "later_revision": "Stage19 量化 LAM1163 enrichment；最终仍保持 provisional。",
    },
    {
        "stage": 19,
        "title": "State15 cross-patient audit",
        "question": "State15 是否只是 LAM1163 这个大患者造成的？",
        "inputs": "frozen State15；candidate pool composition；author availability；patient pseudobulk/LOPO",
        "method": "patient composition baseline、patient-matched comparison、7-patient LOPO、remove-LAM1163 sensitivity；author unavailable 标成 not_assayed。",
        "outputs": "state15_patient_composition.csv；author_annotation_availability.csv；LOPO；state15_without_LAM1163.csv",
        "checkpoint": "LAM1163 占 candidate 9.2971%、State15 63.5%，enrichment 6.8301；去除后保留 73 cells、LAMCORE median 0.4683；author labels 仅 GSE135851 可检验。",
        "impact": "把 State15 定位为患者富集但 profile 可保留的中间结果。",
        "later_revision": "为 Stage20 的 anchor-centered geometry 提供组成警示。",
    },
    {
        "stage": 20,
        "title": "State15-centered latent geometry",
        "question": "以 State15 为中心，周围细胞的 identity/program 是否沿 latent distance 变化？",
        "inputs": "X_scVI；State15 200；candidate 5,378+boundary 16,883；normal remote",
        "method": "latent-only distance、k=30 neighbors、distance bins、identity/lineage mapping、State16 co-expression、patient/dataset gradients。",
        "outputs": "state15_centered_manifold.csv；distance bins/gradients；State16 audits；stage20 report",
        "checkpoint": "初始 checkpoint supports_lam_centered_transcriptional_manifold；full LAMCORE 近端→远端 0.2152→0.1556，4/4 dataset rho<0，8/12 patient rho<0。",
        "impact": "提出 State15-centered manifold 假设，并把验证对象从单个 state 扩展到邻域。",
        "later_revision": "Stage21 明确削弱 pooled/global manifold；Stage22 改为局部分支 candidate。",
    },
    {
        "stage": 21,
        "title": "Independent State15 manifold validation",
        "question": "去掉 anchor、自身 gate/scVI feature overlap 和 composition 后，gradient 是否仍成立？",
        "inputs": "Stage20 distance；非-State15 22,061；777-gene LAMCORE；candidate-only null pool",
        "method": "LAMCORE full/no_gate/outside_scVI/independent；patient-adjusted regression、dataset/patient replication、500 matched fake anchors、State16/boundary auxiliary。",
        "outputs": "gradient_models.csv；dataset/patient gradients；matched_anchor_null.csv；stage21 report",
        "checkpoint": "554 genes independent；candidate-only independent slope -0.015466，null empirical two-sided p=0.001996；但 pooled rank rho +0.075665、near/far medians 0.1083/0.1443，患者异质；checkpoint state15_lam_rich_gradient_but_not_robust_manifold。",
        "impact": "不再支持稳健统一 global manifold，但保留局部/方向性结构的可能。",
        "later_revision": "Stage22 将问题改成 local branch decomposition。",
    },
    {
        "stage": 22,
        "title": "State15 local branch decomposition",
        "question": "State15 周围哪些局部方向真正保留 LAM identity？",
        "inputs": "同一 X_scVI；22,261 candidate+boundary；Stage21 scores；frozen State labels",
        "method": "k=30 不重新 Leiden；1–3 hop；自动选择 ≥10 direct cells 且 ≥2 patients 的 external states；near/mid/far、patient consistency、500 matched null。",
        "outputs": "branch_candidates.csv；branch evidence/gradients；boundary assignments；stage22 report",
        "checkpoint": "历史版本：State16/12/20/7 入选；旧的 zero-centered matched-null p=0.025948 曾将 State16 标为 transition candidate。该版本已被 Stage22 修正分析 supersede。",
        "impact": "最终主线从“统一 manifold”收窄到 State15 周围的局部 branch/adjacency 审计；修正后不再升级 State16 为 transition candidate。",
        "later_revision": "Stage23 仅可视化，不再改变该 checkpoint。",
    },
    {
        "stage": 23,
        "title": "Latent-space visualization",
        "question": "如何把 Stage15–22 已有结构直观看清，而不新增结论？",
        "inputs": "X_scVI；Stage20/21/22 frozen tables；State15 anchor",
        "method": "全局/局部 2D static plots；3D UMAP、3D local graph、3D PCA；统一 hover/dropdown；不重训、不重聚类。",
        "outputs": "results/stage23_visualization/ 五组 PNG/PDF 与三个 HTML；visualization_manifest.json",
        "checkpoint": "可视化资产完成；2D 优先复用已有 UMAP，3D UMAP/PCA 只用于展示/对照。",
        "impact": "为 Stage24 artifact provenance 和最终报告提供可引用图件。",
        "later_revision": "无；Stage24 不把可视化解释升级为新证据。",
    },
]


def stage22_snapshot() -> dict[str, Any]:
    """Read the current corrected Stage22 outputs for downstream synthesis."""
    directory = PROJECT_ROOT / "results/stage22"
    manifest = read_json(directory / "stage22_manifest.json")
    branches = read_csv(directory / "branch_candidates.csv")
    evidence = read_csv(directory / "branch_evidence_summary.csv")
    boundary = read_csv(directory / "boundary_local_branch_assignment.csv")
    selected = sorted(branches.get("source_state", pd.Series(dtype=str)).dropna().astype(str).unique().tolist())
    selected_text = ", ".join(selected) if selected else "none"
    lines: list[str] = []
    for _, row in evidence.iterrows():
        lines.append(
            f"{row.get('source_state', '')}: slope={row.get('LAMCORE_independent_latent_slope', ''):.6f}, "
            f"raw_p={row.get('matched_null_empirical_p', ''):.6f}, "
            f"BH_q={row.get('matched_null_q_value', ''):.6f}, "
            f"label={row.get('evidence_label', '')}"
        )
    state16 = evidence[evidence.get("source_state", pd.Series(dtype=str)).astype(str).eq("State_16")]
    state16_row = state16.iloc[0].to_dict() if len(state16) else {}
    boundary_unresolved = int((boundary.get("branch_assignment", pd.Series(dtype=str)).astype(str).eq("unresolved")).sum()) if len(boundary) else 0
    snapshot = {
        "manifest": manifest,
        "selected_states": selected,
        "selected_text": selected_text,
        "evidence_lines": "; ".join(lines) if lines else "no branch evidence",
        "state16_slope": float(state16_row.get("LAMCORE_independent_latent_slope")) if pd.notna(state16_row.get("LAMCORE_independent_latent_slope")) else float("nan"),
        "state16_raw_p": float(state16_row.get("matched_null_empirical_p")) if pd.notna(state16_row.get("matched_null_empirical_p")) else float("nan"),
        "state16_q": float(state16_row.get("matched_null_q_value")) if pd.notna(state16_row.get("matched_null_q_value")) else float("nan"),
        "state16_label": str(state16_row.get("evidence_label", "not_available")),
        "boundary_rows": len(boundary),
        "boundary_unresolved": boundary_unresolved,
        "branch_count": len(branches),
    }
    return snapshot


STATE_INTERPRETATIONS: dict[str, dict[str, str]] = {
    "1": {"analogue": "ciliated/airway-like epithelial", "class": "relatively_clear_normal-lineage analogue", "support": "Top DE markers include DNAI1, DNAH9, DCDC1 and other cilia genes; strong ciliated geometry/identity signal in Stage16 diagnostic summary.", "conflict": "Upstream cell_type is unknown; Stage13 programs also contain hormone/HOX and interstitial signals.", "lam": "Not supported as a LAM-core state; likely a normal airway-lineage contaminant within the high-recall pool.", "uncertainty": "Moderate: cell-type label is inferred from expression, not an upstream verified annotation."},
    "2": {"analogue": "undetermined rare substate", "class": "insufficient evidence", "support": "One cell only; mTOR/MDK program signal.", "conflict": "No patient or dataset-level reproducibility; no formal DE support.", "lam": "Cannot establish LAM relationship.", "uncertainty": "Very high; do not assign a biological analogue."},
    "3": {"analogue": "myeloid/inflammatory-like, provisional", "class": "insufficient evidence", "support": "Four cells; macrophage_TREM2_TYROBP and inflammatory programs.", "conflict": "Only two patients and no formal DE/LOO biological support.", "lam": "Not sufficient for LAM interpretation.", "uncertainty": "Very high."},
    "4": {"analogue": "mixed immune–mesenchymal/interstitial", "class": "mixed or uncertain", "support": "LAM-myogenic/uterine-smooth and normal-interstitial programs; top markers include IGLV3-21, TNFRSF17, F13A1 and CD163L1.", "conflict": "Immune/plasma-cell-like markers conflict with a simple LAM interpretation; upstream cell_type unknown.", "lam": "May contain LAM-shared myogenic signal but is not LAM-specific.", "uncertainty": "High; mixed profile and no single decisive lineage assignment."},
    "5": {"analogue": "macrophage/myeloid-like", "class": "relatively clear normal-lineage analogue", "support": "FABP4, APOC1, MARCO, RETN and CCL23; macrophage_TREM2_TYROBP program delta 2.238, the strongest program signal in the atlas.", "conflict": "LAM-myogenic, mTOR and protease programs are also elevated, reflecting shared tissue programs or admixture.", "lam": "Not LAM-core; LAM-associated programs should not override the strong myeloid analogue.", "uncertainty": "Low-to-moderate for the broad analogue; exact macrophage subtype is not resolved."},
    "6": {"analogue": "rare myeloid/AT2-mixed substate", "class": "insufficient evidence", "support": "Three cells; macrophage and IL6/AT2-repair programs.", "conflict": "Only two patients, no formal biological reproducibility.", "lam": "Cannot distinguish LAM biology from a mixed rare state.", "uncertainty": "Very high."},
    "7": {"analogue": "AT2-like alveolar epithelial/repair", "class": "relatively clear normal-lineage analogue", "support": "SFTPC, SFTPA1, SFTPA2, KRT78 and SFTA2; IL6_AT2_repair is the leading program; Stage22 labels its adjacency ordinary lineage adjacency.", "conflict": "Some myogenic/macrophage program overlap is present in the candidate pool.", "lam": "Not a LAM state on current evidence.", "uncertainty": "Low-to-moderate; exact repair-state biology remains unresolved."},
    "8": {"analogue": "undetermined TGFβ/interstitial rare substate", "class": "insufficient evidence", "support": "One cell; TGFbeta_fibroblast, hormone and HOX/PBX signals.", "conflict": "No replication or DE support.", "lam": "Cannot establish LAM relationship.", "uncertainty": "Very high."},
    "9": {"analogue": "mixed repair/interstitial-like", "class": "mixed or uncertain", "support": "IL6_AT2_repair, hormone and TGFβ programs; 406 cells across four datasets.", "conflict": "Only five patients, no formal patient-level DE support and weak independent biological evidence.", "lam": "May reflect tissue-repair programs rather than LAM identity.", "uncertainty": "High."},
    "10": {"analogue": "rare myeloid/repair mixed state", "class": "insufficient evidence", "support": "Four cells; macrophage and IL6/AT2-repair programs.", "conflict": "Single dataset and no patient support.", "lam": "No reliable LAM conclusion.", "uncertainty": "Very high."},
    "11": {"analogue": "rare fibroblast/HOX-hormone-like state", "class": "insufficient evidence", "support": "Three cells; HOX/PBX, hormone and TGFβ programs.", "conflict": "Single dataset and no biological replication.", "lam": "No reliable LAM conclusion.", "uncertainty": "Very high."},
    "12": {"analogue": "endothelial/lymphatic endothelial-like", "class": "relatively clear normal-lineage analogue", "support": "MMRN1, TM4SF18, SELE, STAB2, CCL21 and FLT4; Stage22 dominant adjacent direction is endothelial.", "conflict": "Some LAM-myogenic/normal-interstitial programs overlap and high-recall candidate selection contributes to inclusion.", "lam": "Not LAM-core; current branch null supports ordinary lineage adjacency.", "uncertainty": "Moderate; endothelial and lymphatic components are not separated into new states."},
    "13": {"analogue": "rare HOX/CORE3-mixed substate", "class": "insufficient evidence", "support": "Six cells; HOX/PBX and CORE3_identity signals.", "conflict": "Only three datasets, no patient-level DE support and mixed programs.", "lam": "A possible LAM-like signal cannot be separated from sampling noise.", "uncertainty": "Very high."},
    "14": {"analogue": "rare LAM-myogenic/contractile-like substate", "class": "insufficient evidence", "support": "Three cells with high LAM-myogenic, uterine-smooth, CORE1 and CORE3 program deltas.", "conflict": "One dataset, no patient support and no formal DE replication.", "lam": "Interesting LAM-like signal, but not a reproducible state claim.", "uncertainty": "Very high."},
    "15": {"analogue": "LAM-rich contractile/mesenchymal candidate; provisional LAM-core anchor", "class": "LAM-associated candidate", "support": "200 cells; 7 patients/4 datasets; LAMCORE median 0.5125 in Stage18; CORE1/CORE3/LAM-myogenic/ECM programs; 49 author-style cells where author labels are available; removal of LAM1163 leaves 73 cells and LAMCORE median 0.4683.", "conflict": "LAM1163 contributes 127/200 with 6.8301-fold composition enrichment; author labels are available only in GSE135851; Stage11 biological reproducibility is 0.374758; Stage18 did not formally upgrade the anchor.", "lam": "Strongest LAM-rich state in this project, but still a provisional reference candidate rather than a diagnostic classifier or formal cross-patient anchor.", "uncertainty": "Moderate-to-high due to patient composition and incomplete independent author annotation."},
    "16": {"analogue": "immune/T-NK-adjacent mixed state; no confirmed transition", "class": "mixed or uncertain", "support": "396 cells; CORE1/CORE3/LAM-myogenic programs; corrected Stage22 local slope=-0.023931, direct empirical p=0.878244 and BH q=0.878244; per-patient slopes were negative in 7/7 represented patients and LOPO slopes in 10/10 omissions, but this direction was not separated from the matched local null.", "conflict": "Only 8 patients contribute State16 cells in the direct 1-hop branch criterion; immune signal is heterogeneous, upstream cell_type is unknown, and Stage11 biological reproducibility is 0.323088.", "lam": "Geometrically adjacent to State15 with a directional LAMCORE decrease, but current corrected evidence does not support a LAM-preserving transition label; ordinary adjacency and mixed biology remain plausible.", "uncertainty": "High; no temporal or lineage-transition evidence."},
    "17": {"analogue": "mesothelial/secretory epithelial-like, uncertain", "class": "mixed or uncertain", "support": "121 cells; ITLN1, CALB2, CPB1, CPA4 and ANXA8 among top markers; mTOR/ECM/inflammatory programs.", "conflict": "Only 3 supported patients, 3 datasets, and the marker profile is not a single clean mesothelial signature.", "lam": "No evidence for LAM-core; shared ECM/inflammatory signals are non-specific.", "uncertainty": "High."},
    "18": {"analogue": "pericyte/VSMC/smooth-muscle-like", "class": "relatively clear normal-lineage analogue", "support": "COX4I2, FOXC2, CASQ2, KCNA5, FHL5 and HIGD1B; high uterine-smooth/LAM-myogenic and CORE1 signal; strong structural stability 0.940033.", "conflict": "LAM and smooth-muscle programs overlap biologically, so ACTA2/myogenic signal alone cannot identify LAM; State18 did not meet the direct-connection branch rule.", "lam": "Important LAM mimic/comparator; current evidence favors ordinary VSMC/pericyte adjacency rather than a LAM branch.", "uncertainty": "Moderate."},
    "19": {"analogue": "undetermined rare interstitial/hormone-like substate", "class": "insufficient evidence", "support": "Two cells; hormone, normal-interstitial and LAF programs.", "conflict": "Two patients, two datasets and no biological replication.", "lam": "Cannot establish LAM relationship.", "uncertainty": "Very high."},
    "20": {"analogue": "fibroblast/lung interstitial-like", "class": "relatively clear normal-lineage analogue", "support": "PI16, SFRP2, SCARA5, MYOC, DPT, C7, COMP and CXCL14; normal-lung-interstitial/LAF/ECM programs; Stage22 ordinary lineage adjacency.", "conflict": "Macrophage/MDK/LAM-myogenic shared programs and candidate enrichment can blur boundaries.", "lam": "Not LAM-core; a principal normal interstitial/fibroblast comparator.", "uncertainty": "Low-to-moderate for broad analogue; fibroblast subtype not resolved."},
}


def stage_frame() -> pd.DataFrame:
    live_stage22 = stage22_snapshot()
    rows = []
    for original_spec in STAGE_SPECS:
        spec = dict(original_spec)
        if int(spec["stage"]) == 22:
            spec.update(
                {
                    "method": "k=30 不重新 Leiden；只使用 State15 局部 1–3-hop；branch eligibility 按 1-hop 患者数；real/null 统一 local scope，并按 patient×dataset、细胞数和五档距离结构匹配；直接经验左右尾、BH-q、每患者 slope 与 LOPO。",
                    "checkpoint": f"修正后仍入选 {live_stage22['selected_text']}；{live_stage22['evidence_lines']}；Stage22 checkpoint={live_stage22['manifest'].get('checkpoint', 'not_available')}。",
                    "impact": "依据校正后的 local geometry 重新评估分支；不再把 raw empirical p 当作唯一证据，State16 原 transition 标签被降级/撤回（若当前 q 不支持）。",
                    "later_revision": "Stage23 已按新 branch/boundary 输出重生成；Stage24 采用本次修正数字，未改变 State15、X_scVI 或 State1–20 标签。",
                }
            )
        rows.append({
            "stage": spec["stage"],
            "title": spec["title"],
            "research_question": spec["question"],
            "inputs": spec["inputs"],
            "method_and_parameters": spec["method"],
            "outputs": spec["outputs"],
            "checkpoint": spec["checkpoint"],
            "why_next": spec["impact"],
            "later_revision_or_scope": spec["later_revision"],
        })
    return pd.DataFrame(rows)


def state_frame() -> pd.DataFrame:
    summary = read_csv(PROJECT_ROOT / "results/stage7/state_consensus_state_summary.csv")
    atlas = read_csv(PROJECT_ROOT / "results/stage13/state_atlas.csv")
    repro = read_csv(PROJECT_ROOT / "results/stage11/state_reproducibility_summary.csv")
    identity = read_csv(PROJECT_ROOT / "results/stage16/new_candidate_by_old_state.csv")
    audit = read_csv(PROJECT_ROOT / "results/stage15/state_identity_summary.csv")
    programs = read_csv(PROJECT_ROOT / "results/stage10/state_program_scores.csv")
    markers = read_csv(PROJECT_ROOT / "results/stage10/state_markers.csv")
    records: list[dict[str, Any]] = []
    for state in range(1, 21):
        sid = str(state)
        row: dict[str, Any] = {"state_id": state, **STATE_INTERPRETATIONS[sid]}
        for source, frame, source_id in [
            ("consensus", summary, "consensus_state"),
            ("atlas", atlas, "state_id"),
            ("repro", repro, "state_id"),
            ("identity", identity, "consensus_state"),
            ("audit", audit, "consensus_state"),
        ]:
            if frame.empty or source_id not in frame.columns:
                continue
            matches = frame[frame[source_id].map(norm_state).eq(sid)]
            if matches.empty:
                continue
            selected = matches.iloc[0]
            for column in selected.index:
                if column in {source_id, "state_id", "consensus_state"}:
                    continue
                key = f"{source}_{column}"
                value = selected[column]
                row[key] = "" if pd.isna(value) else value
        program_rows = programs[programs["state_id"].map(norm_state).eq(sid)] if not programs.empty else pd.DataFrame()
        if not program_rows.empty:
            top = program_rows.sort_values("delta_state_minus_rest", ascending=False).head(5)
            row["top_program_deltas"] = "; ".join(
                f"{r.program_name}={float(r.delta_state_minus_rest):.3f}" for r in top.itertuples()
            )
        marker_rows = markers[markers["state_id"].map(norm_state).eq(sid)] if not markers.empty else pd.DataFrame()
        if not marker_rows.empty:
            top_markers = marker_rows[(marker_rows["direction"] == "up") & (marker_rows["padj"] < 0.05)]
            top_markers = top_markers.sort_values(["log2FoldChange", "padj"], ascending=[False, True]).head(8)
            row["top_DE_markers"] = ", ".join(
                f"{r.gene}({float(r.log2FoldChange):.2f})" for r in top_markers.itertuples()
            )
        records.append(row)
    return pd.DataFrame(records)


def external_paths() -> list[Path]:
    candidates: list[Path] = []
    for path in [
        PROJECT_ROOT.parent / "LAM-Cell-Research",
        PROJECT_ROOT.parent / "data-temp",
        PROJECT_ROOT / "data" / "upstream",
    ]:
        candidates.append(path)
    known_relative = [
        "data/processed/reproduction_core/GSE135851_core_reproduction.h5ad",
        "data/processed/external/GSE190260.h5ad",
        "data/processed/external/GSE217108.h5ad",
        "data/processed/external/GSE302356.h5ad",
        "results/program_discovery/candidate_pool_labels.csv",
        "results/program_discovery/external_GSE190260/candidate_pool_labels.csv",
        "results/program_discovery/external_GSE217108/candidate_pool_labels.csv",
        "results/program_discovery/external_GSE302356/candidate_pool_labels.csv",
        "config/known_lam_programs.yaml",
        "config/signatures.yaml",
        "manifests/donor_registry.yaml",
        "data/raw/reference/LAM_core_signature_genes.csv",
    ]
    records: list[Path] = []
    for root in candidates[:2]:
        records.extend(root / item for item in known_relative)
    records.extend([
        PROJECT_ROOT / "data/raw/reference/LAM_core_signature_genes.csv",
        PROJECT_ROOT.parent / "data-temp/LAM_core_signature_genes.csv",
    ])
    for inventory_path in [
        PROJECT_ROOT / "results/stage1_6/input_inventory.csv",
        PROJECT_ROOT / "results/stage1_6/upstream_annotation_manifest.csv",
    ]:
        inventory = read_csv(inventory_path)
        for column in ["h5ad_path", "candidate_path", "path"]:
            if column not in inventory.columns:
                continue
            for value in inventory[column].dropna().astype(str):
                if value.strip():
                    records.append(Path(value))
    return records


def artifact_frame() -> pd.DataFrame:
    paths: set[Path] = set()
    stage24_output = PROJECT_ROOT / "results" / "stage24_final"
    for pattern in ["[0-9][0-9]_*.py", "config/**/*", "manifests/**/*", "environment/**/*", "reports/**/*", "results/**/*", "data/processed/**/*", "data/interim/**/*", "README.md", "RESEARCH_PLAN.md", "pyrightconfig.json"]:
        paths.update(
            p for p in PROJECT_ROOT.glob(pattern)
            if p.is_file() and not p.is_relative_to(stage24_output)
        )
    rows: list[dict[str, Any]] = []
    for path in sorted(paths):
        script_match = re.match(r"(\d{2})_", path.name)
        stage_match = re.search(r"stage(\d+)(?:_6)?", str(path))
        if script_match:
            stage = str(int(script_match.group(1)))
        elif stage_match:
            stage = stage_match.group(1)
        elif path.name in {"README.md", "RESEARCH_PLAN.md", "pyrightconfig.json"}:
            stage = "project"
        else:
            stage = "shared"
        rows.append({
            "stage": stage,
            "artifact_role": "project_file",
            "relative_path": rel(path),
            "absolute_path": str(path.resolve()),
            "exists": True,
            "bytes": path.stat().st_size,
            "format": path.suffix.lower().lstrip(".") or "file",
            "rows_or_type": file_shape(path)[0],
            "columns": file_shape(path)[1],
            "sha256_if_small": sha256_small(path),
            "notes": "Existing artifact indexed without modification.",
        })
    for path in external_paths():
        rows.append({
            "stage": "upstream",
            "artifact_role": "declared_upstream_input",
            "relative_path": str(path),
            "absolute_path": str(path.resolve()),
            "exists": path.exists(),
            "bytes": path.stat().st_size if path.is_file() else "",
            "format": path.suffix.lower().lstrip(".") or "path",
            "rows_or_type": file_shape(path)[0] if path.is_file() else "missing_or_directory",
            "columns": file_shape(path)[1] if path.is_file() else "",
            "sha256_if_small": sha256_small(path),
            "notes": "External source candidate from project config / input inventory; not copied or modified by Stage 24.",
        })
    return pd.DataFrame(rows).drop_duplicates(subset=["absolute_path", "artifact_role"])


def findings_frame() -> pd.DataFrame:
    live_stage22 = stage22_snapshot()
    branch_observation = (
        f"After corrected direct-hop eligibility/local scope/null calibration, selected branches are {live_stage22['selected_text']}; "
        f"State16 label={live_stage22['state16_label']} (raw p={live_stage22['state16_raw_p']:.6f}, BH q={live_stage22['state16_q']:.6f})."
    )
    rows = [
        ("F01", "candidate gate", "5,238/5,378 candidates entered through marker-combo support; 1,443 had only 1-UMI support.", "Stage15", "root_cause_evidence.csv", "methodological finding", "The original high-recall gate is not a specific LAM classifier; this explains normal-lineage states in the pool.", "confirmed by read-only audit"),
        ("F02", "alias audit", "FIGF/VEGFD duplicate-pass cells=0 and alias-corrected loss=0.", "Stage15", "root_cause_evidence.csv", "negative result", "Duplicate counting is not the observed source of candidate inflation.", "confirmed"),
        ("F03", "upstream annotation", "External author-style fields are present as fields but not assayed; only GSE135851 has real positive labels.", "Stage19", "author_annotation_availability.csv", "data availability", "Zeros in three external datasets must not be interpreted as author-negative.", "confirmed"),
        ("F04", "identity calibration", "GSE190260 CORE3-positive references have 0 core recovery and 0.708077 core+boundary recovery; score shift -6.610707.", "Stage17", "root_cause_by_dataset.csv", "cross-dataset anomaly", "Dataset-aware calibration/dropout/competing penalty must be considered before any future gate revision.", "confirmed; no gate change here"),
        ("F05", "patient composition", "LAM1163 is 9.2971% of candidate pool but 63.5% of State15; enrichment 6.8301.", "Stage19", "state15_patient_composition.csv", "composition warning", "State15 is not an ordinary pooled state; its profile needs patient-aware interpretation.", "confirmed"),
        ("F06", "State15 sensitivity", "After removing LAM1163, 73 State15 cells remain and LAMCORE median is 0.4683.", "Stage19", "state15_without_LAM1163.csv", "supportive sensitivity", "The State15 profile is not fully explained by LAM1163, but the remaining sample is small.", "confirmed"),
        ("F07", "manifold revision", "Stage20 initially supported a pooled State15-centered manifold; Stage21 found independent evidence but not a robust unified manifold.", "Stages20-21", "stage20_manifold_report.md; stage21_manifold_validation_report.md", "superseded interpretation", "Use Stage21 as the later qualification while retaining Stage20 as historical checkpoint.", "preserved, not reconciled"),
        ("F08", "local branch", branch_observation, "Stage22", "branch_evidence_summary.csv", "corrected structural interpretation", "The corrected matched-null/FDR analysis no longer supports the prior State16 transition label; retain local adjacency as exploratory and do not infer a transition.", "current; prior label withdrawn"),
        ("F09", "state nomenclature", "Stage6 has 12 LAM-only grid reference clusters, Stage7 has 20 consensus states, and the old full cohort has 33 clusters.", "Stages6-7", "stage6_checkpoint.json; state_consensus_state_summary.csv", "terminology risk", "These are distinct clustering objects and must not be called interchangeable.", "clarified in glossary"),
        ("F10", "cell type metadata", "state_by_cell_type.csv reports unknown for all 5,378 cells.", "Stage7", "state_by_cell_type.csv", "data gap", "Human-cell analogues in Stage24 are expression/program interpretations, not verified upstream cell_type labels.", "carried as limitation"),
        ("F11", "normal reference", "Normal reference is available in Stage12/18/20 scope but is auxiliary and never defines LAM state count.", "Stages12,18,20", "normal_validation.csv; normal_remote_summary.csv", "scope boundary", "Normal-like/disease-distinct wording remains comparative, not a reclassification.", "clarified"),
        ("F12", "formal LAMCORE timing", "Stage16 recorded formal 777-gene signature unavailable; Stage18/21 used it after it appeared in data-temp.", "Stages16-18", "identity_gate_report.md; state15_anchor_summary.json", "historical input change", "Do not read Stage16 unavailable as proof that the signature never existed; it was unavailable at that run.", "preserved as chronology"),
        ("F13", "DE evidence", "Pathway enrichment and regulon outputs are explicit not_available placeholders.", "Stage10", "state_pathway_enrichment.csv; state_regulon_summary.csv", "method gap", "No pathway/regulon conclusion is included in final state analogue calls.", "carried as limitation"),
        ("F14", "State16 interpretation", "State16 has 80 raw-count LAM/immune coexpressing cells in Stage20, but only 2 LAM-high/immune-high cells.", "Stage20", "state16_lam_immune_coexpression.csv", "technical/biological ambiguity", "Do not call State16 a doublet state; use mixed/transitional wording and retain technical audit caveat.", "carried"),
        ("F15", "small states", "States 2,3,6,8,10,11,13,14,19 have no supported patients in Stage11 summary.", "Stage11", "state_reproducibility_summary.csv", "negative result", "Their labels are retained for completeness but not used for strong human analogues.", "confirmed"),
        ("F16", "unresolved boundary", f"{live_stage22['boundary_unresolved']} of {live_stage22['boundary_rows']} local 1–3-hop boundary cells remain unresolved in Stage22 local branch assignment; farther boundary cells are not projected.", "Stage22", "boundary_local_branch_assignment.csv", "negative/uncertain result", "Boundary is an evidence-ranking cohort, not a forced new LAM class; report local scope explicitly.", "corrected scope"),
    ]
    columns = ["finding_id", "topic", "observation", "stage", "source_file", "finding_type", "implication", "disposition"]
    return pd.DataFrame(rows, columns=columns)


def audit_frame() -> pd.DataFrame:
    live_stage22 = stage22_snapshot()
    rows = [
        ("A01", "Stage20 checkpoint says supports_lam_centered_transcriptional_manifold; Stage21 says state15_lam_rich_gradient_but_not_robust_manifold.", "high", "Stage20 report; Stage21 report", "Retain both chronologically; final synthesis uses Stage21 qualification and Stage22 local branch result.", "historical evolution"),
        ("A02", "Stage6 12 clusters, Stage7 20 consensus states, and full-cohort 33 clusters coexist.", "high", "stage6_checkpoint.json; state_consensus_state_summary.csv", "Add explicit object definitions; do not state that the project has only 12 or only 20 total clusters.", "terminology clarification"),
        ("A03", "Stage16 says formal 777 signature unavailable; Stage18 says available.", "high", "stage16/identity_gate_report.md; stage18/state15_anchor_summary.json", "Explain availability changed when data-temp file arrived; Stage16 was not rerun.", "chronology"),
        ("A04", "External author-style fields are false/present but not_assayed.", "high", "stage19/author_annotation_availability.csv", "Use not_assayed, never author-negative; enrichment is only formally assessed in GSE135851.", "definition correction"),
        ("A05", "Stage21 independent score has negative patient-adjusted slope but positive pooled Spearman and near/far medians not monotonically decreasing.", "high", "stage21/gradient_models.csv; non_state15_distance_gradient.csv", "Report as scope/shape dependence and patient heterogeneity; do not algebraically reconcile as if they were the same estimand.", "statistical qualification"),
        ("A06", f"Stage22 corrected State16 raw empirical p={live_stage22['state16_raw_p']:.6f}, BH q={live_stage22['state16_q']:.6f}; the prior p=0.025948 came from the superseded scope/null method.", "high", "stage22/branch_evidence_summary.csv; branch_patient_lopo.csv", "Withdraw the prior State16 transition label; current evidence is ordinary-lineage adjacency under corrected local matched-null/FDR analysis.", "method correction and conclusion withdrawal"),
        ("A07", "All state_by_cell_type values are unknown, while reports sometimes use human cell analogues.", "medium", "stage7/state_by_cell_type.csv; Stage13/Stage24 interpretation", "Use analogue/inferred lineage terminology, not verified cell_type.", "terminology clarification"),
        ("A08", "Stage13 novel_or_unexplained is false for all rows, but this does not mean all states are biologically explained.", "medium", "stage13/state_atlas.csv", "Treat this field as the atlas generation flag, not as proof of identity or absence of uncertainty.", "interpretation boundary"),
        ("A09", "Stage12 normal and Stage20 normal_remote use different named scopes.", "low", "stage12/state_auxiliary_summary.csv; stage20/normal_remote_summary.csv", "Describe normal as auxiliary comparison and preserve the stage-specific scope.", "scope clarification"),
        ("A10", "Stage10 pathway/regulon files exist but contain not_available placeholders.", "medium", "stage10/state_pathway_enrichment.csv; state_regulon_summary.csv", "Do not imply these analyses were completed.", "missing evidence"),
        ("A11", "Stage15 state summaries show FIGF in marker combinations, but duplicate audit is zero.", "medium", "stage15/state_identity_summary.csv; alias_audit_by_state.csv", "Explain that FIGF can appear as the original marker label while no same-cell FIGF+VEGFD double counting was observed.", "alias wording"),
        ("A12", "Stage18/19 anchor language remains provisional, not a formal anchor artifact.", "medium", "stage18/state15_anchor_report.md; state15_anchor_summary.json", "Final report uses LAM-core reference-anchor candidate.", "strength downgrade"),
    ]
    return pd.DataFrame(rows, columns=["audit_id", "issue", "severity", "evidence", "stage24_disposition", "audit_type"])


def glossary() -> str:
    return """# 术语和结论演变

| 术语 | 本项目最终使用方式 | 演变/边界 |
|---|---|---|
| `pool_high_confidence` | Stage 6–23 主线 `lam_candidate`，固定 5,378 cells | 不能与 broad/unrestricted 做并集；数量不足时不回退 |
| `pool_broad_lam_like` | upstream annotation；与 high 不重叠部分构成 boundary | 只作连续性/边界辅助，不自动成为 LAM |
| `pool_unrestricted_lam` | 审计标签 | 等同所有 condition=LAM 细胞，不能进入默认模型 cohort |
| candidate | 高召回、可供 State Modeling 检查的输入集合 | 不等于已证实 LAM；Stage15 证明 marker-combo gate 特异性不足 |
| consensus state | Stage7 的 20 个 frozen state labels | 不等同 Stage6 的 12 个 LAM-only grid clusters，也不等同 full-cohort 33 clusters |
| State15 | 200-cell frozen LAM-rich candidate state | 从“待验证 anchor”变成“provisional reference-anchor candidate”，未升级为正式 classifier |
| State16 | State15 邻近的 396-cell frozen state | 修正后仅保留为局部几何邻接/混合状态候选；不再标为 `LAM_to_lineage_transition_candidate` |
| manifold | Stage20 提出的全局 State15-centered gradient hypothesis | Stage21 削弱为非稳健统一 manifold；Stage22 保留局部分支 candidate |
| ordinary lineage adjacency | 与 State15 几何相邻但 matched-null/gradient 不支持 LAM branch | 当前用于 State12、20、7 |
| author-style | 上游真实作者逐细胞标签 | 只有 GSE135851 `available`；其余三个 dataset 是 `not_assayed` |
| formal LAMCORE | 777-gene formal signature score | Stage16 运行时 unavailable，Stage18 后在 data-temp 可用；不能倒写 Stage16 |
| structural stability | LOO/partition recovery 等结构维度 | 与 biological reproducibility 分开报告 |
| biological reproducibility | patient-aware DE/program/profile 的复现维度 | 不由 cluster 数或单个 p 值替代 |
| human cell analogue | 基于 marker/program/latent/normal 对照的解释 | 因 upstream `cell_type` 全 unknown，不是已验证人体细胞注释 |
| evidence vs hypothesis | 直接结果、解释、后续 candidate 分开 | Stage24 不把 provisional/uncertain 升级为事实 |

## 结论状态变化

`Stage 6 GO` → `Stage 7 20-state consensus` → `Stage 15 gate-specificity problem` → `Stage 18 State15 provisional anchor` → `Stage 20 global manifold hypothesis` → `Stage 21 non-robust global manifold` → `Stage 22 corrected local branch audit; State16 transition label withdrawn`.

这条链不是矛盾需要消除，而是研究问题在新证据下被逐步收窄。
"""


def state_markdown(states: pd.DataFrame) -> str:
    lines = [
        "# 20 个 consensus state 的生物学解释",
        "",
        "以下 analogue 是基于现有 DE、program、LAMCORE、LOO、normal/boundary 和 latent 邻域的综合解释，不是 upstream 已验证 `cell_type`。由于 `state_by_cell_type.csv` 的 5,378 个细胞均为 `unknown`，所有命名都保留 evidence boundary。",
        "",
    ]
    for _, row in states.iterrows():
        sid = int(row["state_id"])
        lines.extend([
            f"## State {sid} — {row['analogue']}",
            "",
            f"- 解释类别：{row['class']}。",
            f"- 支持证据：{row['support']}",
            f"- 冲突证据：{row['conflict']}",
            f"- 与 LAM 的关系：{row['lam']}",
            f"- 不确定性：{row['uncertainty']}",
            f"- 细胞/覆盖/结构证据：{row['consensus_cells']} cells；{row['consensus_patients']} patients；{row['consensus_datasets']} datasets；structural={row['repro_structural_stability']}；biological={row['repro_biological_reproducibility']}。",
            f"- 表达/程序摘要：top DE markers={row['top_DE_markers'] if pd.notna(row['top_DE_markers']) else 'not available'}；top program deltas={row['top_program_deltas'] if pd.notna(row['top_program_deltas']) else 'not available'}。",
            "",
        ])
    return "\n".join(lines)


def state_table_markdown(states: pd.DataFrame) -> str:
    columns = [
        "state_id", "consensus_cells", "consensus_patients", "consensus_datasets",
        "analogue", "class", "identity_LAM_core_fraction", "identity_non_LAM_like_fraction",
        "repro_structural_stability", "repro_biological_reproducibility", "atlas_boundary_connectivity",
        "atlas_normal_distance", "atlas_patient_coverage", "top_DE_markers", "top_program_deltas",
    ]
    available = [c for c in columns if c in states.columns]
    table = states[available].copy()
    return markdown_table(table)


def source_materials(states: pd.DataFrame, stages: pd.DataFrame, findings: pd.DataFrame, audits: pd.DataFrame) -> str:
    stage_table = markdown_table(stages[["stage", "title", "research_question", "checkpoint", "why_next", "later_revision_or_scope"]])
    finding_table = markdown_table(findings)
    audit_table = markdown_table(audits)
    state_table = state_table_markdown(states)
    live_stage22 = stage22_snapshot()
    stage22_evidence = live_stage22["evidence_lines"]
    stage22_short = (
        f"修正后的 Stage22 仍选择 {live_stage22['selected_text']}；其 State16 为 {live_stage22['state16_label']} "
        f"（raw empirical p={live_stage22['state16_raw_p']:.6f}，BH q={live_stage22['state16_q']:.6f}）。"
    )
    return f"""# LAM-State-Modeling：Stage 24 最终报告原材料包

> 这是面向最终写作的证据材料汇总，不是对历史 artifact 的覆盖。Stage 24 只读取已有结果并生成本目录文件；不重训 scVI、不重新聚类、不修改 candidate gate、State15、State1–20 或 branch artifact。历史矛盾和结论变化按时间保留。

## 1. 研究目标与目标演变

最初问题是：能否在跨数据集、跨患者的 LAM high-confidence candidate 中，用深度 latent representation 找到可重复的 LAM state structure。阶段 1–6 先完成输入继承、外部 QC、NMF/PCA baseline、scVI 和 LAM-only GO/NO-GO。Stage6 的 GO 只说明 high-confidence candidate 内存在值得继续检查的 latent structure，并不证明这些 cell 都是 LAM。

Stage7–13 将目标变成：哪些结构跨参数、患者和数据集稳定，哪些 state 能得到 patient-aware biological support。随后 Stage15 暴露出原 candidate gate 的高召回/低特异性，Stage16–17 改为独立 identity audit，而不是把新 gate 偷换进主模型。Stage18–19 冻结 State15，验证它是否能作为 LAM-core reference anchor；结果支持“LAM-rich、患者丰度异质、去掉 LAM1163 后 profile 仍保留”的 provisional interpretation，但没有正式升级。

Stage20 进一步提出 State15-centered global manifold。Stage21 在排除 anchor、scVI HVG/gate overlap 和 composition matched null 后保留了一部分独立 gradient，但否定了“稳健统一 global manifold”的强表述。{stage22_short} State12、20、7 的标签和 State15 局部关系也以修正后的结果为准。Stage23 仅把这些结果可视化。因而 Stage24 的最终科学问题不是“20 个 cluster 是否都是 LAM”，而是：一个高置信度候选池中，State15 是否代表最强 LAM-rich core，以及其周围是否存在有限、局部且方向依赖的延伸。

## 2. 数据和分析范围

主要数据集为 GSE135851（core reproduction）、GSE190260、GSE217108、GSE302356（external converted AnnData）。Stage1 inventory 记录原始转换对象大小：分别为 30,708、39,979、12,396、23,759 cells；pool_high_confidence 分别为 535、1,564、1,075、2,681，合计 5,855 个 upstream high flags。经过继承/QC/共同 gene universe 后，主线 Stage6–23 使用 5,378 个 high-confidence candidate cells、12 patients、4 datasets。

Stage1–6 inventory 中的 pool counts 为：

| dataset | AnnData cells | genes | high | broad | unrestricted |
|---|---:|---:|---:|---:|---:|
| GSE135851 | 30,708 | 63,677 | 535 | 3,940 | 23,228 |
| GSE190260 | 39,979 | 33,694 | 1,564 | 7,811 | 39,979 |
| GSE217108 | 12,396 | 36,601 | 1,075 | 3,855 | 12,396 |
| GSE302356 | 23,759 | 38,224 | 2,681 | 8,727 | 23,759 |

外部补充 QC 后 pass cells 为 GSE190260 38,701、GSE217108 6,941、GSE302356 21,771；GSE135851 直接继承原 reproduction baseline `qc_pass`。Stage12 normal auxiliary cohort 为 32,977 cells；boundary 为 16,883 cells。State15 为冻结 200 cells、7 patients、4 datasets；State16 为 396 cells。normal/boundary 从未参与 State6–7 的核心 state 数量定义。

候选池语义固定为：`lam_candidate = pool_high_confidence`；`boundary = pool_broad_lam_like AND NOT pool_high_confidence`；`pool_unrestricted_lam` 只审计。由于 unrestricted 近似所有 condition=LAM cell，三层不能做并集。

## 3. 核心方法和关键参数

- 数据继承：优先读取 `../LAM-Cell-Research`，否则 data-temp/data/upstream；不重新 GEO 转换；patient/donor mapping 唯一来源为 `donor_registry.yaml`。
- Alias：导入后、gene 去重前执行 `FIGF → VEGFD`；Stage15 审计显示同一细胞 FIGF+VEGFD duplicate-pass=0。
- QC：core 继承原 `qc_pass`；外部 min_genes=200、min_counts=500、scRNA mt<20%、snRNA mt<10%；doublet 保留、不重跑 caller、不默认删除。
- NMF：counts→library-size normalization→log1p→HVG=4000→top 2000 features→5 components，max_iter=400，最多 12,000 cells。
- scVI：`layer='counts'`、`batch_key='dataset'`、categorical covariates=[]、20 latent、2 layers、128 hidden、200 epochs、early stopping；assay 仅 metadata。后续 Stage7–23 全部复用同一个 `X_scVI`。
- Stage6：LAM-only candidate graph 的参数独立于 preprocess；9 个 configuration 为 n_neighbors 15/30/50 × resolution 0.2/0.4/0.6。
- Stage7：21 raw partitions 先在 configuration 内平均 seed，再 9 configuration 等权；最终 co-assignment 允许一份 5,378×5,378 float32 matrix 和完整 average-linkage。
- Stage8：full-data reference 为 n_neighbors=30、resolution=0.4；LOO overlap 只在 retained cells 上计算，并区分 full-reference→consensus baseline 与 LOO additional loss。
- Stage10：每个 state 独立 `State_k vs Rest_of_LAM`；patient×group pseudobulk；设计 `~ patient_id + group`；正式 DE 至少 3 patients；不使用统一多分类模型。
- Stage18–22：State15 的 200 cells 固定为 anchor；Stage20 distance 完全由 X_scVI；Stage21 的 777 LAMCORE validation 拆为 full/no_gate/outside_scVI/independent；Stage21 matched null 500 次；Stage22 branch null 每条 500 次、局部 k=30、不重新 Leiden。Stage22 修正版用 1-hop 患者数筛选、real/null 相同的 1–3-hop local scope、patient×dataset+距离分箱匹配、直接经验尾部、BH-q 和患者级 LOPO。

## 4. Stage 1–23 完整过程

{stage_table}

每个阶段的完整 artifact 路径、行数/文件大小和是否存在见 `artifact_index.csv`；每一阶段的简短可引用版本见 `stage_summary.md`。

## 5. 最终 state evidence

以下表同时保留数值证据和解释性字段。`analogue` 是 human-cell analogue，不是 upstream verified cell_type。State15、16、12、20、7、5 的详细解读见 `state_human_cell_analogue.md`。

{state_table}

### State15 的核心证据链

1. Unsupervised：Stage7 consensus 中出现固定 200-cell State15。
2. Identity：Stage18 formal 777-gene LAMCORE median=0.5125；normal median=0.0513；Stage13 program correspondence 包含 CORE1、CORE3_identity、LAM_myogenic_contractile 和 ECM_remodeling。
3. Author evidence：总体 enrichment fold=13.9408、one-sided Fisher p=2.229e-36；但真实逐细胞 author annotation 只有 GSE135851 可用，49/50 State15 cells 为 author-style positive，外部三个 dataset 为 not_assayed。
4. Patient composition：LAM1163 在 candidate pool 中 9.2971%，在 State15 中 63.5%，enrichment=6.8301；这排除了“只是因为样本多”的解释。
5. Sensitivity：去除 LAM1163 后剩 73 cells，LAMCORE median=0.4683，仍高于指定 comparators；但样本更小，不能替代独立 replication。
6. Structure：Stage11 structural=0.854111、biological=0.374758、patient direction concordance=0.709647、patient coverage=0.416667；这些支持存在 LAM-rich candidate state，同时说明跨患者生物学证据仍有限。
7. Geometry：Stage20 global gradient 后经 Stage21 独立检验被削弱；Stage22 修正后的 local matched-null/FDR 分析不再支持 State16 transition label。因此 State15 最稳妥的称呼是 `provisional LAM-core reference-anchor candidate`。

### State16 及关键 comparator

- State16：396 cells，4 datasets，Stage22 1-hop=96、直接 1-hop 患者数=8（全 state 覆盖 10 patients）；修正后 local independent LAMCORE slope=-0.023931、直接经验 p=0.878244、BH q=0.878244，当前不再是 `LAM_to_lineage_transition_candidate`，而是 State15 邻接/混合状态的探索性描述。
- State18：pericyte/VSMC/smooth-muscle-like，COX4I2/FOXC2/CASQ2 等支持；结构稳定但未进入 Stage22 direct branch selection，保留为重要 LAM mimic/comparator。
- State20：PI16/SFRP2/SCARA5/DPT/C7/COMP/CXCL14，fibroblast/lung interstitial-like；Stage22 ordinary lineage adjacency。
- State12：MMRN1/CCL21/FLT4 等 endothelial/lymphatic-like；ordinary lineage adjacency。
- State7：SFTPC/SFTPA1/SFTPA2 等 AT2-like；ordinary lineage adjacency。
- State5：FABP4/APOC1/MARCO/RETN 等 macrophage-like，尽管有 LAM-shared mTOR/protease/myogenic programs。

## 6. Candidate gate、State15 和 manifold 的完整证据链

Stage15 的 A/B/C 审计是主线必须保留的校准点：annotation/ID merge field inconsistency=0；FIGF/VEGFD alias duplicate-pass=0；5,238 cells 通过 marker-combo、140 cells 有 author/formal support、formal support alone=0；1,443 marker-combo cells 仅由 1-UMI detections 支持。因此“明显非 LAM-like state 被纳入”主要由 C 类 gate 特异性不足解释，而不是 State Modeling 合并搬错。

Stage16 的独立 continuous gate 不使用现有 20 states 调参；其输出为诊断 artifact，不能回写主线。Stage17 进一步显示 GSE190260 positive reference 的 identity/support dropout 与 competing-lineage penalty/score shift 同时存在：2,117 positive references 中 core recovery=0，core+boundary=0.708077，median final score=-1.006577，relative score shift=-6.610707，median competing penalty=4.950320。因此跨数据集零 author-style 不能当作不支持。

Stage20 初始 pooled result 观察到 full LAMCORE 最近/最远 median 0.2152/0.1556、4/4 dataset rho<0、8/12 patient rho<0，提出 global manifold。Stage21 把 State15 排除为 anchor，使用 22,061 non-State15 cells，并审计 777 formal genes：220 与 scVI HVG 重叠，7 属于旧 gate markers，554 同时不属于二者，729 在表达矩阵可用。candidate-only independent slope=-0.015466（95% CI -0.017451,-0.013482，p=1.410338e-51），500 matched fake-anchor empirical two-sided p=0.001996；但 pooled Spearman rho=+0.075665，非-State15 full scope independent rho=+0.094447，距离箱的独立 score 不显示稳健单调下降，且 patient direction heterogeneous。Stage21 因此把结论降为 `state15_lam_rich_gradient_but_not_robust_manifold`。

Stage22 修正后仍只选出 State16、12、20、7 四个方向，但 real/null 都限制在局部 1–3 hop，并进一步匹配距离结构；各分支结果为：{stage22_evidence}。State16 的旧 transition 标签在直接经验尾部和 BH 校正后不再成立；10,645 的旧 boundary 总体数字也不能直接沿用，当前只报告 {live_stage22['boundary_unresolved']} / {live_stage22['boundary_rows']} 个 local 1–3-hop boundary unresolved。最终不再支持 State15→State16 的 transition candidate；Stage22 仅保留为校正后的 local adjacency diagnostic，而非统一 global manifold 或 temporal trajectory。

## 7. 阴性结果、被否定假设和异常

{finding_table}

不能保留的强表述包括：“所有 20 states 都是 LAM”“State15→State16 是时间转化”“所有 State15 邻近 branch 都是 LAM”“candidate gate 可作诊断 classifier”“Stage20 global manifold 已证实”。这些表述分别被 candidate audit、Stage21/22、author availability 和患者复现结果削弱。

## 8. 文字和逻辑审计

{audit_table}

Stage24 采用的原则是：后来的结果可以限定早期结论，但不删除早期 checkpoint；不同 estimand 的数值不强行调和；缺失的 pathway/regulon 和 upstream cell_type 不补造；formal 777 signature 的可用性按运行时间记录。

## 9. 结论边界与局限

- 数据直接支持：Stage15 是当前 candidate pool 中最 LAM-rich 的 frozen state；它富集 formal LAMCORE/author evidence（仅在可 assay dataset）、并在去除 LAM1163 后保留部分 profile。
- 支持性解释：State15 可能是 LAM-core reference-anchor candidate；State16 与 State15 存在方向性几何邻接，但当前不支持 LAM-preserving/immune-direction transition candidate。
- 尚未证明：State15 是跨患者正式 reference anchor；State15→State16 是时间或谱系转化；存在稳健统一 global manifold；candidate gate 可作为临床/诊断 classifier。
- 患者数量和 composition：State15 只有 7/12 patients，且 LAM1163 enrichment=6.8301；几个小 state 没有 supported patient。
- Dataset heterogeneity：GSE190260 的 identity score shift/dropout 明显；author-style annotation 只在一个 dataset available；scRNA/snRNA、测序深度和转换来源不同。
- Candidate enrichment：主线 high-confidence gate 经 Stage15 证明包含大量普通肺细胞/lineage-like clusters；这既是限制，也让后续 latent space 暴露了 gate 的边界。
- 技术与统计：DE 依赖有限 patient pseudobulk；pathway/regulon 未完成；局部分支 p 值为探索性；distance gradient 不同 scope/estimand 方向不完全一致。
- 生物学验证缺口：没有真实时间维度、空间关系、实验验证、独立 prospective cohort，也没有把 ATAC/spatial 加入本轮验证。

## 10. 未来方向

1. 在独立、跨数据集且 author/formal reference 同步 available 的 cohort 验证 State15 anchor。
2. 对 State16 进行独立实验/空间或正交 modality 验证，区分 LAM-like extension、immune adjacency 与 mixed/doublet profile。
3. 以 patient-aware、dataset-calibrated 的 continuous identity model 重新评估 candidate，而不是恢复任意两个 marker 阳性门槛。
4. 用空间/ATAC/时间序列检验 local branch 是否有真实组织位置、染色质或动态方向。
5. 对普通 lineage comparator（AT2、endothelial、fibroblast、VSMC、macrophage）建立独立参考，减少 candidate enrichment 对 LAM identity 的混淆。
6. 对现有 20 states 的 pathway/regulon 和更大患者数复现补齐，但不在 Stage24 重新训练/聚类。

## 11. 追溯入口

- `stage_index.csv`：Stage1–23 的问题、输入、方法、checkpoint、下一步影响和后续修正。
- `artifact_index.csv`：项目脚本、报告、表格、模型/AnnData 和声明的 upstream paths；包含存在性、大小、行数/类型和小文件 SHA-256。
- `state_human_cell_analogue.csv/md`：20 states 的解释及支持/冲突证据。
- `other_findings_registry.csv/md`：主线外发现、异常、负结果和方法学问题。
- `narrative_audit.csv`：定义变化、数字冲突、历史结论升级/降级和处理方式。
- `appendices/`：方法、state 详细表、统计摘录和 provenance。
"""


def stage_summary_markdown(stages: pd.DataFrame) -> str:
    sections = ["# Stage 1–23 summary", "", "本文件是因果链摘要；完整证据路径见 `artifact_index.csv`，历史冲突见 `narrative_audit.csv`。", ""]
    for row in stages.itertuples(index=False):
        sections.extend([
            f"## Stage {int(row.stage)} — {row.title}",
            "",
            f"**问题**：{row.research_question}",
            f"**输入**：{row.inputs}",
            f"**方法/参数**：{row.method_and_parameters}",
            f"**输出**：{row.outputs}",
            f"**结果/checkpoint**：{row.checkpoint}",
            f"**为什么进入下一步**：{row.why_next}",
            f"**后续修正/边界**：{row.later_revision_or_scope}",
            "",
        ])
    return "\n".join(sections)


def report_text(stages: pd.DataFrame, states: pd.DataFrame, findings: pd.DataFrame, audits: pd.DataFrame) -> str:
    return source_materials(states, stages, findings, audits)


def write_appendices(output: Path, stages: pd.DataFrame, states: pd.DataFrame, artifacts: pd.DataFrame, findings: pd.DataFrame, audits: pd.DataFrame) -> None:
    live_stage22 = stage22_snapshot()
    methods = """# Detailed methods appendix

## Input contract

The project inherited converted AnnData and upstream annotations from `LAM-Cell-Research` where available, with `data-temp` as the supplied fallback. GSE135851 inherited its upstream QC status; the three external AnnData objects received project-specific supplementary QC. The canonical candidate rule remained `pool_high_confidence`; broad-minus-high was boundary; unrestricted was audit-only.

## Matrix separation

NMF used raw counts followed by library-size normalization, log1p, State Modeling HVG selection and NMF. scVI used `layers["counts"]` as raw counts with `batch_key="dataset"` only. `assay` remained metadata. Later stages used the existing `X_scVI` and did not retrain.

## Statistical scope

Stage10 fitted each state independently against same-patient Rest_of_LAM using patient×state pseudobulk and `~ patient_id + group`. Stage8 LOO metrics were computed on retained cells only. Stage21/22 nulls were composition/cell-count matched exploratory controls, not replacements for independent biological replication.

## Stage24 boundary

Stage24 is a read-only synthesis stage. It performs deterministic reading, table joining and document audit. It does not run scanpy clustering, scVI training, DE, or new biological discovery.
"""
    stats = f"""# Statistical evidence appendix

| Evidence item | Value | Source |
|---|---:|---|
| Stage6 high-confidence candidate cells | 5,378 | `results/stage1_6/stage6_checkpoint.json` and grid tables |
| Stage6 LAM-only reference clusters | 12 | `results/stage1_6/stage6_checkpoint.json` |
| Stage6 patient/dataset/assay ARI | 0.152554 / 0.086786 / 0.009094 | `results/stage1_6/stage6_driver_ari.csv` |
| Stage7 consensus states | 20 | `results/stage7/state_consensus_state_summary.csv` |
| Stage10 DE rows | 159,059 | `results/stage10/state_de_results.csv` |
| Stage15 marker-combo candidates | 5,238 | `results/stage15/root_cause_evidence.csv` |
| Stage15 1-UMI-only marker-combo cells | 1,443 | `results/stage15/root_cause_evidence.csv` |
| Stage18 State15 LAMCORE median | 0.5125 | `results/stage18/state15_anchor_report.md` |
| Stage19 LAM1163 composition enrichment | 6.8301 | `results/stage19/state15_patient_composition.csv` |
| Stage21 candidate-only independent slope | -0.015466 | `results/stage21/gradient_models.csv` |
| Stage21 matched-null empirical two-sided p | 0.001996 | `results/stage21/matched_anchor_null.csv` |
| Stage22 State16 independent branch slope | {live_stage22['state16_slope']:.6f} | `results/stage22/branch_evidence_summary.csv` |
| Stage22 State16 corrected matched-null raw p | {live_stage22['state16_raw_p']:.6f} | `results/stage22/branch_evidence_summary.csv` |
| Stage22 State16 corrected BH q | {live_stage22['state16_q']:.6f} | `results/stage22/branch_evidence_summary.csv` |

The Stage21 slope, Spearman rho and binned medians are different estimands/scopes. Stage22 uses the corrected local 1–3-hop scope, distance-structure matched null, direct empirical tails and BH correction; its branch labels are not comparable to the superseded raw-p-value-only version. Their signs and shapes are intentionally reported without forced reconciliation.
"""
    atlas = """# State atlas detailed appendix

This appendix preserves every consensus state, including states with little or no patient-level support. Analogue labels are interpretive and not verified `cell_type` labels.

""" + state_table_markdown(states)
    provenance = """# Artifact provenance appendix

The authoritative machine-readable provenance is `artifact_index.csv`. It lists all project scripts, reports, result tables, AnnData/model files and declared upstream paths, with existence, file size, CSV row/column counts where applicable, and SHA-256 for small files. Large binary artifacts are intentionally not hashed during Stage24 to avoid unnecessary I/O and memory pressure.

Core frozen artifacts include:

- `data/processed/state_model_scvi.h5ad` and `data/processed/scvi_model/`
- `results/stage7/state_consensus_assignments.csv`
- `results/stage7/state_consensus_state_summary.csv`
- `results/stage13/state_atlas.csv`
- `results/stage18/state15_anchor_summary.json`
- `results/stage20/state15_centered_manifold.csv`
- `results/stage21/gradient_models.csv`
- `results/stage22/branch_evidence_summary.csv`
- `results/stage23_visualization/visualization_manifest.json`

Stage24-generated files are confined to `results/stage24_final/`.
"""
    write_text(output / "appendices/methods_detailed.md", methods)
    write_text(output / "appendices/statistics.md", stats)
    write_text(output / "appendices/state_atlas_detailed.md", atlas)
    write_text(output / "appendices/artifact_provenance.md", provenance)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output_dir if args.output_dir.is_absolute() else PROJECT_ROOT / args.output_dir
    output.mkdir(parents=True, exist_ok=True)

    stages = stage_frame()
    states = state_frame()
    findings = findings_frame()
    audits = audit_frame()
    artifacts = artifact_frame()

    write_frame(output / "stage_index.csv", stages)
    write_text(output / "stage_summary.md", stage_summary_markdown(stages))
    write_frame(output / "artifact_index.csv", artifacts)
    write_frame(output / "state_human_cell_analogue.csv", states)
    write_text(output / "state_human_cell_analogue.md", state_markdown(states))
    write_frame(output / "other_findings_registry.csv", findings)
    write_text(output / "other_findings.md", "# Other findings registry\n\n" + markdown_table(findings))
    write_frame(output / "narrative_audit.csv", audits)
    write_text(output / "terminology_glossary.md", glossary())

    full_report = report_text(stages, states, findings, audits)
    write_text(output / "final_project_source_materials.md", full_report)
    write_text(output / "final_project_report.md", full_report)
    live_stage22 = stage22_snapshot()
    write_text(output / "final_project_summary.md", f"""# 最终项目摘要

本项目最终支持的最稳妥结论是：在一个高召回、包含明显普通肺谱系的 LAM candidate pool 中，State15 是当前最 LAM-rich 的 frozen consensus state，并具有 formal LAMCORE、可用数据集中的 author-label enrichment、患者匹配及去除 LAM1163 后仍保留的 profile。它尚未达到独立、均衡跨患者的正式 reference anchor 标准。

Stage20 的 pooled State15-centered gradient 在 Stage21 被限定为“存在部分独立 LAM-rich gradient，但不是稳健统一 global manifold”。Stage22 修正 branch eligibility、local scope、距离匹配、经验尾部和多重比较后，当前选中 {live_stage22['selected_text']}；State16 的当前标签为 {live_stage22['state16_label']}（raw p={live_stage22['state16_raw_p']:.6f}，BH q={live_stage22['state16_q']:.6f}）。因此不再把 State15→State16 称为 transition candidate。没有证明时间转化、诊断 classifier 或所有 candidate state 都是 LAM。

完整材料见 `final_project_source_materials.md`；数字和文件来源见 `artifact_index.csv`、`stage_index.csv` 与 `narrative_audit.csv`。
""")
    write_text(output / "limitations.md", f"""# Limitations

- State15 仅覆盖 7/12 patients，且 LAM1163 composition enrichment=6.8301。
- External author-style annotation is `not_assayed` for GSE190260/GSE217108/GSE302356.
- Upstream `cell_type` is unknown for all consensus cells; Stage24 analogues are inferred.
- Candidate gate is high-recall and includes ordinary lineage states.
- Dataset heterogeneity, marker dropout and GSE190260 score shift limit cross-dataset calibration.
- No time, spatial, prospective cohort or experimental validation was included.
- Pathway enrichment and regulon outputs are unavailable placeholders.
- Stage21 pooled and patient-adjusted gradient estimands differ; corrected Stage22 branch p-values remain exploratory and State16 no longer passes the corrected matched-null/FDR evidence label.
- Stage22 branch selection and null analysis are local 1–3-hop analyses; {live_stage22['boundary_rows']} local boundary cells were projected, while farther boundary cells were intentionally not assigned.
""")
    write_text(output / "future_directions.md", """# Future directions

1. Validate State15 in an independent cohort with synchronized formal and author-level references.
2. Test State16 with spatial, chromatin or experimental evidence as a local adjacency/mixed-state hypothesis; do not assume a transition from the corrected Stage22 result.
3. Develop dataset-calibrated continuous identity evidence rather than restoring an arbitrary two-marker gate.
4. Use independent normal-lineage references to separate LAM identity from VSMC, fibroblast, endothelial, AT2 and myeloid programs.
5. Increase patient-level replication and complete pathway/regulon analyses.
6. Test whether the local State15→State16 direction has spatial or temporal meaning; do not assume a trajectory from the current latent geometry.
""")
    write_appendices(output, stages, states, artifacts, findings, audits)

    expected_core = [
        "data/processed/state_model_scvi.h5ad",
        "results/stage7/state_consensus_assignments.csv",
        "results/stage13/state_atlas.csv",
        "results/stage18/state15_anchor_summary.json",
        "results/stage20/state15_centered_manifold.csv",
        "results/stage21/gradient_models.csv",
        "results/stage22/branch_evidence_summary.csv",
        "results/stage23_visualization/visualization_manifest.json",
    ]
    missing_core = [item for item in expected_core if not (PROJECT_ROOT / item).exists()]
    manifest = {
        "stage": 24,
        "generated_date": date.today().isoformat(),
        "mode": "read_existing_artifacts_and_write_stage24_package_only",
        "project_root": str(PROJECT_ROOT),
        "output_dir": str(output.resolve()),
        "stages_indexed": list(range(1, 24)),
        "stage_count": 23,
        "state_count": int(len(states)),
        "core_artifacts_frozen": True,
        "no_scvi_training": True,
        "no_reclustering": True,
        "no_candidate_gate_change": True,
        "no_state15_change": True,
        "no_branch_change": True,
        "stage22_corrected_downstream_regenerated": ["results/stage23_visualization", "results/stage24_final"],
        "stage22_prior_state16_transition_label_withdrawn": True,
        "stage22_current_checkpoint": live_stage22["manifest"].get("checkpoint", "not_available"),
        "history_conflicts_preserved": True,
        "expected_core_artifacts_missing": missing_core,
        "artifact_index_rows": int(len(artifacts)),
        "generated_files": sorted(
            [rel(path) for path in output.rglob("*") if path.is_file()]
            + [rel(output / "stage24_manifest.json")]
        ),
        "primary_frozen_inputs": expected_core,
    }
    write_text(output / "stage24_manifest.json", json.dumps(manifest, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
