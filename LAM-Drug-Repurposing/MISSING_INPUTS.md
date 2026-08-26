# 当前缺口与最小补充方案

本项目已经把可以在本地自行完成的部分先推进了。剩余缺口主要来自外部扰动数据库或作者处理结果，单靠当前本地文件无法严谨补出。

## 已经自行补上的部分

- GSE104335：从完整 `GSE104335_RAW.tar` 中抽取了 9 个基因表达 CEL 和 9 个已处理的 HTA-2_0 `sst-rma-gene-full` CHP 文件，直接解析 CHP、完成 HTA 2.0 gene mapping 和 limma contrasts，脚本为 `scripts/analyze_gse104335_chp.R`。
- GSE302356：根据正式论文和可访问的预印本正文，建立了可追溯的操作性状态面板 `data/processed/GSE302356/paper_state_marker_panels.csv`，并将其加入原始矩阵评分。它不是作者发布的正式 `cell_id → state` 标签，LAMCORE3 尤其只有 shared-core/translation surrogate，不能当成唯一标记。
- GSE302356：当前 4 个已下载样本已经同时输出 residual、功能模块和五类状态面板的 cell/spot score。

## 仍然真正缺少的输入

### 0. CLUE API（已停止，接入层不再作为执行入口）

CLUE 官方已公告 clue.io 网站及其工具自 2026-01-31 退役。因此 `scripts/clue_api.py` 和 `results/cmap/clue_query_plan.json` 仅保留为历史记录，不再需要 API key，也不再把在线 CLUE 查询作为研究依赖。当前改用 GEO 中的公开 LINCS/CMap Level 5 数据做本地计算。

### 1. LINCS/CMap 原始 Level 5 数据（已下载并校验）

GSE92742 和 GSE70138 的 Level 5 矩阵及配套 metadata 已移动到：

```text
data/raw/LINCS/GSE92742/
data/raw/LINCS/GSE70138/
```

文件通过 GEO 提供的 SHA512SUMS 和 gzip integrity 检查。详细路径、大小和 SHA512 记录在 `manifests/LINCS_download_manifest.json`。目前尚未将 GCTX 矩阵转换为候选过滤脚本需要的长表，因此下一步是本地读取 Level 5、整理 perturbation metadata，并生成：

最小可用文件放在：

```text
data/processed/LINCS/perturbation_signatures.csv
```

至少需要三列：

```text
perturbation,gene,perturbation_score
```

其中 `perturbation_score` 应是药物/基因扰动对该基因的有方向效应；药物、target KD/KO 最好保留在同一文件并用不同的 `perturbation` 名称区分。若能同时提供以下三个小表，筛选会更可靠：

```text
data/processed/LINCS/target_perturbation_concordance.csv
data/processed/LINCS/generic_cytotoxicity.csv
data/processed/LINCS/exposure_feasibility.csv
```

它们的最低列分别是：

```text
perturbation,target_concordance
perturbation,generic_cytotoxicity_score
perturbation,exposure_feasible
```

原始数据已经完整下载；不需要再从 CLUE 下载或提供 API key。若后续本地处理遇到磁盘或 GCTX 读取限制，再考虑构建按 perturbation 子集的派生文件，但不改变原始文件。

### 2. GSE302356 作者状态标签或补充 marker 表

当前操作性面板可以继续做探索，但若要把结论升级为“富集于 LAMCORE2”或“来自 LAF-seed”，最好补充作者的 processed `h5ad`/Seurat 对象、cell metadata，或 supplementary marker table。需要至少保留：

```text
cell_id,state,donor_id,modality
```

如果只能下载一个文件，优先下载作者处理后的 cell metadata/marker table，而不是再次下载完整原始多组学归档。正式论文 DOI 为 [10.1183/13993003.02049-2025](https://doi.org/10.1183/13993003.02049-2025)。

### 3. GSE104335 已经解决

GSE104335 的归档实际包含 `HTA-2_0` 的 `sst-rma-gene-full` CHP 文件。项目已直接解析这 9 个已处理的 gene-level CHP、映射 23,910 个基因，并完成 limma 对比；不再需要用户补充文件。由于该数据没有 `shSRPK2 + rapamycin` 组合组，它仍然只能做 SRPK2/rapamycin 机制比较，不能单独证明联合治疗效应。

## 用户最值得提供的东西

按优先级：

1. 一份 LINCS/CMap 精简结果表，至少覆盖药物扰动和对应 target KD/KO；
2. GSE302356 作者的 state metadata 或 supplementary marker table；
3. 如果不希望下载 LINCS，提供一个候选药物清单及靶点、人体暴露/常用浓度信息，我可以先做机制互补、状态表达和 generic cytotoxicity 的人工可追溯筛选。

在这些输入到位前，当前研究仍可以继续做发现、机制解释和人类程序映射，但不会把探索性 score 误写成正式状态复现，也不会强行生成 0–5 个联合药物。
