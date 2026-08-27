# 当前缺口与最小补充方案

本项目已经把可以在本地自行完成的部分先推进了。剩余缺口主要来自外部扰动数据库或作者处理结果，单靠当前本地文件无法严谨补出。

## 已经自行补上的部分

- GSE104335：从完整 `GSE104335_RAW.tar` 中抽取了 9 个基因表达 CEL 和 9 个已处理的 HTA-2_0 `sst-rma-gene-full` CHP 文件，直接解析 CHP、完成 HTA 2.0 gene mapping 和 limma contrasts，脚本为 `scripts/analyze_gse104335_chp.R`。
- GSE302356：根据正式论文和可访问的预印本正文，建立了可追溯的操作性状态面板 `data/processed/GSE302356/paper_state_marker_panels.csv`，并将其加入原始矩阵评分。它不是作者发布的正式 `cell_id → state` 标签，LAMCORE3 尤其只有 shared-core/translation surrogate，不能当成唯一标记。
- GSE302356：当前 4 个已下载样本已经同时输出 residual、功能模块和五类状态面板的 cell/spot score。
- LINCS/CMap：GSE92742/GSE70138 Level 5 GCTX 与配套 metadata 已下载、校验并完成本地 WTCS/NCS、候选排名和 cross-phase/cross-release recurrence 分析。当前主候选为 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 任一环境命中的 258 条记录/66 个去重药物；不再依赖在线查询服务或 API。

## 仍然真正缺少的输入

### 1. GSE302356 作者状态标签或补充 marker 表

当前操作性面板可以继续做探索，但若要把结论升级为“富集于 LAMCORE2”或“来自 LAF-seed”，最好补充作者的 processed `h5ad`/Seurat 对象、cell metadata，或 supplementary marker table。需要至少保留：

```text
cell_id,state,donor_id,modality
```

如果只能下载一个文件，优先下载作者处理后的 cell metadata/marker table，而不是再次下载完整原始多组学归档。正式论文 DOI 为 [10.1183/13993003.02049-2025](https://doi.org/10.1183/13993003.02049-2025)。

## 已经解决的历史缺口

### GSE104335 已经解决

GSE104335 的归档实际包含 `HTA-2_0` 的 `sst-rma-gene-full` CHP 文件。项目已直接解析这 9 个已处理的 gene-level CHP、映射 23,910 个基因，并完成 limma 对比；不再需要用户补充文件。由于该数据没有 `shSRPK2 + rapamycin` 组合组，它仍然只能做 SRPK2/rapamycin 机制比较，不能单独证明联合治疗效应。

## 用户最值得提供的东西

目前最有价值的补充输入是：

1. GSE302356 作者的 state metadata 或 supplementary marker table；
2. RET、BTK、NTRK3、MKNK1/2 等方向的选择性药理或 CRISPR/shRNA 扰动结果；
3. 候选药物的人体暴露、有效浓度和 generic cytotoxicity 资料。

在这些输入到位前，当前研究仍可以继续做发现、机制解释和人类程序映射，但不会把探索性 score 误写成正式状态复现，也不会强行生成 0–5 个联合药物。
