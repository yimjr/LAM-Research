# GSE135851 与 GSE302356：人体 LAM 映射

## 研究问题

候选药物所对应的 residual 或稳定模块，是否存在于人体 LAMCORE、LAF 或 ECM niche 状态？

## 数据与方法

- GSE135851：`scripts/map_human_states.py`，结果位于 `results/human_mapping/`；
- GSE302356：`scripts/analyze_gse302356_raw.py` 和 `scripts/analyze_lincs_deconvolution_human.py`；
- 主要方法：rank/module score、样本或患者层面的 enrichment 和跨模态比较。

## 当前状态

GSE135851 当前快照主要包含 candidate/other 标签，因此只能作为初步人类 LAM-like 支持。GSE302356 已对部分 scRNA-seq、Visium HD 和 Visium 样本进行模块评分，并使用 paper-derived LAMCORE/LAF operational panel。

当前结果支持继续检查 LAMCORE1/2/3、LAF-seed/LAF-niche 与 ECM、stress、myogenic 和代谢模块的关系，但不能把不同样本、模态或患者自动解释成配对多组学验证。

## 限制

LAMCORE3 缺乏独特 marker，当前使用 shared-core + translation/low-activity surrogate。正式 cell label、配对 ATAC/空间 metadata 和完整状态注释仍是后续补充内容。
