# Translation analysis stage

本目录表示一个独立的后续验证环节：在候选药物生成已转向直接匹配 TSC2-loss disease signature 之后，用 GSE277844 检查 TSC2 loss 相关的 translation abnormality 是否与 GSE179044 的 residual 程序相连，并检查翻译靶向处理是否使重叠基因的 KO-vs-WT 翻译差异向 WT 靠近。它不是重新生成 LINCS 候选，也不把 plastic residual 与 hydrogel residual 强行合并。

## 研究问题

```text
GSE277844
  → TSC2-loss translation-up/down genes
  → 分别比较 ordinary/plastic residual 与 hydrogel residual
  → 检查 RMC-6272、eFT-508 是否降低重叠基因的翻译差异
```

普通 residual 与 hydrogel residual 是两条独立比较线。分析不要求它们先取交集；只有在需要描述共同成员时，才把交集作为辅助信息。

## 输入

- `data/raw/GSE277844/GSE277844_raw_counts.txt.gz`：GSE277844 的 total 与 polysome 原始 counts。
- `results/tables/GSE179044_factorial_contrasts.csv`：GSE179044 的 TSC2-loss、rapamycin residual 和环境相关 contrast。

原始文件不随版本库提供。远端重建时可先运行：

```bash
./.venv/bin/python scripts/download_resumable.py \\
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277844/suppl/GSE277844_raw_counts.txt.gz \\
  data/raw/GSE277844/GSE277844_raw_counts.txt.gz
./.venv/bin/python scripts/analyze_gse277844_translation_residuals.py
```

GSE277844 是人源神经祖细胞的 TSC2 等基因模型，和 LAM 培养模型并非同一实验体系。因此这里得到的是跨模型的 translation–residual 程序联系，不是 LAM 的直接复现。

## 分析顺序

1. 从匹配的 polysome/total 样本估计 translation efficiency，并提取 TSC2-null 相对 WT 的 translation-up/down 基因。
2. 将这些基因分别与以下 residual 集合比较：plastic persistent residual、hydrogel persistent residual、rapamycin 条件下的 hydrogel residual，以及 hydrogel-specific residual。
3. 对每个独立 residual 类别分别汇总 RMC-6272 与 eFT-508 的处理效果。
4. 将 residual overlap 的恢复率与全部 translation-abnormal genes 的恢复率比较；同时保留去掉 overlap 后的非重叠背景作为敏感性比较。
5. 逐基因比较 RMC-6272 与 eFT-508，区分两种药都恢复、仅一种药恢复、两种药都未恢复的基因。
6. 以 `distance_to_WT_reduced`、signed residual ratio、方向和处理后效应量作为探索性证据，形成后续机制假说。

对固定的 13 个共同恢复基因，功能注释再单独使用 GO Biological Process、GO Cellular Component、Reactome 和 MSigDB Hallmark。输出同时保留逐基因落入的完整条目、主题级重复情况和带 FDR 的描述性富集，重点检查 ECM、actin/cytoskeleton、adhesion、Rho GTPase、focal adhesion、Hippo/YAP/TAZ、mechanotransduction、TGF-β 和 migration。

## 方法边界

- 翻译效应使用透明的 conditional polysome-vs-total model：`log2 CPM Poly ~ centered log2 CPM Total + genotype`。这与 anota2seq 的翻译成分相似，但不是官方 anota2seq 的完整复现。
- 默认使用 published-style FDR 0.15 选择 translation-up/down；效应量保留用于排序和解释。
- residual 的 `_q10` 集合另加绝对效应量 ≥ 0.5 与 moderated FDR ≤ 0.10。
- 处理后的“恢复”只表示 KO-vs-WT 翻译距离变小，不等于已经证明药物特异性机制或恢复为正常生理状态。
- 背景比较同时报告全部 selected translation genes 和带有绝对基线效应量门槛的子集；后者用于避免小基线效应导致的距离判定不稳定。Fisher 检验仅作小样本描述，不能把基因视为完全独立重复。
- RMC-6272 组重复数较少；所有药物结果均应视为探索性，需要基因/模块级和独立扰动验证。
- 13 基因功能注释中的关键词主题是解释性归纳，不等同于新的统计检验；四类基因集之间有大量重叠，不能把跨库重复简单相加。

## 文件位置

- 所有代码继续放在 `scripts/`，主脚本为 `scripts/analyze_gse277844_translation_residuals.py`。
- 可复用的结构化结果写入 `data/processed/translation_analysis/`。
- 背景恢复率比较写入 `GSE277844_translation_residual_recovery_background_comparison.csv`；逐基因药物一致性写入 `GSE277844_translation_residual_drug_gene_concordance.csv` 和 `GSE277844_translation_residual_drug_concordance_summary.csv`。
- 固定 13 基因的功能注释写入 `GSE277844_hydrogel_translation_core_*` 文件，阶段报告为 `reports/07_translation_analysis/GSE277844_hydrogel_translation_core_functional_annotation.md`。
- 稳定阶段报告写入 `reports/07_translation_analysis/`。
- 输入、参数、统计口径和输出记录写入 `manifests/GSE277844_translation_residual_analysis.json`。
- `hypotheses/` 保存本环节形成的机制假说卡；目前先保留目录说明，待 gene/module-level 解释完成后再写入具体假说。

远端重新获取项目后，需要重新运行相应脚本生成被 Git 忽略的运行结果。
