# LINCS：候选药物生成

## 输入

- 疾病签名：`results/signatures/GSE179044_cmap_query_signatures.csv`；
- GSE92742/GSE70138 Level 5 GCTX：`data/processed/LINCS/gctx/`；
- 两个 release 的 signature、gene 和 perturbation metadata：`data/raw/LINCS/`。

## 方法

脚本：`scripts/analyze_lincs_local.py`。

每个 LINCS signature 在完整 BING gene-space 上排序，然后计算 weighted KS ES、WTCS、reversal_WTCS、published NCS 和 weighted-correlation sensitivity analysis。GSE92742 与 GSE70138 分别分析，之后只做 cross-phase/cross-release recurrence 比较，不称为独立生物学复现。

相同实际 gene set 的 query-size 会去重；environment-dependent escape 因可用基因很少，只作为 exploratory query。

## 输出

主要输出位于 `results/candidates/`，包括 candidate ranking、context/perturbation summary、recurrence、positive-control sanity check，以及后续 92 行/29 个候选药物表。

## 当前结论边界

LINCS 结果用于产生候选和机制线索，不等于治疗建议。generic cytotoxicity、人体状态、target perturbation、人体暴露和外部实验模型证据仍需继续整合，因此当前不生成 Tier 1 联合治疗候选。
