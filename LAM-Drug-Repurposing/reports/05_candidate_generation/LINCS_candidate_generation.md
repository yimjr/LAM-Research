# LINCS：候选药物生成

## 本阶段的研究转折

项目最初希望从 sirolimus/rapamycin 后仍存在的 residual 中寻找 reversal 药物，并要求 plastic 与 hydrogel 的结果能够共同支持同一方向。实际比较发现，两个环境之间稳定、共同的 reversal 方向很少，因此不再把共同 residual reversal 作为候选生成的必要条件。

本阶段改为分别使用 `tsc2_loss_plastic` 和 `tsc2_loss_hydrogel` 的 TSC2-loss disease signature 做本地 LINCS 匹配。这里的 TSC2-loss signature 是 WT 与 TSC2-null 的整体表达差异，而不是只使用 TSC2 单个基因。sirolimus/rapamycin 仍保留为背景和机制参照，residual/reversal 结果则在候选生成之后用于解释药物作用程序。

## 输入

- 疾病签名：`results/signatures/GSE179044_cmap_query_signatures.csv`；
- GSE92742/GSE70138 Level 5 GCTX：`data/processed/LINCS/gctx/`；
- 两个 release 的 signature、gene 和 perturbation metadata：`data/raw/LINCS/`。

## 方法

脚本：`scripts/analyze_lincs_local.py`。

每个 LINCS signature 在完整 BING gene-space 上排序，然后计算 weighted KS ES、WTCS、reversal_WTCS、published NCS 和 weighted-correlation sensitivity analysis。GSE92742 与 GSE70138 分别分析，之后只做 cross-phase/cross-release recurrence 比较，不称为独立生物学复现。候选筛选允许 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 任一环境进入，不要求两个环境先取共同 reversal 集合。

相同实际 gene set 的 query-size 会去重；environment-dependent escape 因可用基因很少，只作为 exploratory query。

## 输出

主要输出位于 `results/candidates/`，包括 candidate ranking、context/perturbation summary、recurrence、positive-control sanity check，以及按 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 筛选的 258 行/66 个候选药物表。原有 92 行/29 个 `tsc2_loss_plastic` 单环境结果保留为历史子集，便于比较候选范围扩大前后的变化。

## 当前结论边界

LINCS 结果用于产生候选和机制线索，不等于治疗建议。generic cytotoxicity、人体状态、target perturbation、人体暴露和外部实验模型证据仍需继续整合，因此当前不生成 Tier 1 联合治疗候选。
