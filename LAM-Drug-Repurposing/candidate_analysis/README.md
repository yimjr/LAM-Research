# Candidate analysis stage

本目录表示“候选药物生成之后的研究环节”，不是代码目录，也不是大型数据目录。

## 本阶段在整个项目中的位置

候选生成最初希望寻找能够同时 reversal plastic/hydrogel residual 的药物，但两个环境之间稳定共同的 reversal 方向很少，因此没有把共同 residual reversal 设为候选生成的硬门槛。本阶段承接的是直接基于 TSC2-loss disease signature 得到的候选集合：药物在 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 中命中即可进入当前后续分析，再结合两个 LINCS release 的 recurrence、基因程序和靶点证据进行解释。

这意味着本阶段要回答的是：

> 这些药物为什么能够匹配 TSC2-loss 状态？它们共同影响哪些基因/模块？这些模块是特异的 TSC2-loss 关联程序，还是 generic stress/cytotoxicity？药物靶点是否能提出可检验的机制？

## 输入

候选生成结果来自 `results/candidates/`，主要包括：

- `LINCS_candidate_ranking.csv`；
- `LINCS_cross_dataset_recurrence.csv`；
- `tsc2_loss_plastic_or_hydrogel_replicated_concordant_compounds.csv`；
- 当前 258 行候选记录；
- 当前 66 个去重候选药物。

原有的 `tsc2_loss_plastic` 单环境 92 行/29 个药物结果保留为历史子集，不作为合并 plastic/hydrogel 筛选的当前主表。

当前合并候选的靶点与 LINCS 基因程序结果已由独立新脚本生成，分别写入 `data/processed/candidate_analysis/drug_targets/`、`data/processed/candidate_analysis/programs/`、`data/processed/candidate_analysis/audit/` 和 `data/processed/candidate_analysis/validation/`，文件名以 `tsc2_loss_plastic_or_hydrogel_replicated_concordant_` 开头。原 29 药物的同类文件不受影响。

目前已分别使用 `tsc2_loss_plastic` 和 `tsc2_loss_hydrogel` 疾病面板分析同一批 66 个药物，并生成独立的 plastic-vs-hydrogel 比较结果。比较重点是两个 panel 各自的 reversal gene set、交集、Jaccard 和 panel-specific reversal genes；共同基因上的 LINCS effect 来自同一套 perturbation，只作为结构性可比性诊断。

## 分析顺序

```text
66 个 TSC2-loss 匹配候选药物
  → 每个 LINCS release 内稳健汇总 drug × gene response
  → reversal/mimic 基因集合和跨 release 稳定模块
  → generic stress / cytotoxicity 去卷积
  → ChEMBL / PubChem / BindingDB 靶点证据聚合
  → RET、BTK 等外部参考验证
  → 人体 LAM 状态映射
  → 机制假说和最小验证实验
```

## 文件职责

- `hypotheses/`：候选分析阶段形成的机制假说卡；
- `data/processed/candidate_analysis/`：可复用的结构化中间数据；
- `reports/06_candidate_analysis/`：稳定、可阅读的阶段报告；
- `scripts/`：所有执行代码；
- `manifests/`：输入、参数、输出和分析审计记录。

候选表是本阶段的输入，不是本阶段已经证明的治疗结果。原有 92 行/29 药物是 `tsc2_loss_plastic` 单环境历史子集；当前主分析是 plastic/hydrogel 任一环境命中的 258 行/66 药物集合。任何候选都需要结合人体状态、generic cytotoxicity、遗传扰动、暴露和外部模型证据继续判断，因此本阶段不自动生成 Tier 1 联合治疗候选。
