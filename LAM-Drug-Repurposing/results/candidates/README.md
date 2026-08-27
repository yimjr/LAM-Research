# Candidate generation results

这里保存候选药物生成和筛选阶段的结果，不保存候选药物后续的基因程序、stress、靶点和外部验证分析。

候选生成经历过一个重要调整：最初考虑以 sirolimus residual 的共同 reversal 作为入口，但 plastic 与 hydrogel 之间稳定共同的 reversal 方向很少，因此当前候选不要求两个环境先取 residual reversal 交集。主表改为保留 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 任一 TSC2-loss disease signature 命中的 concordant compound，再交给 `candidate_analysis/` 做后续解释。

主要输入来自 `results/signatures/GSE179044_cmap_query_signatures.csv` 和本地 GSE92742/GSE70138 LINCS Level 5 数据。主要输出包括 LINCS ranking、context/perturbation summary、cross-release recurrence、positive-control sanity check，以及按 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 筛选的 258 行/66 个候选药物表。原有 92 行/29 个 `tsc2_loss_plastic` 单环境结果保留为历史子集。

候选表之后的分析位于 `data/processed/candidate_analysis/`，阶段说明和机制假说位于 `candidate_analysis/`，正式阶段报告位于 `reports/`。
