# Candidate generation results

这里保存候选药物生成和筛选阶段的结果，不保存候选药物后续的基因程序、stress、靶点和外部验证分析。

主要输入来自 `results/signatures/GSE179044_cmap_query_signatures.csv` 和本地 GSE92742/GSE70138 LINCS Level 5 数据。主要输出包括 LINCS ranking、context/perturbation summary、cross-release recurrence、positive-control sanity check，以及 92 行/29 个候选药物表。

候选表之后的分析位于 `data/processed/candidate_analysis/`，阶段说明和机制假说位于 `candidate_analysis/`，正式阶段报告位于 `reports/`。
