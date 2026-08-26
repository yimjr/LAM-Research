# GSE84476 与 GSE104335：机制支持

## GSE84476：STAT3

- 脚本：`scripts/analyze_mechanisms.py`；
- 输出：`results/mechanisms/GSE84476_gene_level_log2_tpm_contrasts.csv`；
- 任务：比较 STAT3 knockdown、rapamycin 及其上下文差异。

STAT3 target engagement 和 rapamycin response 已在 gene level 整理，但该数据集主要承担机制支持角色，不是 TSC2 × rapamycin 的结构性复现。

## GSE104335：SRPK2

- 脚本：`scripts/analyze_gse104335_chp.R`、`scripts/analyze_mechanisms.py`；
- 输出：`results/mechanisms/GSE104335_gene_level_contrasts.csv`、`GSE104335_cross_dataset_summary.csv` 和 `GSE104335_hydrogel_specific_overlap.csv`；
- 任务：比较 SRPK2 knockdown 与 rapamycin 对 ECM、代谢和 stress 轴的影响。

NNMT、COL8A1、DCN、FBLN5 等结果提示 SRPK2 可能接触 GSE179044 的 ECM–metabolic axis，但由于没有 `shSRPK2 + rapamycin` 组，不能从该数据估计联用 interaction。

## 解释边界

STAT3 和 SRPK2 结果用于机制解释和候选优先级判断，不替代完整因子设计或独立药理/遗传验证。
