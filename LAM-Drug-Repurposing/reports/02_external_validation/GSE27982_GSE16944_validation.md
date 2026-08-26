# GSE27982 与 GSE16944：外部支持

## 研究问题

GSE179044 的 TSC2-loss、rescue/persistence 和 rapamycin-insensitive 程序能否在其他实验中得到支持？

## GSE27982

- 数据：mouse MEF 的完整 Tsc2 × rapamycin 2×2 设计；
- 脚本：`scripts/analyze_external.py`；
- 输出：`results/tables/GSE27982_external_response.csv`；
- 主要任务：验证 TSC2-loss、rapamycin 后 residual/persistence 和 genotype-dependent response。

TSC2-loss 的跨模型方向一致性强于普遍 residual 程序。由于该实验处于低血清条件，G×R 不能自动解释为 escape。

## GSE16944

- 数据：LAM-like/TSC2 模型；
- 脚本：`scripts/analyze_support.py`；
- 输出：`results/tables/GSE16944_historical_support.csv` 和 `GSE16944_module_support.csv`；
- 主要任务：MMP2、ECM/invasion 等经典 rapamycin-insensitive program 的历史支持。

该数据缺少 TSC2-restored + rapamycin，因此不能计算正式的 KO_rapamycin − WT_rapamycin 或完整 G×R。

## 结论边界

这些数据提供外部支持，但不构成同一实验设计下的完整 residual/escape 复现。MMP2 更适合作为正交 assay control 和机制线索。
