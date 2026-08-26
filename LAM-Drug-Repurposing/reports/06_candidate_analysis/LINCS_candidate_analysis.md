# 候选药物后续分析

## 输入

本阶段以 `results/candidates/` 中的候选生成结果为输入，重点使用 `LINCS_candidate_ranking.csv`、`LINCS_cross_dataset_recurrence.csv`、92 行候选表和 29 个去重药物表。

## 分析顺序

```text
候选药物
  → reversal/mimic 基因响应
  → 跨 release 稳定模块
  → generic stress/cytotoxicity 去卷积
  → PubChem/ChEMBL/BindingDB 靶点整理
  → RET/BTK 等外部参考
  → 人体 LAM 映射
  → 机制假说
```

## 当前观察

已得到 reversal-only 和 mimic-only 的跨 release 稳定模块。strongest reversal module 具有 glycolysis/hypoxia 与 ECM organization/remodeling 特征；mimic module 同时包含 stress/proteostasis 信号，因此需要先完成 generic stress 去卷积。

Lestaurtinib 的 RET/NTRK 轴和 QL-X-138 的 BTK/MNK 轴目前是待验证机制，不应根据单一全局相关性直接确定靶点。PubChem 与 BindingDB 作为证据层，不直接等同于主要作用靶点。

## 输出位置

- 结构化数据：`data/processed/candidate_analysis/`；
- 阶段报告：本目录；
- 机制假说：`candidate_analysis/hypotheses/`；
- 详细过程记录：`research_log/`。

## 下一步

继续完成 generic cytotoxicity 去卷积、稳定模块的人体映射，以及 RET、BTK、NTRK3 的选择性药理或 CRISPR/shRNA 验证。当前不生成 Tier 1 联合治疗候选。
