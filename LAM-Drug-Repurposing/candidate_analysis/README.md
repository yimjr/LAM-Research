# Candidate analysis stage

本目录表示“候选药物生成之后的研究环节”，不是代码目录，也不是大型数据目录。

## 输入

候选生成结果来自 `results/candidates/`，主要包括：

- `LINCS_candidate_ranking.csv`；
- `LINCS_cross_dataset_recurrence.csv`；
- 92 行候选记录；
- 29 个去重候选药物。

## 分析顺序

```text
候选表
  → 药物—基因响应和稳定模块
  → generic stress / cytotoxicity 去卷积
  → PubChem/ChEMBL/BindingDB 靶点整理
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

本阶段不自动生成 Tier 1 联合治疗候选。任何候选都需要结合人体状态、generic cytotoxicity、遗传扰动、暴露和外部模型证据继续判断。
