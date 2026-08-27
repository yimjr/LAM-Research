# 候选药物后续分析

本阶段承接的是直接匹配 TSC2-loss disease signature 得到的候选集合，而不是要求药物先同时 reversal plastic/hydrogel residual。两个环境之间共同、稳定的 reversal 方向很少，因此当前主集合允许药物在 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 任一环境中命中，再开展基因程序、靶点和机制分析。

当前主分析为 plastic/hydrogel 合并后的 66 个去重药物；原 `tsc2_loss_plastic` 单环境 29 药物分析作为历史子集保留。新增的 66 药物分析使用独立脚本和独立输出前缀，详见 `tsc2_loss_plastic_or_hydrogel_replicated_concordant_LINCS_gene_program_analysis.md`。
新增的 `tsc2_loss_hydrogel` 独立面板及其与 plastic 面板的比较，详见 `tsc2_loss_plastic_or_hydrogel_replicated_concordant_hydrogel_panel_LINCS_gene_program_analysis.md` 和 `tsc2_loss_plastic_or_hydrogel_replicated_concordant_plastic_vs_hydrogel_comparison.md`。

## 输入

当前主分析以 `results/candidates/` 中的 258 行/66 药物候选结果为输入，重点使用 `LINCS_candidate_ranking.csv`、`LINCS_cross_dataset_recurrence.csv` 和 plastic/hydrogel 合并候选表；92 行/29 药物表只作为 `tsc2_loss_plastic` 单环境历史子集保留。

## 分析顺序

```text
TSC2-loss 匹配候选药物
  → 每个 release 内稳健汇总 drug × gene response
  → reversal/mimic 基因响应和稳定模块
  → generic stress/cytotoxicity 去卷积
  → ChEMBL/PubChem/BindingDB 靶点整理
  → RET/BTK 等外部参考
  → 人体 LAM 映射
  → 机制假说
```

## 当前观察

已得到 reversal-only 和 mimic-only 的跨 release 稳定模块。strongest reversal module 具有 glycolysis/hypoxia 与 ECM organization/remodeling 特征；mimic module 同时包含 stress/proteostasis 信号，因此需要先完成 generic stress 去卷积。这里的 reversal/mimic 是对 TSC2-loss disease panel 的药物映射方向，不能直接等同于“逆转了 sirolimus residual”。

Lestaurtinib 的 RET/NTRK 轴和 QL-X-138 的 BTK/MNK 轴目前是待验证机制，不应根据单一全局相关性直接确定靶点。PubChem 与 BindingDB 作为证据层，不直接等同于主要作用靶点。

## 输出位置

- 结构化数据：`data/processed/candidate_analysis/`；
- 阶段报告：本目录；
- 机制假说：`candidate_analysis/hypotheses/`；
- 详细过程记录：`research_log/`。

## 下一步

继续完成 generic cytotoxicity 去卷积、稳定模块的人体映射，以及 RET、BTK、NTRK3 的选择性药理或 CRISPR/shRNA 验证。当前不生成 Tier 1 联合治疗候选。
