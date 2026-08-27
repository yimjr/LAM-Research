# 阶段报告

这里保存各研究环节的稳定总结，不替代 `research_log/` 中按时间记录的观察和新想法。

每份报告尽量说明：研究问题、数据集、输入文件、执行脚本、方法、输出、当前发现、限制和下一步。

## 推荐阅读顺序

```text
01 discovery
  → 02 external validation
  → 03 mechanisms / 04 human mapping
  → 05 candidate generation（研究转向：直接匹配 TSC2-loss 状态）
  → 06 candidate analysis（66 个药物及其 gene program/靶点）
  → 07 translation analysis（GSE277844 后续验证）
```

其中 01–04 主要解释疾病状态和 residual 的来源，05 负责生成候选，06 负责解释候选，07 是独立的 translation–residual 验证线。各线的结果不能直接合并为已确认的联合治疗证据。

```text
01_discovery              GSE179044 核心发现
02_external_validation    GSE27982/GSE16944 外部支持
03_mechanisms             STAT3/SRPK2 机制
04_human_mapping          人体 LAM 状态映射
05_candidate_generation   LINCS 候选生成
06_candidate_analysis     候选表之后的深入分析
07_translation_analysis   翻译程序与 residual 分开比较
```
