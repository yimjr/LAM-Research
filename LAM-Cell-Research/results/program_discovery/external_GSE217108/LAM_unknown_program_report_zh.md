# LAMCORE 未知状态程序发现报告：GSE217108

> 本报告记录的是 GSE217108 中的程序发现候选，不是独立验证完成后的新生物学结论。

运行模式：`configured_full`。当前使用配置中的完整 rank、seed 和 donor-wise 分析参数。

## 当前完成内容

- 输入：12396 个细胞、36601 个基因；
- 候选池：高置信 1075，宽松 3855，unrestricted guardrail 12396；
- 同时运行 pooled NMF 与 donor-wise 独立 NMF，并用核心基因重叠建立初步 meta-program；
- 已知程序只做事后比较，没有在主分析前回归掉；
- CORE3 使用 identity、深度校正后的低活性和 translation enrichment 三部分评分。

## 当前结果如何解释

目前所有程序仍是候选。外部 GSE190260、GSE217108 和 GSE302356 已转换为可分析 AnnData，并已完成完整参数运行；结果仍需正交模态验证。 即使完整 RNA 分析已完成，也仍需正交模态和 PatientID 独立性检查。

## 运行摘要

- `high_confidence`：选择 rank=6、seed=0，factor stability=1.000。
- `broad_lam_like`：选择 rank=5、seed=2，factor stability=1.000。
- `unrestricted_lam`：选择 rank=6、seed=2，factor stability=1.000。

## 下一步

1. 下载并检查外部公开处理后矩阵；
2. 按 PatientID 进行 donor-wise 独立发现和 meta-program matching；
3. 对跨 donor 候选做 doublet、assay、深度、已知程序和 leave-one-donor-out 验证；
4. 用 GSE217108 的 ATAC、GSE302356 的 ATAC/空间和蛋白数据提高证据等级；
5. 只有在身份、独立 donor 和正交证据同时支持时，才升级为高可信研究线索。
