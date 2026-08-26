# LAMCORE 未知状态程序发现报告：GSE135851_core_reproduction

> 本报告记录的是 GSE135851_core_reproduction 中的程序发现候选，不是独立验证完成后的新生物学结论。

运行模式：`fast_smoke`。当前参数用于流程验收，不作为最终模块数和稳定性结论。

## 当前完成内容

- 输入：30708 个细胞、63677 个基因；
- 候选池：高置信 535，宽松 3940，unrestricted guardrail 23228；
- 同时运行 pooled NMF 与 donor-wise 独立 NMF，并用核心基因重叠建立初步 meta-program；
- 已知程序只做事后比较，没有在主分析前回归掉；
- CORE3 使用 identity、深度校正后的低活性和 translation enrichment 三部分评分。

## 当前结果如何解释

目前所有程序仍是候选。外部 GSE190260、GSE217108 和 GSE302356 已转换为可分析 AnnData；本报告是快速 smoke run，不能替代完整参数和独立性分析。 即使外部数据已可运行，也不能仅凭 smoke run 计算最终独立验证等级。

## 运行摘要

- `broad_lam_like`：选择 rank=3、seed=1，factor stability=1.000。
- `high_confidence`：选择 rank=3、seed=0，factor stability=1.000。
- `unrestricted_lam`：选择 rank=3、seed=0，factor stability=1.000。

## 下一步

1. 下载并检查外部公开处理后矩阵；
2. 按 PatientID 进行 donor-wise 独立发现和 meta-program matching；
3. 对跨 donor 候选做 doublet、assay、深度、已知程序和 leave-one-donor-out 验证；
4. 用 GSE217108 的 ATAC、GSE302356 的 ATAC/空间和蛋白数据提高证据等级；
5. 只有在身份、独立 donor 和正交证据同时支持时，才升级为高可信研究线索。
