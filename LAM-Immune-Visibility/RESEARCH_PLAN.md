# LAMCORE 细胞状态依赖的免疫可见性

## 目标

比较 LAMCORE 的抗原相关表达、抗原加工/呈递、IFN 反应和免疫逃逸状态，检验这些状态是否与 stem-like/inflammatory 程序以及 CD8 T、NK、TREM2/TYROBP macrophage 状态相关。

## 数据范围

第一轮只读复用：

- GSE135851 核心肺部数据；
- GSE190260、GSE217108、GSE302356 外部数据；
- GSE122960 正常肺参考；
- GSE302356 空间矩阵；
- GSE179044 rapamycin retention 结果；
- 已有 LAM-Drug-Repurposing 人类状态映射结果。

不下载 FASTQ 或新的补充数据。

## 身份保护

主分析复用既有 candidate pool：`pool_high_confidence`；`pool_broad_lam_like` 和 `pool_unrestricted_lam` 只用于敏感性分析和 guardrail。PMEL、MLANA、TYRP1、DCT、MITF、GPNMB 等候选抗原不进入 identity-protected panel。

## 状态定义

非零原始 count 即表示“已检出”。低表达只表示已检出后处于预先定义的低表达区间：

```text
not_assayed    基因不在当前矩阵/panel
not_detected   基因存在，但原始 count = 0
detected_low   原始 count > 0，处于低表达区间
detected_high  原始 count > 0，处于高表达区间
```

低/高表达区间按 `dataset × assay × gene` 的非零表达分布计算，并写入运行 manifest。非零观测过少时只报告 `detected`，不强行划分 low/high。

## 首轮结果

1. 患者级 antigen/presentation/evasion 状态；
2. PMEL/MLANA 与 HLA/B2M/TAP 的共存关系；
3. stem-like/inflammatory 与免疫可见性模块的关联；
4. LAMCORE 与 CD8/NK/macrophage 状态的患者级对应；
5. GSE135851、GSE190260、GSE217108、GSE302356 的方向性重复；
6. rapamycin retention 与候选抗原排序。

