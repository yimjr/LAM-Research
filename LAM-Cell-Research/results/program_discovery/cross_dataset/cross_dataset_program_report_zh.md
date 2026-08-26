# 跨数据集程序比较报告

> 本报告比较不同运行中候选程序的 top-50 gene overlap，不把 overlap 直接当作生物学复现。

## 当前结果

- 共比较 73 对程序；
- Jaccard ≥ 0.15 的匹配：2 对；
- 来自不同 PatientID 集合的强匹配：0 对；
- 含同一 PatientID 的匹配：2 对，主要涉及 GSE217108 与 GSE302356 的 LAM32。

当前没有发现达到该阈值、且来自不同 PatientID 集合的稳定 meta-program。这个结果只能说明当前候选定义、特征选择和 top-gene 匹配下没有强信号，不能证明不存在新的跨患者程序。

## 最高重叠匹配

| pool | 数据集 | 程序 | 数据集 | 程序 | Jaccard | PatientID 关系 |
|---|---|---|---|---|---:|---|
| broad_lam_like | GSE217108 | program_1 | GSE302356 | program_5 | 0.220 | same_patient_overlap_present |
| broad_lam_like | GSE217108 | program_4 | GSE302356 | program_4 | 0.220 | same_patient_overlap_present |
| high_confidence | GSE190260 | program_1 | GSE302356 | program_2 | 0.136 | different_patient_sets |
| high_confidence | GSE217108 | program_3 | GSE302356 | program_1 | 0.111 | same_patient_overlap_present |
| unrestricted_lam | GSE190260 | program_4 | GSE302356 | program_5 | 0.099 | different_patient_sets |
| broad_lam_like | GSE190260 | program_1 | GSE302356 | program_2 | 0.087 | different_patient_sets |
| unrestricted_lam | GSE217108 | program_6 | GSE302356 | program_4 | 0.075 | same_patient_overlap_present |
| broad_lam_like | GSE217108 | program_2 | GSE302356 | program_1 | 0.075 | same_patient_overlap_present |
| high_confidence | GSE217108 | program_5 | GSE302356 | program_5 | 0.064 | same_patient_overlap_present |
| broad_lam_like | GSE135851_core | program_2 | GSE190260 | program_3 | 0.053 | different_patient_sets |

## 下一步

1. 在 donor 内独立发现结果上进行更稳健的 meta-program matching，而不是只比较 pooled NMF；
2. 对候选程序进行 rank-based score、已知状态解释比例和 leave-one-donor-out 验证；
3. 使用 GSE217108 ATAC、GSE302356 ATAC/空间数据检查正交支持；
4. 将 LAM32 的同患者跨 assay 重复与真正不同患者的复现分开报告。
