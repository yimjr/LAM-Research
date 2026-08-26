# 研究线索卡 05：LAMCORE 的肺部获得性程序

## 当前观察

在保守的 LAM1–LAM4 LAMCORE 候选细胞中，肺适应程序相对于正常肺参考的 donor-level 中位数差值约为 0.667；相对于小鼠正常子宫谱系参考也高约 0.667。平滑肌谱系程序仍然存在，但比小鼠子宫参考高约 0.367，而 ECM interaction 程序低约 1.216。

## 初步解释

这与“LAMCORE 保留平滑肌/谱系特征，同时获得部分肺部适应和 ECM 交互特征”的模型相容。它提示肺 LAM 细胞可能不是简单复制子宫谱系，而是在肺环境中重新组合了原有谱系和肺部适应程序。

## 重要限制

正常肺是多个正常肺 donor 的混合参考；小鼠子宫只有一个参考对象，且这里使用保守的 gene-symbol overlap，不是完整的正式人鼠 ortholog 映射。该比较不能证明细胞起源、迁移方向或因果适应机制。

## 证据等级

探索性假说。需要独立 PatientID、细胞类型匹配的子宫/肺参考、ATAC 或空间定位来升级。

## 可检验预测

1. 肺 LAMCORE 中的肺适应程序应在独立肺 donor 中重复，而不应只由 LAM3/LAM4 驱动；
2. 肺适应程序的核心基因应与 ECM/淋巴管/缺氧环境存在可定位的空间关系；
3. 与平滑肌谱系相同的基因应表现为“保留”，而肺适应独有部分应在子宫谱系参考中较弱。

## 文件

`results/adaptation/adaptation_cell_scores.csv`、`adaptation_donor_program_scores.csv`、`adaptation_contrasts.csv`。
