# 研究线索卡 04：LAMCORE–protease 空间生态位

## 当前观察

在 GSE302356 的 Visium、Visium HD 和 Xenium 中，LAMCORE-like module score 较高的空间单位，与 protease module score 呈正向邻域富集。LAMCORE–protease 的单位级 Spearman 相关约为：Visium 0.818、Visium HD 0.402、Xenium 0.895；top-decile 邻域富集比约为 6.66、3.96 和 7.44。

## 解释边界

这是跨三种空间技术方向一致的候选空间关系，不是已证实的细胞通信，也不能证明蛋白酶由 LAMCORE 细胞本身产生。Visium 是 spot-level，Visium HD 是 bin/segment-level，Xenium 是 targeted cell-level；原始单位、score 和 p 值没有合并。Xenium 未检测到某个 antiprotease 基因只能记为 panel-unobserved。

## 可能机制

LAMCORE、免疫细胞、成纤维细胞和 ECM 可能共同形成局部蛋白水解/基质重塑生态位，参与囊性组织破坏。需要进一步拆解具体来源：LAMCORE 相关 CTSK/MMP、免疫细胞蛋白酶、LAF/成纤维细胞 ECM，以及 antiprotease 是否空间分离。

## 稳健性与证据等级

- PatientID：LAM4（Visium LAM20 + Xenium LAM19）和 LAM3（Visium HD LAM18）；不是两个独立患者验证。
- 技术证据：三种空间技术方向一致，属于正交技术支持。
- QC 门：空间数据不直接等同于 140→85 的单细胞 QC 门；后续若生成细胞级 Hypothesis Card，仍需进行 baseline/strict-QC 检查。
- 当前等级：探索性假说，尚不是高可信新程序。

## 可检验预测

1. LAMCORE-rich 区域应同时出现 MMP/CTSK 活性和 ECM 断裂或重塑标记；
2. protease 来源应可分解为 LAMCORE、免疫细胞和/或 LAF 贡献，而不是单一细胞来源；
3. antiprotease 与 protease 的空间关系可能比单纯表达均值更能区分病灶区域。

## 文件

`results/spatial/GSE302356/` 下的三个 `*_unit_scores.csv`、三个 `*_co_localization.csv` 和各自 manifest。
