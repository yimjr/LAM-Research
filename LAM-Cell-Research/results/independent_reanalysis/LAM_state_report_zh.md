# LAM 细胞状态研究报告

- 数据集：`GSE135851`
- 纳入分析的细胞数：**30,077**
- counts 数据契约中的基因数：**63,677**
- 正式 LAMCORE signature 基因数：**777**；在本矩阵中可匹配：**777**
- 正式 LAMCORE 候选识别方法：`cluster_mean_score`
- 正式 LAMCORE 候选细胞数：**2,305**
- 四基因探索版候选细胞数：**146**
- 两种候选集合的重叠细胞数：**96**

## 这次真正复现了什么

本轮从 LAM Cell Atlas 官方页面下载了 777 个不重复的 LAMCORE 基因，并在 GSE135851 细胞上用完整基因表重新计算分数。分数采用 Scanpy 的 control-gene module score，输入是标准化并 log1p 后的表达矩阵。因此，本轮复现了正式基因集合和透明的再评分过程，但不是对原论文细胞标签或分类阈值的逐字节复制。
ACTA2、PMEL、FIGF、MLANA 四基因只保留作探索性对照，不能替代正式 signature。

## 结果解释边界

本分析属于探索性研究；单细胞层面的统计不作为供体级证据。
只有四位 LAM 供体都满足预先设定的供体级规则，才能把某个状态称为跨供体稳定。正式 signature 候选细胞数达到最低要求的供体为 4/4；四基因对照为 3/4。由于候选细胞本身是用正式 LAMCORE signature 选出来的，与 signature 重叠的程序分数升高不构成独立验证；重叠关系记录在 `lamcore_signature_program_overlap.csv`。

## 正式 signature 的供体级表达程序一致性

```text
     程序  可比较LAM供体数  候选细胞较高的供体数  方向一致比例  通过可观察供体规则  通过严格3/4规则
收缩/平滑肌样          4           4    1.00       True       True
  ECM重塑          4           4    1.00       True       True
     增殖          4           1    0.25      False      False
  应激/缺氧          4           4    1.00       True       True
炎症/免疫响应          4           4    1.00       True       True
```

## 正式 signature 与探索版对照

LAM2 中正式 signature 候选细胞为 28 个，四基因探索版为 3 个；Donor1 中正式 signature 候选细胞为 1292 个。这说明候选定义必须明确报告，任何一种候选集合都不能直接当作细胞类型标签。

## 稳健性检查

阈值敏感性、宽松/主分析/严格 QC 敏感性以及按 assay 分层的汇总由 `scripts/robustness_tests.py` 生成。这些检查用于评估当前候选定义是否对参数敏感，不能替代独立队列验证。

## 结论

按照正式 777 基因候选规则，以下程序满足预设的供体方向 3/4 规则：收缩/平滑肌样、ECM重塑、应激/缺氧、炎症/免疫响应；未满足的程序：增殖。这说明在当前操作性定义下，部分表达程序具有可重复方向，但不能据此证明存在新的 LAMCORE 亚型。候选比例在供体之间差异很大，LAM2 只有 28 个正式候选，而且候选定义与部分程序存在 signature 基因重叠，因此仍需独立数据验证。
Guo 等人的原始研究使用其自己的分类流程，在 LAM2 中报告没有 LAMCORE 细胞。因此，本轮出现 28 个 LAM2 候选并不等于推翻原结果：本轮使用了发表的基因集合，但采用了新的、透明的 Scanpy 分数和 cluster 阈值，并没有复刻原始细胞级分类阈值。

## 数据与运行输出

- 正式 signature 文件：`data/raw/reference/LAM_core_signature_genes.csv`（SHA-256 `c96cf4c684e6aac0976775fc3ce4dc894e13d264007e035f29950bb95d7fe9d6`）
- 处理后的 AnnData：`data/processed/GSE135851_lam_states.h5ad`
- 结果表：`results/tables`
- 图形：`results/figures`
- 官方来源：https://research.cchmc.org/pbge/lunggens/lungDisease/lamcore_query.html
