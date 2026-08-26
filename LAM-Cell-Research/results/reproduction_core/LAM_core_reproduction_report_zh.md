# LAMCORE 核心肺部复现：Phase 1 基线

## 结论定位

本结果是依据作者公开 R 脚本建立的独立重实现，不是严格的逐版本 Seurat 复现。原因是作者脚本混用 Seurat 2/3，而当前主流程使用 AnnData/Python；所有差异记录在 `method_deviation_table.csv`。

- 作者报告的 LAMCORE 细胞数：约 125 个，来自 LAM1、LAM3、LAM4。
- 本项目作者风格 marker/cluster 重实现候选数：**140**。
- 各供体候选数：`{"LAM1": 31, "LAM2": 4, "LAM3": 84, "LAM4": 21}`。
- LAM2 候选数：**4**；该数值只能说明本操作性规则下的 marker 候选，不能直接等同为论文定义的 LAMCORE。

## 方法边界

候选细胞先由 PMEL、ACTA2、ESR1、FIGF/VEGFD、CTSK、MLANA 等已知特征结合作者风格图聚类定位；777-gene signature 只在候选确定以后作为一致性检查。因此没有用由原始 LAMCORE 细胞总结出的 777 genes 反过来定义它们。

QC 仅表示“在处理后矩阵允许范围内恢复下游 QC”；FASTQ、Cell Ranger、初始 barcode/cell calling 和 empty droplet 判断没有被恢复。doublet 分数和预测已记录，但 Phase 1 没有据此删除细胞。

## 下一步

核心基线已具备进入第二阶段针对性稳健性验证和第三阶段新发现探索的条件；两者并行推进。若后续 R/Seurat 运行能进一步缩小方法差异，则会更新为更严格的复现版本。
