# LAMCORE 针对性稳健性验证

## 目的

第二阶段针对具体候选结果进行质量过滤，而不是把所有方法比较做成项目终点。Phase 1 候选标签固定后，分别比较 doublet 是否去除、QC 宽严、聚类种子/分辨率、777 module score、rank-based score、assay 分层和 leave-one-donor-out。

- Phase 1 不去除 doublet：30708 个细胞，候选 140 个。
- 去除预测 doublet：30077 个细胞，候选 140 个。
- 预测 doublet 总数：631。
- 宽松/严格 QC 候选数：140 / 85。

这些结果回答的是候选群对技术处理是否敏感，不等同于独立 donor 验证，也不设置机械的 3/4 donor 通过门槛。
