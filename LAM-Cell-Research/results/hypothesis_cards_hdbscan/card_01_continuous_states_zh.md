# LAMCORE 可能包含连续表达状态

## 分类：探索性方法结果（同一数据集，尚无独立 LAM donor 验证）

## 观察到了什么
在作者风格 marker 候选细胞中，连续状态程序在多个 LAM donor 中呈现方向一致的候选信号：contractile, ecm_remodeling, hormone_related, stress_hypoxia, proliferative。具体差值见 `lamcore_state_programs_by_donor.csv`。

## HDBSCAN 结果
主设定在 140 个候选细胞中识别出 2 个密度簇，分配 28 个细胞；112 个细胞被标记为噪声，覆盖 2 个 donor。当前主设定中的非噪声细胞只来自 LAM3 和 LAM4；LAM1 与 LAM2 没有被分配到密度簇。因此，HDBSCAN 没有支持一个覆盖四个 donor 的稳定离散状态；噪声标签表示状态空间中未形成足够密度的区域，不表示细胞没有生物学意义。

## donor、细胞和 pathway
单位是 LAM1–LAM4 donor；涉及 contractile、ECM remodeling、stress/hypoxia、inflammatory、hormone-related、metabolic 和 mTOR-related 程序。

## 稳健性
候选群本身来自已知 marker + 作者风格图聚类；Phase 2 已比较 doublet、QC、聚类种子/分辨率、777 module score、rank-based score、assay 分层和 leave-one-donor-out。由于同一批数据用于发现与评估，仍不能称独立验证。

## 替代解释
候选定义与部分 contractile/ECM marker 有重叠；组织处理、assay、细胞周期、应激和供体差异都可能贡献信号。

## 下一步验证
在独立 LAM donor、空间转录组、蛋白或 snATAC 数据中，先测试这些连续程序是否仍出现在相同 marker-defined LAMCORE-like 细胞中，再做 TSC2/mTOR、ECM 和淋巴管相关实验。

## 新颖性/可信度/优先级
新颖性：中；可信度：低到中；推荐优先级：中。当前结果更支持连续状态和局部高密度结构的进一步研究，不支持新的跨 donor 亚型结论。
