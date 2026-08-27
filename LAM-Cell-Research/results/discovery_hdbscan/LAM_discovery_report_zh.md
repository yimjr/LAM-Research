# LAMCORE 新生物学探索（Phase 3）

## 结论定位

本阶段从已经建立的作者风格 marker 候选出发，分析连续表达程序、LAM2 差异和微环境关联。结果用于生成研究线索，不把算法结构直接称为新亚型或真实细胞通信。

候选细胞数：140；候选供体：LAM1, LAM2, LAM3, LAM4。

## HDBSCAN 结果

主设定在 140 个候选细胞中识别出 2 个密度簇，分配 28 个细胞；112 个细胞被标记为噪声，覆盖 2 个 donor。
主设定的非噪声细胞只来自 LAM3 和 LAM4；LAM1 和 LAM2 均未进入密度簇。该结果不支持当前 140 个候选细胞中存在覆盖四个 donor 的稳定离散状态；应继续把状态视为连续程序或局部状态结构，并结合 donor-wise 分析验证。

## 已生成的主要结果

- `lamcore_state_programs_by_donor.csv`：每个 donor 的状态程序差值；
- `lamcore_state_heterogeneity_hdbscan.csv`：HDBSCAN 参数敏感性分析；
- `lamcore_state_hdbscan_cluster_summary.csv`：HDBSCAN 聚类摘要；
- `candidate_microenvironment_associations.csv`：候选表达关联，不等同真实通信；
- `external_validation_status.json`：GSE122960 正常肺和 GSE118180 小鼠子宫的处理状态；
- `results/hypothesis_cards_hdbscan/`：三张中英文研究线索卡。

当前最值得继续追踪的是 LAMCORE 内部连续状态，其次是 LAM2 的弱/异质信号以及与淋巴管和 ECM 的候选关联。独立 LAM donor 验证仍未完成。
