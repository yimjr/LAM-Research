# 最终项目摘要

本项目最终支持的最稳妥结论是：在一个高召回、包含明显普通肺谱系的 LAM candidate pool 中，State15 是当前最 LAM-rich 的 frozen consensus state，并具有 formal LAMCORE、可用数据集中的 author-label enrichment、患者匹配及去除 LAM1163 后仍保留的 profile。它尚未达到独立、均衡跨患者的正式 reference anchor 标准。

Stage20 的 pooled State15-centered gradient 在 Stage21 被限定为“存在部分独立 LAM-rich gradient，但不是稳健统一 global manifold”。Stage22 修正 branch eligibility、local scope、距离匹配、经验尾部和多重比较后，当前选中 State_12, State_16, State_20, State_7；State16 的当前标签为 ordinary_lineage_adjacency（方向性 left-tail p=0.439122，left-tail BH q=1.000000；two-sided p=0.878244，two-sided BH q=0.878244）。因此不再把 State15→State16 称为 transition candidate。没有证明时间转化、诊断 classifier 或所有 candidate state 都是 LAM。

完整材料见 `final_project_source_materials.md`；数字和文件来源见 `artifact_index.csv`、`stage_index.csv` 与 `narrative_audit.csv`。
