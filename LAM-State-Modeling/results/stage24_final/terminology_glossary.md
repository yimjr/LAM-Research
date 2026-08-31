# 术语和结论演变

| 术语 | 本项目最终使用方式 | 演变/边界 |
|---|---|---|
| `pool_high_confidence` | Stage 6–23 主线 `lam_candidate`，固定 5,378 cells | 不能与 broad/unrestricted 做并集；数量不足时不回退 |
| `pool_broad_lam_like` | upstream annotation；与 high 不重叠部分构成 boundary | 只作连续性/边界辅助，不自动成为 LAM |
| `pool_unrestricted_lam` | 审计标签 | 等同所有 condition=LAM 细胞，不能进入默认模型 cohort |
| candidate | 高召回、可供 State Modeling 检查的输入集合 | 不等于已证实 LAM；Stage15 证明 marker-combo gate 特异性不足 |
| consensus state | Stage7 的 20 个 frozen state labels | 不等同 Stage6 的 12 个 LAM-only grid clusters，也不等同 full-cohort 33 clusters |
| State15 | 200-cell frozen LAM-rich candidate state | 从“待验证 anchor”变成“provisional reference-anchor candidate”，未升级为正式 classifier |
| State16 | State15 邻近的 396-cell frozen state | 修正后仅保留为局部几何邻接/混合状态候选；不再标为 `LAM_to_lineage_transition_candidate` |
| manifold | Stage20 提出的全局 State15-centered gradient hypothesis | Stage21 削弱为非稳健统一 manifold；Stage22 保留局部分支 candidate |
| ordinary lineage adjacency | 与 State15 几何相邻但 matched-null/gradient 不支持 LAM branch | 当前用于 State16、State12、State20、State7；State16 虽有患者内 slope/LOPO 方向一致的探索性特点，但未通过方向匹配的 left-tail q 检验 |
| author-style | 上游真实作者逐细胞标签 | 只有 GSE135851 `available`；其余三个 dataset 是 `not_assayed` |
| formal LAMCORE | 777-gene formal signature score | Stage16 运行时 unavailable，Stage18 后在 data-temp 可用；不能倒写 Stage16 |
| structural stability | LOO/partition recovery 等结构维度 | 与 biological reproducibility 分开报告 |
| biological reproducibility | patient-aware DE/program/profile 的复现维度 | 不由 cluster 数或单个 p 值替代 |
| human cell analogue | 基于 marker/program/latent/normal 对照的解释 | 因 upstream `cell_type` 全 unknown，不是已验证人体细胞注释 |
| evidence vs hypothesis | 直接结果、解释、后续 candidate 分开 | Stage24 不把 provisional/uncertain 升级为事实 |

## 结论状态变化

`Stage 6 GO` → `Stage 7 20-state consensus` → `Stage 15 gate-specificity problem` → `Stage 18 State15 provisional anchor` → `Stage 20 global manifold hypothesis` → `Stage 21 non-robust global manifold` → `Stage 22 corrected local branch audit; State16 transition label withdrawn`.

这条链不是矛盾需要消除，而是研究问题在新证据下被逐步收窄。
