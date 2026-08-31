# 阶段 1–6 与 Step 7–13 研究边界

## Stage 24：最终项目材料包

Stage 24 是项目冻结后的只读整理阶段。开始后不修改既有 state、State15、candidate gate、scVI、branch 或任何 Stage 1–23 核心 artifact；不重新训练、不重新聚类、不扩展新的生物学问题。脚本 `24_finalize_project.py` 盘点 Stage 1–23 的脚本、reports、results、关键 AnnData/model、配置和声明的 upstream 路径，并在 `results/stage24_final/` 生成最终报告原材料。

输出包括 `stage24_manifest.json`、`stage_index.csv`、`stage_summary.md`、`artifact_index.csv`、20 个 state 的 `state_human_cell_analogue.csv/md`、`other_findings_registry.csv/md`、`narrative_audit.csv`、`terminology_glossary.md`、`final_project_source_materials.md`、最终报告摘要、局限/未来方向以及 methods、state atlas、statistics、provenance 附录。历史冲突按时间保留，不以 Stage24 覆盖旧报告。运行：`$PY 24_finalize_project.py`。

目标是判断跨数据集、跨患者的高置信度 LAM candidate 是否存在可重复的 latent disease-state structure。

本轮只完成：输入盘点与继承、外部 QC、统一 gene universe、PCA/UMAP/Leiden、State Modeling 专用 NMF、scVI、Stage 6 Go/No-Go。已有 NMF/program/state 仅在模型训练完成后进行 post-hoc interpretation。

`lam_candidate` 固定等于 `pool_high_confidence`；`pool_broad_lam_like AND NOT pool_high_confidence` 是 boundary；`pool_unrestricted_lam` 不得扩张核心 cohort。

Go 需要高置信度 candidate 内存在至少两个 latent clusters、至少一个 cluster 跨两个独立 patient 出现，且 patient/dataset/assay 不构成主要分离轴。No-Go 只包括无内部结构、结构只来自单一患者、或结构主要由 patient/dataset/assay 驱动。latent state 与旧状态不对应不属于 No-Go，而是 novel/unexplained 的后续研究候选。

## Step 7–13

Step 7–13 不重训 scVI，严格使用既有 `X_scVI` 的 5,378 个 high-confidence LAM cells。Step 7 对 9 个 configuration 等权：每个 configuration 内先平均 seed co-assignment，再平均 configuration；固定的 seed 为 20260829–20260833，`n_neighbors=30` 的三个 configuration 各有五个 seed，其余 configuration 各有一个基础 partition。最终 consensus 允许形成一份完整 float32 co-assignment 和 average-linkage hierarchy，不预设 consensus cluster 数。

Step 8 以 full-data `n_neighbors=30, resolution=0.4` 为 reference，先保存 reference→consensus correspondence；每个 patient 和 dataset 的 LOO 先定义 retained cell set，再将 full consensus、full reference 和 LOO clustering 都限制到该集合，计算 baseline、recovery 与 additional loss。matching 允许 split/merge。

Step 9 分析 consensus state 的 latent distance、PAGA/graph connectivity、grid split/merge 和 boundary transitions；parent/substate 仅作解释，不改写 consensus。Step 10 对每个 state 单独构造 `State_k` 与同患者 `Rest_of_LAM` 的 pseudobulk，拟合 `~ patient_id + group`；不建立所有 state 合并的多分类或统一二元模型。低于最少患者支持时只输出描述性 pseudobulk。

Step 11 输出 patient×state 证据以及连续的 structural stability、biological reproducibility。Step 12 单独评估 boundary 与 optional normal reference，不把它们加入 state 数量。Step 13 汇总 atlas 与 hypothesis candidates，不设置硬性 confidence 门槛，也不以 cluster 数量或单个 DE p 值替代患者级复现证据。全流程单进程、顺序运行，并保留 Stage 1–6 和 scVI artifacts 不被覆盖。

## Stage 16：LAM candidate identity gate 重构

Stage 16 独立于 Step 7–13，且本轮不重跑 Step 7–13、不重训 scVI、不重新聚类。脚本 `16_rebuild_lam_identity_gate.py` 从 prepared AnnData 中所有 `condition=LAM` 细胞开始，使用连续 module/signature score 重建 candidate identity evidence：PMEL/MLANA/MITF 与 LAMCORE/CORE2/CORE3 为 identity anchors，ACTA2/ESR1/VEGFD/CTSK 为 supportive evidence，ciliated、AT2、myeloid/macrophage、endothelial/lymphatic endothelial、fibroblast、mesothelial、pericyte/VSMC 为 competing-lineage evidence。`FIGF` 全部 canonicalize 为 `VEGFD`；不再用任意两个 marker 的 `expression > 0` 规则作为主要判据。

阈值只使用独立 upstream LAM positive reference、normal/control 与明确 competing-lineage reference 校准，并分别执行四个 dataset 的 leave-one-dataset-out 验证。旧 consensus states 只在规则冻结后用于 post-hoc 诊断，不参与阈值选择。pericyte/VSMC 只在 identity anchors 弱时产生 conditional penalty，不能仅凭 ACTA2/MYH11 排除 LAM。Stage 16 输出到 `results/stage16/`，包括 `cell_identity_evidence.csv`、`identity_score_by_dataset.csv`、`reference_calibration.csv`、`leave_one_dataset_out_validation.csv`、`new_candidate_assignment.csv`、`new_candidate_by_old_state.csv` 与 `identity_gate_report.md`。可选 formal 777-gene LAMCORE reference 缺失时记录为 unavailable，不伪造输入，也不阻断核心审计。
## Stage 17：跨数据集 identity calibration 审计

Stage 17 只读取 Stage 16 的 `cell_identity_evidence.csv`、`reference_calibration.csv`、`leave_one_dataset_out_validation.csv`，并从 prepared AnnData 的 `layers["counts"]` 为 Stage 16 positive-reference 细胞补充 raw-count marker detection。重点审计 GSE190260 的 upstream CORE3-positive 细胞为何未被 core gate 恢复；不修改 Stage 16 candidate assignment、不重训 scVI、不重新聚类，也不使用旧 consensus state 调参。

审计内容包括：PMEL/MLANA/MITF、ACTA2/ESR1/VEGFD/CTSK、CORE2/CORE3、LAMCORE 777-gene 补充分数、competing-lineage penalty 分解、测序深度与 marker dropout、raw/within-dataset-z/within-dataset-percentile 三种反事实 score、四个 dataset 的 LODO recovery，以及 dataset-level root-cause attribution。新放入 `data-temp/LAM_core_signature_genes.csv` 的 formal reference 在 Stage 17 中按清理后的 `Gene` 列读取为 777 个基因，但由于它晚于 Stage 16 artifact 到达，只作为 supplemental audit score，不回填历史 Stage 16 score。

输出目录为 `results/stage17/`，包含 `positive_reference_failures.csv`、`component_scores_by_dataset.csv`、`marker_detection_by_dataset.csv`、`competing_lineage_penalty_audit.csv`、`counterfactual_calibration.csv`、`lodo_recovery.csv`、`root_cause_by_dataset.csv` 和 `identity_calibration_audit.md`。Stage 17 不据此自动冻结新 gate；GSE190260 的跨数据集阈值问题需在单独变更中决定是否重跑 Stage 16。

## Stage 18：State 15 LAM-core reference anchor 验证

Stage 18 将当前 consensus 中的 State 15 作为固定验证对象。其细胞 ID 集合必须恰为 200 个，并在 `state15_anchor_summary.json` 中保存排序后 ID 的 SHA-256；本阶段不重跑任何聚类、不重训 scVI、不修改 Stage 16 candidate gate，也不覆盖 Stage 1–17 的 artifact。

验证使用当前可用的 `LAM_core_signature_genes.csv` 777-gene formal LAMCORE，分别汇总 State 15、其余 19 个 consensus state、boundary 和 normal/control 的连续 score。另行计算 State 15 的 author-style enrichment（overall、dataset、patient 分层及 Fisher exact）、PMEL/MLANA/MITF、ACTA2/ACTG2/MYH11、VEGFD/CTSK、EMX2/HOXA11、ESR1、CORE1/2/3、mTOR、hormone、ECM/protease、HOX/PBX 及 competing-lineage profile，并与 State 18（pericyte/VSMC）、State 20（fibroblast）、State 12（endothelial）、State 7（AT2）、State 5（macrophage）逐一比较。

患者层面使用 `patient × State15` 与 `patient × comparator` 的 raw-count pseudobulk；latent 分析只读取既有 `state_model_scvi.h5ad` 的 `X_scVI`，记录 State 15 周围的 state/boundary/normal 邻域和距离—identity 梯度。normal 与 boundary 只用于辅助验证，不参与 State 15 的定义或 state 数量。Stage 18 输出到 `results/stage18/`，主要文件为 `state15_marker_profile.csv`、`state15_lamcore_summary.csv`、`state15_author_enrichment.csv`、`state15_vs_comparators.csv`、`state15_patient_pseudobulk.csv`、`state15_patient_consistency.csv`、`state15_latent_neighbors.csv`、`state15_latent_distance_gradient.csv` 和 `state15_anchor_report.md`。最终是否升级为正式 reference anchor 只能依据综合的 formal LAMCORE、独立 author evidence、跨患者/数据集复现、normal-lineage 区分和 latent continuity 判断；本阶段不设置单项自动升级规则。

## Stage 19：State 15 跨患者 identity calibration 审计

Stage 19 固定使用当前 5,378 个 high-confidence consensus candidate，先计算患者在 candidate pool 与 State 15 中的组成比例及 enrichment，避免把单个患者的细胞数直接等同于患者富集。分析对象仍是已冻结的 State 15，不重新聚类、不重训 scVI、不修改 gate。

author-style evidence 按上游文件的真实可用性解释：GSE135851 有真实逐细胞阳性标签，因此在该数据集内计算 State 15 enrichment 和 Fisher exact；GSE190260、GSE217108、GSE302356 的同名字段来自外部转换时的全 False 初始化，标记为 `not_assayed`，不能解释为 author-negative。

对每个含 State 15 的患者输出同一套 LAMCORE、CORE1/2/3、PMEL/MLANA/MITF、VEGFD/CTSK、ACTA2、ESR1、HOX/PBX、LAM myogenic、ECM 及 competing-lineage profile；再做 patient-matched `State15 vs 该患者其他 candidate`、7 个患者逐一 LOPO 的 pseudobulk/cosine/Pearson/latent-centroid 比较，以及删除 LAM1163 后的 State 15 与 State 18/20/12/7/5 对照。LOPO 与敏感性分析只比较现有标签和 profile，不产生新 candidate 或新 cluster。

输出目录为 `results/stage19/`，包括 `state15_patient_composition.csv`、`author_annotation_availability.csv`、`state15_author_enrichment_assayed.csv`、`state15_patient_profiles.csv`、`state15_patient_pseudobulk_profiles.csv`、`state15_patient_matched_comparison.csv`、`state15_lopo_validation.csv`、`state15_without_LAM1163.csv`、`stage19_manifest.json` 和 `state19_cross_patient_audit.md`。本次结果归为“患者富集但去除最大患者后生物学 profile 仍保留”的中间情形：LAM1163 的 State 15 composition enrichment 为 6.8301 倍，但其移除后仍有 73 个 State 15 细胞，且 LAMCORE 与患者匹配/LOPO 证据保持；这支持 State 15 是真实但患者丰度异质的候选状态，而不是仅由 LAM1163 产生的 cluster。该结论仍不替代正式 reference-anchor 升级所需的独立跨数据集 author evidence。

## Stage 21：State 15-centered manifold 独立验证

Stage 21 将 State 15 的 200 个细胞冻结为 reference anchor，验证对象严格改为其余 22,061 个 Stage 20 cells；复用 Stage 20 已计算的 `distance_to_state15` 和 `state_model_scvi.h5ad` 中的 `X_scVI`，不重训、不重聚类、不修改 candidate gate。对 777 个 formal LAMCORE genes 审计其与 scVI 4000 HVG、旧 gate marker 和当前表达矩阵的重叠，并构造 `LAMCORE_full`、`LAMCORE_no_gate`、`LAMCORE_outside_scVI`、`LAMCORE_independent` 四套 score。主要 gradient 结果同时提供全非-State15与非-State15 candidate-only 两个 scope，回归使用 dataset-standardized distance，患者作为 fixed effect；另以 20 个 quantile smooth bins 展示连续形状。

Stage 21 对四个 dataset 和 12 个 patient 重复 gradient，并显式保留 Stage 20 的非负患者、`n_State15`、distance range 和 gradient class。composition-matched null 固定进行 500 次：从非-State15 candidate pool 按 State15 的 patient×dataset 组成抽取 200 个假 anchor，比较真实 candidate-only independent slope 与 null slope distribution。State 16 仅作近/中/远 profile，boundary 仅作 evidence ranking，k=30 connectivity 用与 Stage 20 相同的 `X_scVI` scope 汇总，不产生新的 state/candidate。

实际结果：777 个 formal genes 中 220 个属于 scVI 4000 HVG，7 个属于旧 gate marker，554 个同时不属于二者，729 个在当前表达矩阵可用。真实 candidate-only `LAMCORE_independent` patient-adjusted slope 为 -0.01547，500 次 matched-anchor null 的经验双侧 p=0.001996；但 pooled non-State15 rank rho 为 +0.0757，最近与最远距离箱的 independent LAMCORE median 分别为 0.1083 和 0.1443，且患者方向异质。因此 Stage 21 checkpoint 为 `state15_lam_rich_gradient_but_not_robust_manifold`：独立 LAM evidence 在 State 15 外仍存在，且不是任意组成匹配假 anchor 都能产生的结果，但目前不能称为稳健统一 manifold。

## Stage 22：State 15 局部分支分解

Stage 22 专门分析 State 15 周围的局部方向，不再讨论全局 manifold 是否成立。输入固定继承 Stage 20/21：22,261 个 cells、State15 的 200 个 anchor、原始 20 维 `X_scVI`、Stage 21 的四套 LAMCORE score、CORE1/2/3 和 competing-lineage score。使用同一 candidate+boundary scope 的 k=30 kNN，取无向化图上的最短 1–3 hop；不重训、不重聚类、不修改 candidate gate 或任何 State 1–20 标签。

外部既有 state 只有在 1-hop 至少 10 个细胞且至少覆盖 2 个患者时自动进入 branch candidates，boundary 单独保留为局部方向投射。每条 branch 做 latent distance 与 graph-hop 的 near/mid/far profile，计算 independent LAMCORE、CORE1/2/3、LAM phenotype 和 competing lineage 的 Spearman/patient-adjusted slope，并对每个有足够细胞的患者重复 near/far 方向。每条 branch 另外从 1–3 hop 的非 branch 区域按 patient×dataset 和细胞数匹配抽取 500 次 null，报告真实 branch slope、null median 和经验 p。

实际结果选出 4 条主要方向：State16（1/2/3-hop=96/234/65，10 patients，4 datasets）、State12（22/165/327）、State20（22/320/384）和 State7（11/87/308）；State18 未达到直接连接筛选条件。State16 的 independent LAMCORE near/mid/far 为 0.1820/0.1128/0.1146，7 个患者有 ≥10 个 State16 cells，其中 5 个呈近端到远端下降；其 branch matched-null 经验双侧 p=0.02595，并被标记为 `LAM_to_lineage_transition_candidate`。State12、State20、State7 的 null p 分别为 0.67465、0.31537、1.0，当前标记为 `ordinary_lineage_adjacency`。Boundary 投射中 State16/12/20/7 分别获得 762/1700/841/2935 个方向 assignment，10,645 个保持 unresolved，不产生新的 LAM 标签。Stage 22 checkpoint 为 `local_branched_lam_manifold_candidate`：当前最值得继续验证的是 State15→State16 的局部分支，其他方向暂不称为 LAM branch。

## Stage 20：State 15-centered transcriptional manifold

Stage 20 冻结现有 consensus State 15 的 200 个细胞，不重新聚类、不重训 scVI、不修改 candidate gate。使用 `state_model_scvi.h5ad` 中现有的 `X_scVI`，主分析对象为 5,378 个 high-confidence candidates 加 16,883 个继承的 boundary cells；normal reference 只作为最后的远端对照，不参与 State 15 邻域图。`stage20_manifest.json` 记录冻结细胞 ID 的 SHA-256、latent artifact、`X_scVI` 维度以及 cohort scope。

距离轴完全由 latent geometry 构成：每个细胞计算到 State 15 reference 的 nearest、mean-5、mean-15 和 centroid distance，并在 candidate+boundary 上建立 k=30 邻居图，计算 State 15 neighbor fraction。再按 0–10%、10–20%、20–40%、40–60%、60–80%、80–100% 距离分箱，描述当前 State 1–20 和 boundary 组成；这些标签不用于重新定义任何状态。

在几何轴冻结后映射 LAMCORE/CORE1/2/3、melanocytic、LAM-support、myogenic、HOX/PBX、hormone、ECM/protease/mTOR 和 competing-lineage scores；单独审计 State 16 的 LAM/immune 共表达与技术复杂度，执行 patient- 和 dataset-level distance-gradient consistency，并将 boundary 投射到最近 candidate state。normal 只输出 `normal_remote_summary.csv`。Stage 20 的初次 checkpoint 为：pooled LAMCORE 从最近距离分箱的 0.2152 降至最远分箱的 0.1556，4/4 dataset 的 distance~LAMCORE rho 为负，12 个 patient 中 8 个可估计 rho 为负，因此支持一个以 State 15 为中心的跨数据集转录梯度，但仍需结合后续 patient-aware/identity 验证解释其是否为统一 manifold。

输出目录为 `results/stage20/`，包括 `state15_cell_distances.csv`、`state15_distance_bins.csv`、`state15_identity_gradient.csv`、`state15_lineage_gradient.csv`、`state16_cell_audit.csv`、`state16_lam_immune_coexpression.csv`、`state16_doublet_audit.csv`、`patient_gradient_consistency.csv`、`dataset_gradient_consistency.csv`、`boundary_state15_projection.csv`、`state15_centered_manifold.csv`、`stage20_manifest.json` 和 `stage20_manifold_report.md`。运行脚本为 `20_state15_centered_manifold.py`。

## Stage 23：State 15 latent-space visualization

Stage 23 是展示阶段，不新增生物学判定。它固定读取 Stage 20 的 22,261-cell candidate+boundary 主表、Stage 21 的四套 LAMCORE/independent scores、Stage 22 的 k=30 局部图/branch/null 结果以及 `state_model_scvi.h5ad` 的 20-dimensional `X_scVI`。State 15 的 200 个细胞只作为冻结的 reference anchor；不重训 scVI、不重新 Leiden/consensus、不修改 candidate gate 和 State 1–20 标签。Stage 21 没有把 anchor 写入 validation score table，因此脚本在内存中用完全相同的 Stage 21 score modules 补齐 anchor；已有 Stage 20 program scores作为缺失 fallback，但不写回任何 artifact。

静态图包括：全局 latent UMAP（state 与 independent LAMCORE 两 panel）、State15 的 1–3 hop kNN 局部图（state 与 independent LAMCORE 两张图）、排除 State15 的 distance×independent LAMCORE pooled/patient facets、State16 near/mid/far program heatmap，以及 State16/12/20/7 的 branch matched-null slope 对照。交互 HTML 包括从 `X_scVI` 单进程计算的 3D UMAP、局部 3D kNN graph 和 3D PCA 正交模型；hover 统一提供 cell、patient、dataset、state、branch、hop、distance、independent LAMCORE、CORE1/CORE3、T/NK、VSMC/pericyte 和 candidate/boundary 信息，dropdown 支持 State、LAMCORE、patient、dataset、candidate/boundary 和 branch 着色。

2D 直接复用 `state_model_scvi.h5ad` 的 `obsm['X_umap']`（缺失时才从 `X_scVI` 计算），3D UMAP 使用 `umap-learn`、`n_neighbors=30`、单进程；PCA 仅作对照。为控制交互文件体积，局部 kNN 边采用同一 X_scVI 几何确定性重建，默认最多绘制 100,000 条边。输出目录为 `results/stage23_visualization/`，并由 `visualization_manifest.json` 记录输入来源、冻结 anchor 哈希、embedding 来源、局部节点/边数和“无模型改变”约束。

## Stage 19：State 15 跨患者 identity calibration 审计

Stage 19 固定使用当前 5,378 个 high-confidence consensus candidate，先计算患者在 candidate pool 与 State 15 中的组成比例及 enrichment，避免把单个患者的细胞数直接等同于患者富集。分析对象仍是已冻结的 State 15，不重新聚类、不重训 scVI、不修改 gate。

author-style evidence 按上游文件的真实可用性解释：GSE135851 有真实逐细胞阳性标签，因此在该数据集内计算 State 15 enrichment 和 Fisher exact；GSE190260、GSE217108、GSE302356 的同名字段来自外部转换时的全 False 初始化，标记为 `not_assayed`，不能解释为 author-negative。

对每个含 State 15 的患者输出同一套 LAMCORE、CORE1/2/3、PMEL/MLANA/MITF、VEGFD/CTSK、ACTA2、ESR1、HOX/PBX、LAM myogenic、ECM 及 competing-lineage profile；再做 patient-matched `State15 vs 该患者其他 candidate`、7 个患者逐一 LOPO 的 pseudobulk/cosine/Pearson/latent-centroid 比较，以及删除 LAM1163 后的 State 15 与 State 18/20/12/7/5 对照。LOPO 与敏感性分析只比较现有标签和 profile，不产生新 candidate 或新 cluster。

输出目录为 `results/stage19/`，包括 `state15_patient_composition.csv`、`author_annotation_availability.csv`、`state15_author_enrichment_assayed.csv`、`state15_patient_profiles.csv`、`state15_patient_pseudobulk_profiles.csv`、`state15_patient_matched_comparison.csv`、`state15_lopo_validation.csv`、`state15_without_LAM1163.csv`、`stage19_manifest.json` 和 `state19_cross_patient_audit.md`。本次结果归为“患者富集但去除最大患者后生物学 profile 仍保留”的中间情形：LAM1163 的 State 15 composition enrichment 为 6.8301 倍，但其移除后仍有 73 个 State 15 细胞，且 LAMCORE 与患者匹配/LOPO 证据保持；这支持 State 15 是真实但患者丰度异质的候选状态，而不是仅由 LAM1163 产生的 cluster。该结论仍不替代正式 reference-anchor 升级所需的独立跨数据集 author evidence。
