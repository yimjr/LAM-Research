# Stage 1–23 summary

本文件是因果链摘要；完整证据路径见 `artifact_index.csv`，历史冲突见 `narrative_audit.csv`。

## Stage 1 — 输入盘点与上游继承清单

**问题**：四套 AnnData、四套 candidate/program 结果和共享配置是否可追溯、可用？
**输入**：LAM-Cell-Research / data-temp / data/upstream；config 与 manifests
**方法/参数**：只读检查路径、字段、counts layer、三个 pool 层级，并递归发现同类 upstream 结果。
**输出**：results/stage1_6/input_inventory.{csv,json}; upstream_annotation_manifest.csv
**结果/checkpoint**：四套 AnnData 和四套 candidate_pool_labels 均 ready；normal 为可选输入。
**为什么进入下一步**：确定直接继承转换后的 AnnData 和 upstream annotations，不重做 GEO 转换。
**后续修正/边界**：无；后续阶段继续使用这份继承边界。

## Stage 2 — 继承与准备

**问题**：能否把 upstream annotation、registry 和 raw counts 安全映射到统一分析表？
**输入**：四套上游 AnnData、candidate_pool_labels、core3/program 结果、donor_registry.yaml
**方法/参数**：按 cell ID/sample key 合并，canonicalize FIGF→VEGFD，固定 high-confidence candidate 和 boundary，验证 counts。
**输出**：data/interim/inherited；results/stage1_6/upstream_inheritance.json
**结果/checkpoint**：lam_candidate 只等于 pool_high_confidence；unrestricted 仅审计。
**为什么进入下一步**：为 QC/harmonization 建立统一 prepared 输入，同时把旧标签留在 upstream 命名空间。
**后续修正/边界**：Stage 16 的新 gate 只作诊断，不替换这里冻结的主线 candidate pool。

## Stage 3 — QC 与 gene harmonization

**问题**：在不重复核心 QC 的前提下，能否为外部数据补齐本项目 QC 并建立共同 gene universe？
**输入**：继承后的四套 AnnData；core qc_pass；外部 raw counts
**方法/参数**：GSE135851 继承原 qc_pass；外部使用 min_genes=200、min_counts=500、scRNA mt<20%、snRNA mt<10%；保留 doublet。
**输出**：state_model_prepared.h5ad；qc_summary.csv
**结果/checkpoint**：外部 QC 后：GSE190260 38,701、GSE217108 6,941、GSE302356 21,771 cells pass。
**为什么进入下一步**：后续 NMF/scVI 共享 alias-corrected gene universe；GSE135851 不被二次过滤。
**后续修正/边界**：Stage 15–17 发现的是 upstream candidate gate 问题，不是本阶段 QC 合并错误。

## Stage 4 — PCA/NMF baseline

**问题**：在 State Modeling 自己的输入矩阵上，能否建立不依赖旧 state 的 baseline？
**输入**：state_model_prepared.h5ad；counts layer
**方法/参数**：counts→library-size normalization→log1p→HVG(4000)→PCA/UMAP/Leiden/NMF；NMF 5 components、top 2000 features、max_iter 400、最多 12000 cells。
**输出**：state_model_baseline.h5ad；baseline_cluster_summary.csv；nmf_cell_scores.csv；nmf_top_genes.csv
**结果/checkpoint**：获得独立 baseline；旧 program/state 只保留作 post-hoc 对照。
**为什么进入下一步**：明确 NMF 不使用 scaled/PCA/scVI latent，scVI 也不使用 NMF 矩阵。
**后续修正/边界**：没有用 baseline cluster 数直接定义最终 20 states。

## Stage 5 — scVI latent model

**问题**：只用 dataset 做 batch correction 时，是否能构建可复用 latent space？
**输入**：state_model_prepared.h5ad 的 layers[counts]
**方法/参数**：scVI layer=counts、batch_key=dataset、无 categorical covariate、n_latent=20、n_layers=2、n_hidden=128、max_epochs=200、early stopping；CUDA 优先。
**输出**：state_model_scvi.h5ad；scvi_model/model.pt；scvi_training_manifest.json
**结果/checkpoint**：55,238 cells 的 200-epoch 模型完成；assay 留在 metadata，不进入 covariate。
**为什么进入下一步**：Stage 6–23 均复用同一 X_scVI，不再重训。
**后续修正/边界**：无；后续所有结构分析都明确排除重新训练。

## Stage 6 — LAM-only latent structure checkpoint

**问题**：high-confidence LAM candidate 内是否存在跨患者共享、且不由 patient/dataset/assay 主导的结构？
**输入**：X_scVI；5,378 pool_high_confidence cells；boundary/normal 仅辅助
**方法/参数**：LAM-only kNN→Leiden；独立网格 n_neighbors 15/30/50×resolution 0.2/0.4/0.6；保留全体 33-cluster 标签作审计，不用于 LAM 内部判定。
**输出**：stage6_parameter_grid.csv；stage6_checkpoint.json；grid/cluster/ARI 表
**结果/checkpoint**：GO；12 个 LAM-only reference clusters、12 个跨至少 2 patients 的 qualified clusters；patient ARI 0.152554、dataset ARI 0.086786、assay ARI 0.009094。
**为什么进入下一步**：允许进入 consensus；也暴露出单一参数和固定 preprocess 参数不应混用。
**后续修正/边界**：Stage 7 将 9 个 grid configuration 等权并形成 20 个 consensus states；12 与 20 是不同对象，不是覆盖。

## Stage 7 — Consensus stability

**问题**：哪些 LAM-only 结构跨 grid/seed 稳定？
**输入**：同一 X_scVI 的 5,378 candidate cells；9 configurations、21 raw partitions
**方法/参数**：每个 configuration 内先平均 seed co-assignment，再让 9 configurations 等权；最终生成 5,378×5,378 float32 co-assignment 并做完整距离 average-linkage。
**输出**：state_consensus_assignments.csv；state_stability_summary.csv；cluster_matching_across_grid.csv；coassignment_matrix.npz；consensus_dendrogram.npz
**结果/checkpoint**：形成 20 个当前 consensus states；seed stability 单独记录。
**为什么进入下一步**：把 Stage 6 的单一参考 partition 转为可追踪的 20-state consensus。
**后续修正/边界**：20-state labels 在 Stage 8–23 冻结；后续 hierarchy 只作解释。

## Stage 8 — Leave-one-out robustness

**问题**：去掉一个 patient/dataset 后，consensus state 是否仍可恢复？
**输入**：Stage 7 consensus；X_scVI；5,378 candidate cells
**方法/参数**：full-data 30/0.4 reference；LOO 前将 reference、consensus 和 LOO cluster 都限制在 retained cells；允许 split/merge matching。
**输出**：loo_runs.csv；loo_cluster_matches.csv；loo_state_summary.csv；full_reference_consensus_matches.csv
**结果/checkpoint**：产生连续的 baseline Jaccard、LOO recovery 和 additional loss。
**为什么进入下一步**：为每个 state 的结构稳定性提供 patient/dataset 支持，而不是只看 full partition。
**后续修正/边界**：Stage 11 将这些连续指标与患者级生物学复现合并。

## Stage 9 — Hierarchy and continuum

**问题**：20 states 是平级离散状态，还是存在 parent/substate/连续连接？
**输入**：consensus labels；X_scVI；boundary transitions
**方法/参数**：state latent distance、connectivity/PAGA 等价汇总、split/merge tree、boundary transitions；不重新聚类。
**输出**：state_distance_matrix.csv；state_connectivity.csv；state_split_merge_tree.csv；boundary_state_transitions.csv
**结果/checkpoint**：为 atlas 提供描述性 parent/substate 结构。
**为什么进入下一步**：避免把 resolution 改变产生的层级现象硬写成 20 个互不相干状态。
**后续修正/边界**：Stage 22 进一步把全局问题收窄到 State15 周围局部分支。

## Stage 10 — Per-state biology and DE

**问题**：每个 consensus state 的表达/程序差异是否有患者级证据？
**输入**：state_model_prepared.h5ad full gene universe/raw counts；consensus state
**方法/参数**：每个 State_k 单独对同患者 Rest_of_LAM 做 patient×group pseudobulk；设计 ~ patient_id + group；至少 3 patients 才正式 DE；每 state 独立 FDR。
**输出**：state_de_results.csv（159,059 rows）；state_markers.csv；state_pseudobulk_counts.csv；state_program_scores.csv；pathway/regulon placeholders
**结果/checkpoint**：20 个独立 state-vs-rest 分析路径，不存在 all-state 多分类模型。
**为什么进入下一步**：提供 state-level markers/program evidence，同时承认 pathway/regulon 当前 unavailable。
**后续修正/边界**：Stage 24 不重新做 DE，只整理其支持范围。

## Stage 11 — Patient-level reproducibility

**问题**：细胞层面的 state 是否能转化为跨患者 evidence？
**输入**：Stage 8 LOO；Stage 10 DE/pseudobulk；patient metadata
**方法/参数**：patient×state matrix；cells/fraction/signature/log2FC/direction/dataset coverage；保留 structural stability 与 biological reproducibility 两维。
**输出**：patient_state_matrix.csv；state_reproducibility_summary.csv
**结果/checkpoint**：State15 structural 0.854111、biological 0.374758；State16 0.719633/0.323088；State18 0.940033/0.408960；State20 0.912545/0.420718。
**为什么进入下一步**：把“cluster 很稳定”与“跨患者 biology 可复现”分开。
**后续修正/边界**：Stage 19 对 State15 的患者组成和最大患者敏感性作专门审计。

## Stage 12 — Boundary and normal auxiliary validation

**问题**：候选边界和正常参考能否帮助解释 state 的外部邻域？
**输入**：X_scVI；boundary 16,883；normal 32,977
**方法/参数**：boundary/normal 邻域与距离；不加入 state 数量，不改变 candidate/state 定义。
**输出**：boundary_validation.csv；normal_validation.csv；state_auxiliary_summary.csv
**结果/checkpoint**：normal 与 boundary 仅为辅助；State15 normal mean distance 3.981922，State18 3.023263，State20 3.240648。
**为什么进入下一步**：为后续 State15 anchor 与 manifold 分析提供参照。
**后续修正/边界**：Stage 20–22 保持 normal remote、boundary projection 的辅助地位。

## Stage 13 — State atlas and hypotheses

**问题**：如何用连续证据汇总 20 states，而不把单一分数当 confidence？
**输入**：Stages 7–12 的 state、DE、program、LOO、auxiliary 结果
**方法/参数**：汇总 structural/biological/coverage/normal/boundary/upstream correspondence；不设硬性 high/medium/low 门槛。
**输出**：state_atlas.csv/json；state_hypothesis_candidates.csv；state_atlas.h5ad
**结果/checkpoint**：形成第一版 20-state atlas 与 10 个 hypothesis candidates。
**为什么进入下一步**：为 Stage15–22 提供固定状态编号和证据背景。
**后续修正/边界**：Stage 24 重新解释这些 labels 的生物学含义，但不改写 atlas artifact。

## Stage 14 — Consensus/upstream merge

**问题**：如何把 frozen consensus 与 upstream annotation 放在同一逐细胞表？
**输入**：Stage 7 consensus；upstream candidate/state/program annotations
**方法/参数**：按 cell ID 一对一合并，保留 upstream 命名空间；不让旧 labels 进入 scVI/NMF。
**输出**：state_consensus_with_upstream_annotations.csv；merge manifest
**结果/checkpoint**：5,378 cells 的合并字段可追溯；Stage 15 audit later confirmed merge inconsistencies=0。
**为什么进入下一步**：支持 Stage 15 candidate identity audit 和后续 post-hoc comparison。
**后续修正/边界**：外部 candidate/state 结果的完整继承在 Stage 1–6 manifest 中保留。

## Stage 15 — Candidate identity audit

**问题**：5,378 candidate 是否主要由过宽的 marker-combo gate 产生？
**输入**：原 candidate_pool_labels；merged consensus；原始 marker_expr/counts
**方法/参数**：ID/字段一对一核对；重算原规则；FIGF/VEGFD duplicate audit；组合、UMI 和 state-level diagnostics；只读。
**输出**：annotation_merge_audit.csv；rule_recalculation_audit.csv；marker_patterns.csv；root_cause_evidence.csv
**结果/checkpoint**：A merge error=0；FIGF/VEGFD duplicate-pass=0；5,238 marker-combo、140 author/formal；1,443 个 marker-combo candidate 仅有 1-UMI 支持。
**为什么进入下一步**：根因定位到 C：任意两个 marker>0 gate 特异性不足，而非 Stage Modeling 搬错。
**后续修正/边界**：Stage 16 提出连续 identity gate，但没有写回主线 candidate pool。

## Stage 16 — Identity gate reconstruction

**问题**：能否用 identity anchors+support+competing lineage 重建更可解释的 gate？
**输入**：所有 condition=LAM cells；连续 module scores；独立正/负参照
**方法/参数**：PMEL/MLANA/MITF+LAMCORE/CORE evidence 为 identity；ACTA2/ESR1/VEGFD/CTSK 为 support；竞争谱系为 penalty；不使用旧 state 调参。
**输出**：cell_identity_evidence.csv；reference_calibration.csv；LODO；new_candidate_assignment.csv
**结果/checkpoint**：LAM_core 208、boundary 65,930、non_LAM_like 24,503；这是独立诊断 gate，不替换 frozen candidate。
**为什么进入下一步**：使 Stage 17 能够专门追踪外部 positive reference 的漏检。
**后续修正/边界**：Stage 17 发现 GSE190260 的 dataset calibration/dropout/penalty 问题；Stage 24 不继续调 gate。

## Stage 17 — Cross-dataset identity calibration audit

**问题**：GSE190260 为什么漏掉 upstream CORE3-positive？
**输入**：Stage 16 artifacts；formal 777-gene signature；raw counts
**方法/参数**：正参考逐细胞 score decomposition、dropout/depth、competing penalty、raw/z/percentile counterfactual、LODO；不生成新 assignment。
**输出**：positive_reference_failures.csv；component_scores_by_dataset.csv；counterfactual_calibration.csv；root_cause_by_dataset.csv
**结果/checkpoint**：GSE190260 2,117 positives：core recovery 0、core+boundary 0.708077；median final score -1.006577、score shift -6.610707、median penalty 4.950320；primary category competing_penalty_only，且 identity/support dropout 均高。
**为什么进入下一步**：证明跨数据集 score scale/marker dropout/penalty 联合作用，不能把外部 0 author-style 当 negative。
**后续修正/边界**：Stage 18 formal LAMCORE anchor validation独立使用 777 genes；不改 Stage16。

## Stage 18 — State15 anchor validation

**问题**：冻结的 State15 能否作为 LAM-core reference anchor？
**输入**：State15 200 cells；777-gene formal LAMCORE；author labels；normal/boundary/comparators
**方法/参数**：LAMCORE score、author enrichment/Fisher、marker/program、patient pseudobulk、comparator 和 latent-neighborhood validation。
**输出**：state15_anchor_summary.json；state15_lamcore_summary.csv；state15_vs_comparators.csv；anchor report
**结果/checkpoint**：State15 LAMCORE median 0.5125；normal 0.0513；overall author enrichment 13.9408、Fisher p=2.229e-36；decision provisional_reference_candidate_not_formally_upgraded。
**为什么进入下一步**：State15 被冻结为后续 reference anchor candidate，而不是正式 classifier。
**后续修正/边界**：Stage19 量化 LAM1163 enrichment；最终仍保持 provisional。

## Stage 19 — State15 cross-patient audit

**问题**：State15 是否只是 LAM1163 这个大患者造成的？
**输入**：frozen State15；candidate pool composition；author availability；patient pseudobulk/LOPO
**方法/参数**：patient composition baseline、patient-matched comparison、7-patient LOPO、remove-LAM1163 sensitivity；author unavailable 标成 not_assayed。
**输出**：state15_patient_composition.csv；author_annotation_availability.csv；LOPO；state15_without_LAM1163.csv
**结果/checkpoint**：LAM1163 占 candidate 9.2971%、State15 63.5%，enrichment 6.8301；去除后保留 73 cells、LAMCORE median 0.4683；author labels 仅 GSE135851 可检验。
**为什么进入下一步**：把 State15 定位为患者富集但 profile 可保留的中间结果。
**后续修正/边界**：为 Stage20 的 anchor-centered geometry 提供组成警示。

## Stage 20 — State15-centered latent geometry

**问题**：以 State15 为中心，周围细胞的 identity/program 是否沿 latent distance 变化？
**输入**：X_scVI；State15 200；candidate 5,378+boundary 16,883；normal remote
**方法/参数**：latent-only distance、k=30 neighbors、distance bins、identity/lineage mapping、State16 co-expression、patient/dataset gradients。
**输出**：state15_centered_manifold.csv；distance bins/gradients；State16 audits；stage20 report
**结果/checkpoint**：初始 checkpoint supports_lam_centered_transcriptional_manifold；full LAMCORE 近端→远端 0.2152→0.1556，4/4 dataset rho<0，8/12 patient rho<0。
**为什么进入下一步**：提出 State15-centered manifold 假设，并把验证对象从单个 state 扩展到邻域。
**后续修正/边界**：Stage21 明确削弱 pooled/global manifold；Stage22 改为局部分支 candidate。

## Stage 21 — Independent State15 manifold validation

**问题**：去掉 anchor、自身 gate/scVI feature overlap 和 composition 后，gradient 是否仍成立？
**输入**：Stage20 distance；非-State15 22,061；777-gene LAMCORE；candidate-only null pool
**方法/参数**：LAMCORE full/no_gate/outside_scVI/independent；patient-adjusted regression、dataset/patient replication、500 matched fake anchors、State16/boundary auxiliary。
**输出**：gradient_models.csv；dataset/patient gradients；matched_anchor_null.csv；stage21 report
**结果/checkpoint**：554 genes independent；candidate-only independent slope -0.015466，null empirical two-sided p=0.001996；但 pooled rank rho +0.075665、near/far medians 0.1083/0.1443，患者异质；checkpoint state15_lam_rich_gradient_but_not_robust_manifold。
**为什么进入下一步**：不再支持稳健统一 global manifold，但保留局部/方向性结构的可能。
**后续修正/边界**：Stage22 将问题改成 local branch decomposition。

## Stage 22 — State15 local branch decomposition

**问题**：State15 周围哪些局部方向真正保留 LAM identity？
**输入**：同一 X_scVI；22,261 candidate+boundary；Stage21 scores；frozen State labels
**方法/参数**：k=30 不重新 Leiden；只使用 State15 局部 1–3-hop；branch eligibility 按 1-hop 患者数；real/null 统一 local scope，并按 patient×dataset、细胞数和五档距离结构匹配；直接经验左右尾、BH-q、每患者 slope 与 LOPO。
**输出**：branch_candidates.csv；branch evidence/gradients；boundary assignments；stage22 report
**结果/checkpoint**：修正后仍入选 State_12, State_16, State_20, State_7；State_16: slope=-0.023931, left_p=0.439122, left_BH_q=1.000000, two_sided_p=0.878244, two_sided_BH_q=0.878244, label=ordinary_lineage_adjacency; State_12: slope=-0.022764, left_p=0.988024, left_BH_q=1.000000, two_sided_p=0.027944, two_sided_BH_q=0.055888, label=ordinary_lineage_adjacency; State_20: slope=0.010504, left_p=1.000000, left_BH_q=1.000000, two_sided_p=0.003992, two_sided_BH_q=0.015968, label=ordinary_lineage_adjacency; State_7: slope=0.001486, left_p=0.896208, left_BH_q=1.000000, two_sided_p=0.211577, two_sided_BH_q=0.282102, label=ordinary_lineage_adjacency；Stage22 checkpoint=ordinary_lineage_adjacency_dominates。
**为什么进入下一步**：依据校正后的 local geometry 重新评估分支；不再把 raw empirical p 当作唯一证据，State16 原 transition 标签被降级/撤回（若当前 q 不支持）。
**后续修正/边界**：Stage23 已按新 branch/boundary 输出重生成；Stage24 采用本次修正数字，未改变 State15、X_scVI 或 State1–20 标签。

## Stage 23 — Latent-space visualization

**问题**：如何把 Stage15–22 已有结构直观看清，而不新增结论？
**输入**：X_scVI；Stage20/21/22 frozen tables；State15 anchor
**方法/参数**：全局/局部 2D static plots；3D UMAP、3D local graph、3D PCA；统一 hover/dropdown；不重训、不重聚类。
**输出**：results/stage23_visualization/ 五组 PNG/PDF 与三个 HTML；visualization_manifest.json
**结果/checkpoint**：可视化资产完成；2D 优先复用已有 UMAP，3D UMAP/PCA 只用于展示/对照。
**为什么进入下一步**：为 Stage24 artifact provenance 和最终报告提供可引用图件。
**后续修正/边界**：无；Stage24 不把可视化解释升级为新证据。
