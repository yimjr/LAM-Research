# LAM-State-Modeling：Stage 24 最终报告原材料包

> 这是面向最终写作的证据材料汇总，不是对历史 artifact 的覆盖。Stage 24 只读取已有结果并生成本目录文件；不重训 scVI、不重新聚类、不修改 candidate gate、State15、State1–20 或 branch artifact。历史矛盾和结论变化按时间保留。

## 1. 研究目标与目标演变

最初问题是：能否在跨数据集、跨患者的 LAM high-confidence candidate 中，用深度 latent representation 找到可重复的 LAM state structure。阶段 1–6 先完成输入继承、外部 QC、NMF/PCA baseline、scVI 和 LAM-only GO/NO-GO。Stage6 的 GO 只说明 high-confidence candidate 内存在值得继续检查的 latent structure，并不证明这些 cell 都是 LAM。

Stage7–13 将目标变成：哪些结构跨参数、患者和数据集稳定，哪些 state 能得到 patient-aware biological support。随后 Stage15 暴露出原 candidate gate 的高召回/低特异性，Stage16–17 改为独立 identity audit，而不是把新 gate 偷换进主模型。Stage18–19 冻结 State15，验证它是否能作为 LAM-core reference anchor；结果支持“LAM-rich、患者丰度异质、去掉 LAM1163 后 profile 仍保留”的 provisional interpretation，但没有正式升级。

Stage20 进一步提出 State15-centered global manifold。Stage21 在排除 anchor、scVI HVG/gate overlap 和 composition matched null 后保留了一部分独立 gradient，但否定了“稳健统一 global manifold”的强表述。修正后的 Stage22 仍选择 State_12, State_16, State_20, State_7；其 State16 为 ordinary_lineage_adjacency （raw empirical p=0.878244，BH q=0.878244）。 State12、20、7 的标签和 State15 局部关系也以修正后的结果为准。Stage23 仅把这些结果可视化。因而 Stage24 的最终科学问题不是“20 个 cluster 是否都是 LAM”，而是：一个高置信度候选池中，State15 是否代表最强 LAM-rich core，以及其周围是否存在有限、局部且方向依赖的延伸。

## 2. 数据和分析范围

主要数据集为 GSE135851（core reproduction）、GSE190260、GSE217108、GSE302356（external converted AnnData）。Stage1 inventory 记录原始转换对象大小：分别为 30,708、39,979、12,396、23,759 cells；pool_high_confidence 分别为 535、1,564、1,075、2,681，合计 5,855 个 upstream high flags。经过继承/QC/共同 gene universe 后，主线 Stage6–23 使用 5,378 个 high-confidence candidate cells、12 patients、4 datasets。

Stage1–6 inventory 中的 pool counts 为：

| dataset | AnnData cells | genes | high | broad | unrestricted |
|---|---:|---:|---:|---:|---:|
| GSE135851 | 30,708 | 63,677 | 535 | 3,940 | 23,228 |
| GSE190260 | 39,979 | 33,694 | 1,564 | 7,811 | 39,979 |
| GSE217108 | 12,396 | 36,601 | 1,075 | 3,855 | 12,396 |
| GSE302356 | 23,759 | 38,224 | 2,681 | 8,727 | 23,759 |

外部补充 QC 后 pass cells 为 GSE190260 38,701、GSE217108 6,941、GSE302356 21,771；GSE135851 直接继承原 reproduction baseline `qc_pass`。Stage12 normal auxiliary cohort 为 32,977 cells；boundary 为 16,883 cells。State15 为冻结 200 cells、7 patients、4 datasets；State16 为 396 cells。normal/boundary 从未参与 State6–7 的核心 state 数量定义。

候选池语义固定为：`lam_candidate = pool_high_confidence`；`boundary = pool_broad_lam_like AND NOT pool_high_confidence`；`pool_unrestricted_lam` 只审计。由于 unrestricted 近似所有 condition=LAM cell，三层不能做并集。

## 3. 核心方法和关键参数

- 数据继承：优先读取 `../LAM-Cell-Research`，否则 data-temp/data/upstream；不重新 GEO 转换；patient/donor mapping 唯一来源为 `donor_registry.yaml`。
- Alias：导入后、gene 去重前执行 `FIGF → VEGFD`；Stage15 审计显示同一细胞 FIGF+VEGFD duplicate-pass=0。
- QC：core 继承原 `qc_pass`；外部 min_genes=200、min_counts=500、scRNA mt<20%、snRNA mt<10%；doublet 保留、不重跑 caller、不默认删除。
- NMF：counts→library-size normalization→log1p→HVG=4000→top 2000 features→5 components，max_iter=400，最多 12,000 cells。
- scVI：`layer='counts'`、`batch_key='dataset'`、categorical covariates=[]、20 latent、2 layers、128 hidden、200 epochs、early stopping；assay 仅 metadata。后续 Stage7–23 全部复用同一个 `X_scVI`。
- Stage6：LAM-only candidate graph 的参数独立于 preprocess；9 个 configuration 为 n_neighbors 15/30/50 × resolution 0.2/0.4/0.6。
- Stage7：21 raw partitions 先在 configuration 内平均 seed，再 9 configuration 等权；最终 co-assignment 允许一份 5,378×5,378 float32 matrix 和完整 average-linkage。
- Stage8：full-data reference 为 n_neighbors=30、resolution=0.4；LOO overlap 只在 retained cells 上计算，并区分 full-reference→consensus baseline 与 LOO additional loss。
- Stage10：每个 state 独立 `State_k vs Rest_of_LAM`；patient×group pseudobulk；设计 `~ patient_id + group`；正式 DE 至少 3 patients；不使用统一多分类模型。
- Stage18–22：State15 的 200 cells 固定为 anchor；Stage20 distance 完全由 X_scVI；Stage21 的 777 LAMCORE validation 拆为 full/no_gate/outside_scVI/independent；Stage21 matched null 500 次；Stage22 branch null 每条 500 次、局部 k=30、不重新 Leiden。Stage22 修正版用 1-hop 患者数筛选、real/null 相同的 1–3-hop local scope、patient×dataset+距离分箱匹配、直接经验尾部、BH-q 和患者级 LOPO。

## 4. Stage 1–23 完整过程

| stage | title | research_question | checkpoint | why_next | later_revision_or_scope |
| --- | --- | --- | --- | --- | --- |
| 1 | 输入盘点与上游继承清单 | 四套 AnnData、四套 candidate/program 结果和共享配置是否可追溯、可用？ | 四套 AnnData 和四套 candidate_pool_labels 均 ready；normal 为可选输入。 | 确定直接继承转换后的 AnnData 和 upstream annotations，不重做 GEO 转换。 | 无；后续阶段继续使用这份继承边界。 |
| 2 | 继承与准备 | 能否把 upstream annotation、registry 和 raw counts 安全映射到统一分析表？ | lam_candidate 只等于 pool_high_confidence；unrestricted 仅审计。 | 为 QC/harmonization 建立统一 prepared 输入，同时把旧标签留在 upstream 命名空间。 | Stage 16 的新 gate 只作诊断，不替换这里冻结的主线 candidate pool。 |
| 3 | QC 与 gene harmonization | 在不重复核心 QC 的前提下，能否为外部数据补齐本项目 QC 并建立共同 gene universe？ | 外部 QC 后：GSE190260 38,701、GSE217108 6,941、GSE302356 21,771 cells pass。 | 后续 NMF/scVI 共享 alias-corrected gene universe；GSE135851 不被二次过滤。 | Stage 15–17 发现的是 upstream candidate gate 问题，不是本阶段 QC 合并错误。 |
| 4 | PCA/NMF baseline | 在 State Modeling 自己的输入矩阵上，能否建立不依赖旧 state 的 baseline？ | 获得独立 baseline；旧 program/state 只保留作 post-hoc 对照。 | 明确 NMF 不使用 scaled/PCA/scVI latent，scVI 也不使用 NMF 矩阵。 | 没有用 baseline cluster 数直接定义最终 20 states。 |
| 5 | scVI latent model | 只用 dataset 做 batch correction 时，是否能构建可复用 latent space？ | 55,238 cells 的 200-epoch 模型完成；assay 留在 metadata，不进入 covariate。 | Stage 6–23 均复用同一 X_scVI，不再重训。 | 无；后续所有结构分析都明确排除重新训练。 |
| 6 | LAM-only latent structure checkpoint | high-confidence LAM candidate 内是否存在跨患者共享、且不由 patient/dataset/assay 主导的结构？ | GO；12 个 LAM-only reference clusters、12 个跨至少 2 patients 的 qualified clusters；patient ARI 0.152554、dataset ARI 0.086786、assay ARI 0.009094。 | 允许进入 consensus；也暴露出单一参数和固定 preprocess 参数不应混用。 | Stage 7 将 9 个 grid configuration 等权并形成 20 个 consensus states；12 与 20 是不同对象，不是覆盖。 |
| 7 | Consensus stability | 哪些 LAM-only 结构跨 grid/seed 稳定？ | 形成 20 个当前 consensus states；seed stability 单独记录。 | 把 Stage 6 的单一参考 partition 转为可追踪的 20-state consensus。 | 20-state labels 在 Stage 8–23 冻结；后续 hierarchy 只作解释。 |
| 8 | Leave-one-out robustness | 去掉一个 patient/dataset 后，consensus state 是否仍可恢复？ | 产生连续的 baseline Jaccard、LOO recovery 和 additional loss。 | 为每个 state 的结构稳定性提供 patient/dataset 支持，而不是只看 full partition。 | Stage 11 将这些连续指标与患者级生物学复现合并。 |
| 9 | Hierarchy and continuum | 20 states 是平级离散状态，还是存在 parent/substate/连续连接？ | 为 atlas 提供描述性 parent/substate 结构。 | 避免把 resolution 改变产生的层级现象硬写成 20 个互不相干状态。 | Stage 22 进一步把全局问题收窄到 State15 周围局部分支。 |
| 10 | Per-state biology and DE | 每个 consensus state 的表达/程序差异是否有患者级证据？ | 20 个独立 state-vs-rest 分析路径，不存在 all-state 多分类模型。 | 提供 state-level markers/program evidence，同时承认 pathway/regulon 当前 unavailable。 | Stage 24 不重新做 DE，只整理其支持范围。 |
| 11 | Patient-level reproducibility | 细胞层面的 state 是否能转化为跨患者 evidence？ | State15 structural 0.854111、biological 0.374758；State16 0.719633/0.323088；State18 0.940033/0.408960；State20 0.912545/0.420718。 | 把“cluster 很稳定”与“跨患者 biology 可复现”分开。 | Stage 19 对 State15 的患者组成和最大患者敏感性作专门审计。 |
| 12 | Boundary and normal auxiliary validation | 候选边界和正常参考能否帮助解释 state 的外部邻域？ | normal 与 boundary 仅为辅助；State15 normal mean distance 3.981922，State18 3.023263，State20 3.240648。 | 为后续 State15 anchor 与 manifold 分析提供参照。 | Stage 20–22 保持 normal remote、boundary projection 的辅助地位。 |
| 13 | State atlas and hypotheses | 如何用连续证据汇总 20 states，而不把单一分数当 confidence？ | 形成第一版 20-state atlas 与 10 个 hypothesis candidates。 | 为 Stage15–22 提供固定状态编号和证据背景。 | Stage 24 重新解释这些 labels 的生物学含义，但不改写 atlas artifact。 |
| 14 | Consensus/upstream merge | 如何把 frozen consensus 与 upstream annotation 放在同一逐细胞表？ | 5,378 cells 的合并字段可追溯；Stage 15 audit later confirmed merge inconsistencies=0。 | 支持 Stage 15 candidate identity audit 和后续 post-hoc comparison。 | 外部 candidate/state 结果的完整继承在 Stage 1–6 manifest 中保留。 |
| 15 | Candidate identity audit | 5,378 candidate 是否主要由过宽的 marker-combo gate 产生？ | A merge error=0；FIGF/VEGFD duplicate-pass=0；5,238 marker-combo、140 author/formal；1,443 个 marker-combo candidate 仅有 1-UMI 支持。 | 根因定位到 C：任意两个 marker>0 gate 特异性不足，而非 Stage Modeling 搬错。 | Stage 16 提出连续 identity gate，但没有写回主线 candidate pool。 |
| 16 | Identity gate reconstruction | 能否用 identity anchors+support+competing lineage 重建更可解释的 gate？ | LAM_core 208、boundary 65,930、non_LAM_like 24,503；这是独立诊断 gate，不替换 frozen candidate。 | 使 Stage 17 能够专门追踪外部 positive reference 的漏检。 | Stage 17 发现 GSE190260 的 dataset calibration/dropout/penalty 问题；Stage 24 不继续调 gate。 |
| 17 | Cross-dataset identity calibration audit | GSE190260 为什么漏掉 upstream CORE3-positive？ | GSE190260 2,117 positives：core recovery 0、core+boundary 0.708077；median final score -1.006577、score shift -6.610707、median penalty 4.950320；primary category competing_penalty_only，且 identity/support dropout 均高。 | 证明跨数据集 score scale/marker dropout/penalty 联合作用，不能把外部 0 author-style 当 negative。 | Stage 18 formal LAMCORE anchor validation独立使用 777 genes；不改 Stage16。 |
| 18 | State15 anchor validation | 冻结的 State15 能否作为 LAM-core reference anchor？ | State15 LAMCORE median 0.5125；normal 0.0513；overall author enrichment 13.9408、Fisher p=2.229e-36；decision provisional_reference_candidate_not_formally_upgraded。 | State15 被冻结为后续 reference anchor candidate，而不是正式 classifier。 | Stage19 量化 LAM1163 enrichment；最终仍保持 provisional。 |
| 19 | State15 cross-patient audit | State15 是否只是 LAM1163 这个大患者造成的？ | LAM1163 占 candidate 9.2971%、State15 63.5%，enrichment 6.8301；去除后保留 73 cells、LAMCORE median 0.4683；author labels 仅 GSE135851 可检验。 | 把 State15 定位为患者富集但 profile 可保留的中间结果。 | 为 Stage20 的 anchor-centered geometry 提供组成警示。 |
| 20 | State15-centered latent geometry | 以 State15 为中心，周围细胞的 identity/program 是否沿 latent distance 变化？ | 初始 checkpoint supports_lam_centered_transcriptional_manifold；full LAMCORE 近端→远端 0.2152→0.1556，4/4 dataset rho<0，8/12 patient rho<0。 | 提出 State15-centered manifold 假设，并把验证对象从单个 state 扩展到邻域。 | Stage21 明确削弱 pooled/global manifold；Stage22 改为局部分支 candidate。 |
| 21 | Independent State15 manifold validation | 去掉 anchor、自身 gate/scVI feature overlap 和 composition 后，gradient 是否仍成立？ | 554 genes independent；candidate-only independent slope -0.015466，null empirical two-sided p=0.001996；但 pooled rank rho +0.075665、near/far medians 0.1083/0.1443，患者异质；checkpoint state15_lam_rich_gradient_but_not_robust_manifold。 | 不再支持稳健统一 global manifold，但保留局部/方向性结构的可能。 | Stage22 将问题改成 local branch decomposition。 |
| 22 | State15 local branch decomposition | State15 周围哪些局部方向真正保留 LAM identity？ | 修正后仍入选 State_12, State_16, State_20, State_7；State_16: slope=-0.023931, raw_p=0.878244, BH_q=0.878244, label=ordinary_lineage_adjacency; State_12: slope=-0.022764, raw_p=0.027944, BH_q=0.055888, label=ordinary_lineage_adjacency; State_20: slope=0.010504, raw_p=0.003992, BH_q=0.015968, label=ordinary_lineage_adjacency; State_7: slope=0.001486, raw_p=0.211577, BH_q=0.282102, label=ordinary_lineage_adjacency；Stage22 checkpoint=ordinary_lineage_adjacency_dominates。 | 依据校正后的 local geometry 重新评估分支；不再把 raw empirical p 当作唯一证据，State16 原 transition 标签被降级/撤回（若当前 q 不支持）。 | Stage23 已按新 branch/boundary 输出重生成；Stage24 采用本次修正数字，未改变 State15、X_scVI 或 State1–20 标签。 |
| 23 | Latent-space visualization | 如何把 Stage15–22 已有结构直观看清，而不新增结论？ | 可视化资产完成；2D 优先复用已有 UMAP，3D UMAP/PCA 只用于展示/对照。 | 为 Stage24 artifact provenance 和最终报告提供可引用图件。 | 无；Stage24 不把可视化解释升级为新证据。 |

每个阶段的完整 artifact 路径、行数/文件大小和是否存在见 `artifact_index.csv`；每一阶段的简短可引用版本见 `stage_summary.md`。

## 5. 最终 state evidence

以下表同时保留数值证据和解释性字段。`analogue` 是 human-cell analogue，不是 upstream verified cell_type。State15、16、12、20、7、5 的详细解读见 `state_human_cell_analogue.md`。

| state_id | consensus_cells | consensus_patients | consensus_datasets | analogue | class | identity_LAM_core_fraction | identity_non_LAM_like_fraction | repro_structural_stability | repro_biological_reproducibility | atlas_boundary_connectivity | atlas_normal_distance | atlas_patient_coverage | top_DE_markers | top_program_deltas |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | 703 | 10 | 4 | ciliated/airway-like epithelial | relatively_clear_normal-lineage analogue | 0.0 | 0.9402560455192034 | 0.9628848873718564 | 0.3893193992572874 | 1.1974789915966386 | 3.2717679925022587 | 0.6666666666666666 | ZBBX(8.84), TTC29(8.51), DNAI1(8.23), ADGB(7.95), DNAH9(7.94), VWA3B(7.92), DCDC1(7.78), ARMC4(7.69) | hormone_related=0.301; HOX_PBX=0.268; SLS_stem_like=0.182; IL6_AT2_repair=0.044; TGFbeta_fibroblast=-0.017 |
| 2 | 1 | 1 | 1 | undetermined rare substate | insufficient evidence | 0.0 | 0.0 | 0.5006857804766124 | 0.0 |  |  | 0.0 |  | mTOR_translation=1.114; MDK_dormancy_persistence=0.619; SLS_stem_like=0.256; hypoxia_stress=0.225; IL6_AT2_repair=0.116 |
| 3 | 4 | 2 | 2 | myeloid/inflammatory-like, provisional | insufficient evidence | 0.0 | 0.25 | 0.3414992897788756 | 0.0 |  |  | 0.0 |  | macrophage_TREM2_TYROBP=0.373; IS_inflammatory=0.281; protease_ECM_niche=0.140; hypoxia_stress=0.127; IL6_AT2_repair=0.055 |
| 4 | 550 | 10 | 4 | mixed immune–mesenchymal/interstitial | mixed or uncertain | 0.029090909090909 | 0.1890909090909091 | 0.7964477901985803 | 0.4653757956746926 | 2.0496054114994364 | 3.294472275988629 | 0.75 | IGLV3-21(5.71), TNFRSF17(4.89), CLEC1B(4.57), SDS(4.24), F13A1(4.20), ITGAD(4.10), CD163L1(4.03), MIR137HG(3.98) | CORE2=0.113; CORE3_identity=0.039; protease_ECM_niche=-0.014; cell_cycle=-0.045; TGFbeta_fibroblast=-0.049 |
| 5 | 864 | 10 | 4 | macrophage/myeloid-like | relatively clear normal-lineage analogue | 0.0046296296296296 | 0.3472222222222222 | 0.950408815502366 | 0.559758193072894 | 1.1874139626352016 | 3.075114050883912 | 0.75 | FABP4(5.90), APOC1(5.40), AGRP(5.39), RETN(5.27), MARCO(5.16), CD52(5.11), CCL23(5.10), IGSF6(4.94) | macrophage_TREM2_TYROBP=2.238; mTOR_translation=0.531; protease_ECM_niche=0.493; MDK_dormancy_persistence=0.279; cell_cycle=0.150 |
| 6 | 3 | 2 | 2 | rare myeloid/AT2-mixed substate | insufficient evidence | 0.0 | 0.0 | 0.4086910368110567 | 0.0 |  |  | 0.0 |  | macrophage_TREM2_TYROBP=0.643; mTOR_translation=0.221; MDK_dormancy_persistence=0.203; protease_ECM_niche=0.122; SLS_stem_like=0.060 |
| 7 | 576 | 11 | 4 | AT2-like alveolar epithelial/repair | relatively clear normal-lineage analogue | 0.0052083333333333 | 0.1128472222222222 | 0.6826819842967636 | 0.3978258755809657 | 1.4860161591050345 | 3.329043134515176 | 0.9166666666666666 | PLA2G1B(4.94), AGTR2(4.92), SFTPC(4.75), SFTPA2(4.75), KRT78(4.70), SFRP5(4.70), SFTA2(4.58), SFTPA1(4.55) | IL6_AT2_repair=0.602; hormone_related=0.064; TGFbeta_fibroblast=-0.020; HOX_PBX=-0.032; SLS_stem_like=-0.033 |
| 8 | 1 | 1 | 1 | undetermined TGFβ/interstitial rare substate | insufficient evidence | 0.0 | 0.0 | 0.501159222532781 | 0.0 |  |  | 0.0 |  | TGFbeta_fibroblast=0.501; hormone_related=0.348; HOX_PBX=0.148; MDK_dormancy_persistence=0.115; cell_cycle=-0.071 |
| 9 | 406 | 5 | 4 | mixed repair/interstitial-like | mixed or uncertain | 0.0 | 0.2315270935960591 | 0.7286135596598886 | 0.0 | 1.5862573099415205 | 3.2736464411223016 | 0.25 |  | IL6_AT2_repair=0.374; hormone_related=0.230; TGFbeta_fibroblast=0.164; IS_inflammatory=0.025; hypoxia_stress=-0.025 |
| 10 | 4 | 1 | 1 | rare myeloid/repair mixed state | insufficient evidence | 0.0 | 1.0 | 0.4467594784634767 | 0.0 |  |  | 0.0 |  | macrophage_TREM2_TYROBP=0.647; IL6_AT2_repair=0.245; TGFbeta_fibroblast=0.126; protease_ECM_niche=0.101; CORE2=0.080 |
| 11 | 3 | 1 | 1 | rare fibroblast/HOX-hormone-like state | insufficient evidence | 0.0 | 0.6666666666666666 | 0.4719788426043272 | 0.0 |  |  | 0.0 |  | HOX_PBX=0.495; hormone_related=0.359; TGFbeta_fibroblast=0.311; LAF_niche=0.186; normal_lung_interstitial=0.132 |
| 12 | 605 | 9 | 4 | endothelial/lymphatic endothelial-like | relatively clear normal-lineage analogue | 0.0066115702479338 | 0.3785123966942149 | 0.681160600023974 | 0.3264022951522951 | 2.027136258660508 | 3.4094918755908146 | 0.4166666666666667 | MMRN1(4.53), TM4SF18(4.35), SELE(4.35), UNC5A(4.26), STAB2(4.15), CCL21(4.05), FLT4(3.95), SCN3B(3.76) | IL6_AT2_repair=0.170; hormone_related=0.154; TGFbeta_fibroblast=0.130; HOX_PBX=0.055; biomarker_VEGFD_PMEL_CCL14_MMP8=-0.024 |
| 13 | 6 | 3 | 3 | rare HOX/CORE3-mixed substate | insufficient evidence | 0.1666666666666666 | 0.3333333333333333 | 0.3388228531216182 | 0.0 |  |  | 0.0 |  | HOX_PBX=0.387; hormone_related=0.266; CORE3_identity=0.251; cell_cycle=0.100; TGFbeta_fibroblast=0.079 |
| 14 | 3 | 1 | 1 | rare LAM-myogenic/contractile-like substate | insufficient evidence | 0.0 | 0.6666666666666666 | 0.4193117513949105 | 0.0 |  |  | 0.0 |  | LAM_myogenic_contractile=1.894; lineage_uterine_smooth_muscle=1.894; CORE1=1.374; MDK_dormancy_persistence=0.875; mTOR_translation=0.765 |
| 15 | 200 | 7 | 4 | LAM-rich contractile/mesenchymal candidate; provisional LAM-core anchor | LAM-associated candidate | 0.07 | 0.92 | 0.8541110574273367 | 0.37475757007007 | 4.0 | 3.981922298007541 | 0.4166666666666667 | HTN3(7.12), MMP11(6.47), PAGE4(6.25), PGM5-AS1(6.25), HOXA11(6.10), EMX2(5.64), PLAT(5.63), SFRP1(5.54) | LAM_myogenic_contractile=2.051; lineage_uterine_smooth_muscle=2.051; CORE1=1.456; CORE3_identity=1.145; ECM_remodeling=1.035 |
| 16 | 396 | 10 | 4 | immune/T-NK-adjacent mixed state; no confirmed transition | mixed or uncertain | 0.3762626262626262 | 0.3661616161616162 | 0.7196327799664204 | 0.3230880484706342 | 2.10932944606414 | 3.615899055775924 | 0.5833333333333334 | PCP4(5.73), HOXC10(5.42), TMEM196(5.32), CD3G(5.16), LINC00906(4.99), OR2L5(4.98), LINC00402(4.95), FABP7(4.94) | LAM_myogenic_contractile=1.286; lineage_uterine_smooth_muscle=1.286; CORE1=0.948; CORE3_identity=0.696; CORE2=0.220 |
| 17 | 121 | 5 | 3 | mesothelial/secretory epithelial-like, uncertain | mixed or uncertain | 0.0 | 1.0 | 0.9131796417279768 | 0.4048718447231132 | 1.107981220657277 | 3.6002033286624484 | 0.25 | TMEM151A(10.01), CALB2(9.95), ITLN1(9.90), CPB1(9.45), IL20(8.81), BNC1(8.39), CPA4(8.19), ANXA8(8.02) | mTOR_translation=1.163; ECM_remodeling=0.857; IS_inflammatory=0.509; normal_lung_interstitial=0.499; hypoxia_stress=0.376 |
| 18 | 174 | 9 | 4 | pericyte/VSMC/smooth-muscle-like | relatively clear normal-lineage analogue | 0.0689655172413793 | 0.5344827586206896 | 0.9400330318823972 | 0.4089603116630143 | 1.3327085285848173 | 3.0232627917980324 | 0.6666666666666666 | AC093390.1(8.55), CASQ2(8.41), KCNA5(8.39), COX4I2(7.38), FHL5(7.36), HIGD1B(7.36), FOXC2(7.30), ATP1A2(6.96) | lineage_uterine_smooth_muscle=2.437; LAM_myogenic_contractile=2.437; CORE1=1.934; mTOR_translation=0.590; ECM_remodeling=0.287 |
| 19 | 2 | 2 | 2 | undetermined rare interstitial/hormone-like substate | insufficient evidence | 0.0 | 1.0 | 0.4063892162606694 | 0.0 |  |  | 0.0 |  | hormone_related=0.642; normal_lung_interstitial=0.596; LAF_niche=0.578; TGFbeta_fibroblast=0.452; ECM_remodeling=0.408 |
| 20 | 756 | 10 | 4 | fibroblast/lung interstitial-like | relatively clear normal-lineage analogue | 0.0 | 0.9907407407407408 | 0.9125450360451768 | 0.4207180402657939 | 1.313299232736573 | 3.2406475835121595 | 0.6666666666666666 | PI16(7.90), SFRP2(7.75), SCARA5(6.78), MYOC(6.77), DPT(6.58), C7(6.37), COMP(6.09), CXCL14(5.97) | normal_lung_interstitial=1.701; LAF_niche=1.607; ECM_remodeling=0.996; MDK_dormancy_persistence=0.396; mTOR_translation=0.375 |

### State15 的核心证据链

1. Unsupervised：Stage7 consensus 中出现固定 200-cell State15。
2. Identity：Stage18 formal 777-gene LAMCORE median=0.5125；normal median=0.0513；Stage13 program correspondence 包含 CORE1、CORE3_identity、LAM_myogenic_contractile 和 ECM_remodeling。
3. Author evidence：总体 enrichment fold=13.9408、one-sided Fisher p=2.229e-36；但真实逐细胞 author annotation 只有 GSE135851 可用，49/50 State15 cells 为 author-style positive，外部三个 dataset 为 not_assayed。
4. Patient composition：LAM1163 在 candidate pool 中 9.2971%，在 State15 中 63.5%，enrichment=6.8301；这排除了“只是因为样本多”的解释。
5. Sensitivity：去除 LAM1163 后剩 73 cells，LAMCORE median=0.4683，仍高于指定 comparators；但样本更小，不能替代独立 replication。
6. Structure：Stage11 structural=0.854111、biological=0.374758、patient direction concordance=0.709647、patient coverage=0.416667；这些支持存在 LAM-rich candidate state，同时说明跨患者生物学证据仍有限。
7. Geometry：Stage20 global gradient 后经 Stage21 独立检验被削弱；Stage22 修正后的 local matched-null/FDR 分析不再支持 State16 transition label。因此 State15 最稳妥的称呼是 `provisional LAM-core reference-anchor candidate`。

### State16 及关键 comparator

- State16：396 cells，4 datasets，Stage22 1-hop=96、直接 1-hop 患者数=8（全 state 覆盖 10 patients）；修正后 local independent LAMCORE slope=-0.023931、直接经验 p=0.878244、BH q=0.878244，当前不再是 `LAM_to_lineage_transition_candidate`，而是 State15 邻接/混合状态的探索性描述。
- State18：pericyte/VSMC/smooth-muscle-like，COX4I2/FOXC2/CASQ2 等支持；结构稳定但未进入 Stage22 direct branch selection，保留为重要 LAM mimic/comparator。
- State20：PI16/SFRP2/SCARA5/DPT/C7/COMP/CXCL14，fibroblast/lung interstitial-like；Stage22 ordinary lineage adjacency。
- State12：MMRN1/CCL21/FLT4 等 endothelial/lymphatic-like；ordinary lineage adjacency。
- State7：SFTPC/SFTPA1/SFTPA2 等 AT2-like；ordinary lineage adjacency。
- State5：FABP4/APOC1/MARCO/RETN 等 macrophage-like，尽管有 LAM-shared mTOR/protease/myogenic programs。

## 6. Candidate gate、State15 和 manifold 的完整证据链

Stage15 的 A/B/C 审计是主线必须保留的校准点：annotation/ID merge field inconsistency=0；FIGF/VEGFD alias duplicate-pass=0；5,238 cells 通过 marker-combo、140 cells 有 author/formal support、formal support alone=0；1,443 marker-combo cells 仅由 1-UMI detections 支持。因此“明显非 LAM-like state 被纳入”主要由 C 类 gate 特异性不足解释，而不是 State Modeling 合并搬错。

Stage16 的独立 continuous gate 不使用现有 20 states 调参；其输出为诊断 artifact，不能回写主线。Stage17 进一步显示 GSE190260 positive reference 的 identity/support dropout 与 competing-lineage penalty/score shift 同时存在：2,117 positive references 中 core recovery=0，core+boundary=0.708077，median final score=-1.006577，relative score shift=-6.610707，median competing penalty=4.950320。因此跨数据集零 author-style 不能当作不支持。

Stage20 初始 pooled result 观察到 full LAMCORE 最近/最远 median 0.2152/0.1556、4/4 dataset rho<0、8/12 patient rho<0，提出 global manifold。Stage21 把 State15 排除为 anchor，使用 22,061 non-State15 cells，并审计 777 formal genes：220 与 scVI HVG 重叠，7 属于旧 gate markers，554 同时不属于二者，729 在表达矩阵可用。candidate-only independent slope=-0.015466（95% CI -0.017451,-0.013482，p=1.410338e-51），500 matched fake-anchor empirical two-sided p=0.001996；但 pooled Spearman rho=+0.075665，非-State15 full scope independent rho=+0.094447，距离箱的独立 score 不显示稳健单调下降，且 patient direction heterogeneous。Stage21 因此把结论降为 `state15_lam_rich_gradient_but_not_robust_manifold`。

Stage22 修正后仍只选出 State16、12、20、7 四个方向，但 real/null 都限制在局部 1–3 hop，并进一步匹配距离结构；各分支结果为：State_16: slope=-0.023931, raw_p=0.878244, BH_q=0.878244, label=ordinary_lineage_adjacency; State_12: slope=-0.022764, raw_p=0.027944, BH_q=0.055888, label=ordinary_lineage_adjacency; State_20: slope=0.010504, raw_p=0.003992, BH_q=0.015968, label=ordinary_lineage_adjacency; State_7: slope=0.001486, raw_p=0.211577, BH_q=0.282102, label=ordinary_lineage_adjacency。State16 的旧 transition 标签在直接经验尾部和 BH 校正后不再成立；10,645 的旧 boundary 总体数字也不能直接沿用，当前只报告 6736 / 11648 个 local 1–3-hop boundary unresolved。最终不再支持 State15→State16 的 transition candidate；Stage22 仅保留为校正后的 local adjacency diagnostic，而非统一 global manifold 或 temporal trajectory。

## 7. 阴性结果、被否定假设和异常

| finding_id | topic | observation | stage | source_file | finding_type | implication | disposition |
| --- | --- | --- | --- | --- | --- | --- | --- |
| F01 | candidate gate | 5,238/5,378 candidates entered through marker-combo support; 1,443 had only 1-UMI support. | Stage15 | root_cause_evidence.csv | methodological finding | The original high-recall gate is not a specific LAM classifier; this explains normal-lineage states in the pool. | confirmed by read-only audit |
| F02 | alias audit | FIGF/VEGFD duplicate-pass cells=0 and alias-corrected loss=0. | Stage15 | root_cause_evidence.csv | negative result | Duplicate counting is not the observed source of candidate inflation. | confirmed |
| F03 | upstream annotation | External author-style fields are present as fields but not assayed; only GSE135851 has real positive labels. | Stage19 | author_annotation_availability.csv | data availability | Zeros in three external datasets must not be interpreted as author-negative. | confirmed |
| F04 | identity calibration | GSE190260 CORE3-positive references have 0 core recovery and 0.708077 core+boundary recovery; score shift -6.610707. | Stage17 | root_cause_by_dataset.csv | cross-dataset anomaly | Dataset-aware calibration/dropout/competing penalty must be considered before any future gate revision. | confirmed; no gate change here |
| F05 | patient composition | LAM1163 is 9.2971% of candidate pool but 63.5% of State15; enrichment 6.8301. | Stage19 | state15_patient_composition.csv | composition warning | State15 is not an ordinary pooled state; its profile needs patient-aware interpretation. | confirmed |
| F06 | State15 sensitivity | After removing LAM1163, 73 State15 cells remain and LAMCORE median is 0.4683. | Stage19 | state15_without_LAM1163.csv | supportive sensitivity | The State15 profile is not fully explained by LAM1163, but the remaining sample is small. | confirmed |
| F07 | manifold revision | Stage20 initially supported a pooled State15-centered manifold; Stage21 found independent evidence but not a robust unified manifold. | Stages20-21 | stage20_manifold_report.md; stage21_manifold_validation_report.md | superseded interpretation | Use Stage21 as the later qualification while retaining Stage20 as historical checkpoint. | preserved, not reconciled |
| F08 | local branch | After corrected direct-hop eligibility/local scope/null calibration, selected branches are State_12, State_16, State_20, State_7; State16 label=ordinary_lineage_adjacency (raw p=0.878244, BH q=0.878244). | Stage22 | branch_evidence_summary.csv | corrected structural interpretation | The corrected matched-null/FDR analysis no longer supports the prior State16 transition label; retain local adjacency as exploratory and do not infer a transition. | current; prior label withdrawn |
| F09 | state nomenclature | Stage6 has 12 LAM-only grid reference clusters, Stage7 has 20 consensus states, and the old full cohort has 33 clusters. | Stages6-7 | stage6_checkpoint.json; state_consensus_state_summary.csv | terminology risk | These are distinct clustering objects and must not be called interchangeable. | clarified in glossary |
| F10 | cell type metadata | state_by_cell_type.csv reports unknown for all 5,378 cells. | Stage7 | state_by_cell_type.csv | data gap | Human-cell analogues in Stage24 are expression/program interpretations, not verified upstream cell_type labels. | carried as limitation |
| F11 | normal reference | Normal reference is available in Stage12/18/20 scope but is auxiliary and never defines LAM state count. | Stages12,18,20 | normal_validation.csv; normal_remote_summary.csv | scope boundary | Normal-like/disease-distinct wording remains comparative, not a reclassification. | clarified |
| F12 | formal LAMCORE timing | Stage16 recorded formal 777-gene signature unavailable; Stage18/21 used it after it appeared in data-temp. | Stages16-18 | identity_gate_report.md; state15_anchor_summary.json | historical input change | Do not read Stage16 unavailable as proof that the signature never existed; it was unavailable at that run. | preserved as chronology |
| F13 | DE evidence | Pathway enrichment and regulon outputs are explicit not_available placeholders. | Stage10 | state_pathway_enrichment.csv; state_regulon_summary.csv | method gap | No pathway/regulon conclusion is included in final state analogue calls. | carried as limitation |
| F14 | State16 interpretation | State16 has 80 raw-count LAM/immune coexpressing cells in Stage20, but only 2 LAM-high/immune-high cells. | Stage20 | state16_lam_immune_coexpression.csv | technical/biological ambiguity | Do not call State16 a doublet state; use mixed/transitional wording and retain technical audit caveat. | carried |
| F15 | small states | States 2,3,6,8,10,11,13,14,19 have no supported patients in Stage11 summary. | Stage11 | state_reproducibility_summary.csv | negative result | Their labels are retained for completeness but not used for strong human analogues. | confirmed |
| F16 | unresolved boundary | 6736 of 11648 local 1–3-hop boundary cells remain unresolved in Stage22 local branch assignment; farther boundary cells are not projected. | Stage22 | boundary_local_branch_assignment.csv | negative/uncertain result | Boundary is an evidence-ranking cohort, not a forced new LAM class; report local scope explicitly. | corrected scope |

不能保留的强表述包括：“所有 20 states 都是 LAM”“State15→State16 是时间转化”“所有 State15 邻近 branch 都是 LAM”“candidate gate 可作诊断 classifier”“Stage20 global manifold 已证实”。这些表述分别被 candidate audit、Stage21/22、author availability 和患者复现结果削弱。

## 8. 文字和逻辑审计

| audit_id | issue | severity | evidence | stage24_disposition | audit_type |
| --- | --- | --- | --- | --- | --- |
| A01 | Stage20 checkpoint says supports_lam_centered_transcriptional_manifold; Stage21 says state15_lam_rich_gradient_but_not_robust_manifold. | high | Stage20 report; Stage21 report | Retain both chronologically; final synthesis uses Stage21 qualification and Stage22 local branch result. | historical evolution |
| A02 | Stage6 12 clusters, Stage7 20 consensus states, and full-cohort 33 clusters coexist. | high | stage6_checkpoint.json; state_consensus_state_summary.csv | Add explicit object definitions; do not state that the project has only 12 or only 20 total clusters. | terminology clarification |
| A03 | Stage16 says formal 777 signature unavailable; Stage18 says available. | high | stage16/identity_gate_report.md; stage18/state15_anchor_summary.json | Explain availability changed when data-temp file arrived; Stage16 was not rerun. | chronology |
| A04 | External author-style fields are false/present but not_assayed. | high | stage19/author_annotation_availability.csv | Use not_assayed, never author-negative; enrichment is only formally assessed in GSE135851. | definition correction |
| A05 | Stage21 independent score has negative patient-adjusted slope but positive pooled Spearman and near/far medians not monotonically decreasing. | high | stage21/gradient_models.csv; non_state15_distance_gradient.csv | Report as scope/shape dependence and patient heterogeneity; do not algebraically reconcile as if they were the same estimand. | statistical qualification |
| A06 | Stage22 corrected State16 raw empirical p=0.878244, BH q=0.878244; the prior p=0.025948 came from the superseded scope/null method. | high | stage22/branch_evidence_summary.csv; branch_patient_lopo.csv | Withdraw the prior State16 transition label; current evidence is ordinary-lineage adjacency under corrected local matched-null/FDR analysis. | method correction and conclusion withdrawal |
| A07 | All state_by_cell_type values are unknown, while reports sometimes use human cell analogues. | medium | stage7/state_by_cell_type.csv; Stage13/Stage24 interpretation | Use analogue/inferred lineage terminology, not verified cell_type. | terminology clarification |
| A08 | Stage13 novel_or_unexplained is false for all rows, but this does not mean all states are biologically explained. | medium | stage13/state_atlas.csv | Treat this field as the atlas generation flag, not as proof of identity or absence of uncertainty. | interpretation boundary |
| A09 | Stage12 normal and Stage20 normal_remote use different named scopes. | low | stage12/state_auxiliary_summary.csv; stage20/normal_remote_summary.csv | Describe normal as auxiliary comparison and preserve the stage-specific scope. | scope clarification |
| A10 | Stage10 pathway/regulon files exist but contain not_available placeholders. | medium | stage10/state_pathway_enrichment.csv; state_regulon_summary.csv | Do not imply these analyses were completed. | missing evidence |
| A11 | Stage15 state summaries show FIGF in marker combinations, but duplicate audit is zero. | medium | stage15/state_identity_summary.csv; alias_audit_by_state.csv | Explain that FIGF can appear as the original marker label while no same-cell FIGF+VEGFD double counting was observed. | alias wording |
| A12 | Stage18/19 anchor language remains provisional, not a formal anchor artifact. | medium | stage18/state15_anchor_report.md; state15_anchor_summary.json | Final report uses LAM-core reference-anchor candidate. | strength downgrade |

Stage24 采用的原则是：后来的结果可以限定早期结论，但不删除早期 checkpoint；不同 estimand 的数值不强行调和；缺失的 pathway/regulon 和 upstream cell_type 不补造；formal 777 signature 的可用性按运行时间记录。

## 9. 结论边界与局限

- 数据直接支持：Stage15 是当前 candidate pool 中最 LAM-rich 的 frozen state；它富集 formal LAMCORE/author evidence（仅在可 assay dataset）、并在去除 LAM1163 后保留部分 profile。
- 支持性解释：State15 可能是 LAM-core reference-anchor candidate；State16 与 State15 存在方向性几何邻接，但当前不支持 LAM-preserving/immune-direction transition candidate。
- 尚未证明：State15 是跨患者正式 reference anchor；State15→State16 是时间或谱系转化；存在稳健统一 global manifold；candidate gate 可作为临床/诊断 classifier。
- 患者数量和 composition：State15 只有 7/12 patients，且 LAM1163 enrichment=6.8301；几个小 state 没有 supported patient。
- Dataset heterogeneity：GSE190260 的 identity score shift/dropout 明显；author-style annotation 只在一个 dataset available；scRNA/snRNA、测序深度和转换来源不同。
- Candidate enrichment：主线 high-confidence gate 经 Stage15 证明包含大量普通肺细胞/lineage-like clusters；这既是限制，也让后续 latent space 暴露了 gate 的边界。
- 技术与统计：DE 依赖有限 patient pseudobulk；pathway/regulon 未完成；局部分支 p 值为探索性；distance gradient 不同 scope/estimand 方向不完全一致。
- 生物学验证缺口：没有真实时间维度、空间关系、实验验证、独立 prospective cohort，也没有把 ATAC/spatial 加入本轮验证。

## 10. 未来方向

1. 在独立、跨数据集且 author/formal reference 同步 available 的 cohort 验证 State15 anchor。
2. 对 State16 进行独立实验/空间或正交 modality 验证，区分 LAM-like extension、immune adjacency 与 mixed/doublet profile。
3. 以 patient-aware、dataset-calibrated 的 continuous identity model 重新评估 candidate，而不是恢复任意两个 marker 阳性门槛。
4. 用空间/ATAC/时间序列检验 local branch 是否有真实组织位置、染色质或动态方向。
5. 对普通 lineage comparator（AT2、endothelial、fibroblast、VSMC、macrophage）建立独立参考，减少 candidate enrichment 对 LAM identity 的混淆。
6. 对现有 20 states 的 pathway/regulon 和更大患者数复现补齐，但不在 Stage24 重新训练/聚类。

## 11. 追溯入口

- `stage_index.csv`：Stage1–23 的问题、输入、方法、checkpoint、下一步影响和后续修正。
- `artifact_index.csv`：项目脚本、报告、表格、模型/AnnData 和声明的 upstream paths；包含存在性、大小、行数/类型和小文件 SHA-256。
- `state_human_cell_analogue.csv/md`：20 states 的解释及支持/冲突证据。
- `other_findings_registry.csv/md`：主线外发现、异常、负结果和方法学问题。
- `narrative_audit.csv`：定义变化、数字冲突、历史结论升级/降级和处理方式。
- `appendices/`：方法、state 详细表、统计摘录和 provenance。
