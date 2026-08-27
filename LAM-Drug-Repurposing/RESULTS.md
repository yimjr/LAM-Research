# 首轮执行结果

项目根目录：`.`

## 研究路线与关键转折

本项目最初从 sirolimus/rapamycin residual 出发：利用 GSE179044 的完整因子设计，区分 TSC2-loss、rapamycin 修复、persistent residual、hydrogel-specific residual 和 environment-dependent response，并希望找到能够同时 reversal plastic 与 hydrogel residual 的药物。

实际比较发现，两个环境之间能够稳定、共同印证的 reversal 方向很少。这个结果不表示 plastic 与 hydrogel 没有共同生物学，而是说明“共同 residual reversal”不适合作为候选生成的硬门槛。因此后续候选生成保留 sirolimus 作为参照和机制背景，改为分别用 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 的整体 TSC2-loss disease signature 进行本地 LINCS 匹配。

当前的研究链条是：

```text
GSE179044 residual / environment 机制发现
        ↓
TSC2-loss plastic 或 hydrogel signature → 本地 LINCS 候选生成
        ↓
258 条候选记录 → 去重为 66 个药物
        ↓
药物 gene program 与靶点聚合
        ↓
GSE277844 translation abnormality 与 residual 的分别比较
        ↓
generic stress、人体 LAM 状态、选择性扰动和药理可行性验证
```

因此，下面的 residual 分析、66 个候选药物及其靶点聚合、translation 验证分别回答不同问题，不能把它们合并成一项“已证明的 residual reversal”证据。

## 已完成的科学分析

### GSE179044

- 读取 16 个样本、59,055 个基因，按完整 WT/TSC2−/− × vehicle/rapamycin × plastic/hydrogel 设计计算了 TSC2-loss、rapamycin residual、G×R response、rapamycin 条件下的 hydrogel-specific residual 和 G×R×E environment-dependent escape。
- `signed_residual_ratio=d1/d0` 使用 `|d0|>=0.5` 门槛；保留 absolute ratio 作为残留幅度。效应量来自饱和因子模型，标准误和 FDR 使用 pooled residual-variance moderation。
- hydrogel 中 ratio-eligible 基因分类为：near-complete rescue 1,056；partial rescue + residual 3,099；persistent residual 668；worsened residual 302；direction reversal 1,385。它们是同一实验内的发现分层，不是独立复现。
- curated mTORC1 module 在 hydrogel residual 中接近零、没有 FDR-supported member；相对保留下来的信号更多出现在选定的 myogenic、ECM/invasion、代谢和 autophagy 基因中。这提示 residual 不只是“mTORC1 没有完全压下去”。
- rapamycin 条件下 hydrogel-vs-plastic residual 的探索性信号包括 NNMT、COL8A1、MIR210HG、SLC40A1、FBLN5、DCN 和 LUM，但方向并不统一：部分基因在 hydrogel 中更残留，部分基因反而更接近 rescue。因此它更像 selective state reshaping，而不是普遍增强 resistance。
- 环境依赖 escape 的信号较少，当前只有 8 个基因通过 effect/FDR 筛选，且多数 q 约为 0.097；目前不支持广泛的 hydrogel-dependent escape program。

### GSE27982

- 完成了独立 mouse MEF Tsc2 × rapamycin 2×2 分析，共 13,064 个映射基因/特征。
- 已输出 TSC2-loss、rapamycin 后 residual、genotype-dependent rapamycin response 和 gated signed ratio。
- 与 GSE179044 选择性 top-feature overlap 的符号一致率：TSC2-loss 0.91，rapamycin residual 0.47，genotype-dependent response 0.82。这个结果支持 TSC2-loss 的可重复性强于一个普遍残留程序。
- 低血清条件下 G×R 未被自动升级为 escape。

### GSE16944

- 完成 MMP2/ECM-invasion 历史支持分析。MMP2 的 KO vehicle − addback vehicle 为 35.12，KO rapamycin − KO vehicle 为 8.69（平台为自定义 Codelink，数值只作该实验内描述）。
- 因缺少 TSC2-restored + rapamycin，本数据集没有被用于正式 D(rapamycin) 或 G×R 复现。

### 人体映射

- 已对 GSE135851 的本地快照按 donor × 状态 pseudobulk，并用跨 profile 的 rank score 映射 residual/escape signatures 和功能模块。
- 当前快照只包含 `candidate`/`other` 标签，不含正式 LAMCORE1/2/3 或 LAF-seed/LAF-niche，因此只作为初步人类 LAM-like 支持，不作为 GSE302356 正式验证。
- 新增 GSE302356 的 LAM3/LAM4 scRNA-seq、LAM18 Visium HD 和 LAM20 Visium 原始矩阵后，已完成样本/模态级 residual 与功能模块评分。LAM20 的 ECM、myogenic、metabolism、stress 和 autophagy 分数整体最高；LAM3/LAM4 更偏 survival/stress；LAM18 HD 的总体分数较低。
- 已根据正式论文中报告的 LAMCORE1/2、LAMCORE3 的“无独特 marker”限制，以及 LAF-seed/LAF-niche 的 marker/program 描述，建立 `data/processed/GSE302356/paper_state_marker_panels.csv` 并加入评分。当前结果是 paper-derived operational enrichment，不是作者提供的正式 cell label；LAMCORE3 只作为 shared-core + translation/low-activity surrogate。
- 加入状态面板后，LAM20 的 LAMCORE/LAF 原始分数整体较高，LAM3/LAM4 的单细胞中 LAF-seed/LAF-niche 的高分尾部更明显；但 LAM18/LAM20 与 LAM3/LAM4 跨平台、跨样本，不能把这种差异解释为状态比例或同一患者的配对验证。

### 翻译程序与 residual 分开比较

- 已下载并分析 GSE277844 的 total/polysome counts，使用 conditional polysome-vs-total model 提取 TSC2-null 相对 WT 的 translation-up/down 程序；当前选择结果为 89 个 translation-up 和 107 个 translation-down 基因，属于跨模型探索性结果。
- GSE179044 的 ordinary/plastic persistent residual 与 hydrogel persistent residual 分别比较，没有要求两者先取交集；另外保留 rapamycin 条件下 hydrogel residual 和 hydrogel-specific residual 作为独立类别。
- translation genes 与 plastic persistent residual 的重叠较少（up 2、down 3）；与 hydrogel residual 的重叠更明显，尤其在 effect/FDR 过滤后的 `hydrogel_residual_q10` 中为 up 7、down 13。`hydrogel_specific_residual_q10` 当前没有重叠成员，不能据此推断不存在环境特异的翻译机制。
- 在重叠基因上，RMC-6272 与 eFT-508 都使 `hydrogel_residual_q10` 类别中 18 个 baseline-effect-eligible 基因里的 15 个 KO-vs-WT translation distance 变小；两者 median signed residual ratio 约为 0.52 和 0.50。这个结果只说明 translation distance 的探索性回缩，不能替代独立重复或证明恢复为正常。
- 背景比较显示，全部 196 个 selected translation-abnormal genes 中，RMC-6272 恢复 133/196（67.9%），eFT-508 恢复 147/196（75.0%）；按 |baseline effect|≥0.5 的可比较子集则分别为 128/176（72.7%）和 137/176（77.8%）。因此 hydrogel residual 的 15/18（83.3%）高于两个药的可比较背景，但优势幅度有限且 overlap 样本较小；eFT-508 在不加效应量门槛时为 15/20（75.0%），与全部背景完全相同。plastic persistent residual 则不呈现统一优势：RMC-6272 为 4/5，eFT-508 为 2/5。
- 逐基因结果显示，`hydrogel_residual_q10` 的 18 个可比较重叠基因中，13 个同时被两种药拉近 WT（CACFD1、CDC42EP3、FBN2、GPC4、GPR27、NFATC4、PNMA2、REEP2、RND3、SERPINE2、SPIN4、WWTR1、ZNF354C），2 个主要由 RMC-6272 支持（FIBIN、HOXC6），2 个主要由 eFT-508 支持（APOL1、TYMS），1 个两者都未支持（ZWINT）。两种药的恢复集合 Jaccard 为 0.765，提示主要作用于同一批基因，但并非完全相同。
- plastic persistent residual 的 5 个可比较基因中，GPR27 和 RNF182 被两种药共同支持，HOXC6 与 NETO1 仅被 RMC-6272 支持，ALCAM 两者均未恢复；该类别数量过小，不能据此判断药物机制差异。
- 当前分析明确使用 anota2seq-like conditional model，而非官方 anota2seq 完整流程；GSE277844 是人源 NPC 模型，不是 LAM 直接复现。

### 固定 13 基因的功能注释

- 已固定筛选对象：`residual_category=hydrogel_residual_q10`、两种药均 `baseline_effect_eligible=True`、且 `both_drugs_distance_reduced=True`，共 13 个基因：CDC42EP3、RND3、SERPINE2、GPR27、WWTR1、FBN2、REEP2、ZNF354C、PNMA2、CACFD1、NFATC4、SPIN4、GPC4。
- 已分别查询 GO Biological Process、GO Cellular Component、Reactome 和 MSigDB Hallmark，并输出每个基因的完整功能条目，而不是只报告一个 GO 富集表。
- 多基因重复出现的解释性主题包括 ECM（SERPINE2/FBN2）、cell adhesion（RND3/SERPINE2/GPC4）、Rho GTPase（CDC42EP3/RND3）、TGF-β（WWTR1/FBN2）和 migration（CDC42EP3/SERPINE2）。actin cytoskeleton、focal adhesion、Hippo/YAP/TAZ 和 mechanotransduction 目前主要由单个基因支持，暂不作为多基因模块结论。
- 13 个基因较小，term-level FDR 仅作描述性参考；当前更重要的是多个基因是否指向同一主题，不能因单个主题 FDR 不显著而直接否定。

### 候选药物

- 已生成 8 类默认 CMap 查询 signature；direction reversal 只进入诊断表，不进入默认查询。
- generic cytotoxicity、target KD/KO concordance、相关状态表达、无关疾病 promiscuity、人体暴露和 sirolimus 互补性均已纳入候选过滤接口。
- GSE92742/GSE70138 的 Level 5 GCTX 与配套 metadata 已从 GEO 下载、SHA512 校验并移入 `data/raw/LINCS/`；本项目的 connectivity、WTCS/NCS 和候选分析均由 `scripts/analyze_lincs_local.py` 在本地完成。
- 两套 LINCS release 只用于 `cross-phase/cross-release recurrence`，不作为在线提交结果或独立生物学复现；本阶段不依赖在线查询服务。
- 本地 LINCS/CMap 计算已完成：两套 release 共构建 21 个去重 query；GSE92742 输出 9,483,222 条 signature-query 结果，GSE70138 输出 2,337,006 条结果。评分使用完整 10,174-gene BING space、weighted KS/WTCS、published NCS 和 weighted-correlation sensitivity analysis。
- 已生成 context-level、normalized perturbation-level、cross-phase/cross-release recurrence 和 positive-control sanity 表。recurrence 分类为：`replicated_concordant` 283、`replicated_discordant` 29、`replication_available_but_weak` 22,410、`replication_not_available` 599,970。未测到不被当作复现失败；两边都测到但方向冲突才标记 discordant。
- mTOR/rapamycin panel 已作为 biological sanity check 保留；它没有被设为算法硬验收条件。当前仍不生成 Tier 1 候选，所有候选只保留为后续人体状态、generic cytotoxicity、target concordance、暴露和外部模型过滤的输入。

### 新增机制输入

- GSE84476 已用 GENCODE v24 transcript map 汇总到 gene-level，并生成 102/103 cell context 下的 STAT3 knockdown、rapamycin 和两者比较表。
- GSE104335 的 9 个基因表达样本已从 HTA-2_0 `sst-rma-gene-full` CHP 直接解析为 23,910 个 gene-level features，并用 3 个生物学重复/组做 limma moderated contrasts。这个路径使用的是归档内已经由 Expression Console 生成的 RMA gene-level 值，而不是把 CHP 与重新运行的 CEL-RMA 混在一起。
- 在 shGFP 细胞中，rapamycin 后 COL8A1、NNMT、DCN、ACTA2 和 MMP2 上升，而 HMGCS1、NUPR1、HSPA5 下降；这与“rapamycin 后 ECM/代谢残留并不等于单纯 mTORC1 未抑制”的主线相容，但仍是机制支持，不是 TSC2×rapamycin 因子复现。
- SRPK2 knockdown 强烈降低 SRPK2，并降低 COL3A1、NNMT、LUM、FBLN5；同时升高 DDIT3/ATF4。尤其 NNMT 与 COL8A1 相对 rapamycin 的方向相反，提示 SRPK2 可能触及 GSE179044 中 hydrogel residual 的 ECM–metabolic axis，但它不是简单的 rapamycin mimic。由于没有 `shSRPK2 + rapamycin` 组，不能从本数据声称联用效应。
- 在同时满足 GSE104335 rapamycin response 与 GSE179044 hydrogel-specific residual 的 FDR<0.05 条件下，共有 19 个重叠基因，其中 11 个方向一致、8 个方向相反。NNMT 和 COL8A1 是最清楚的同向组合；FAP、SLC40A1 等相反方向的例子说明该轴不是普遍 ECM resistance，而是状态选择性重塑。
- `GSE104335_cross_dataset_summary.csv` 和 `GSE104335_hydrogel_specific_overlap.csv` 已把这些机制对比与 GSE179044 的 hydrogel residual、hydrogel-specific residual 和 escape contrast 对齐，便于后续筛选真正跨模型的候选。

### TSC2-loss signature 直接匹配得到的候选药物

- 按 `contrast ∈ {tsc2_loss_plastic, tsc2_loss_hydrogel}`、`perturbation_class=compound`、`cross_phase_status=replicated_concordant` 筛得 258 行；去除 dataset/query-size 重复后为 66 个唯一药物。
- 258 行中 194 行为 `reversal_direction`，64 行为 `mimic_direction`。因此 `replicated_concordant` 表示两个 LINCS release 的方向类别一致，不表示全部是候选治疗药。
- 66 个唯一药物中，54 个只呈 reversal 方向，12 个只呈 mimic 方向；当前不据此生成 Tier 1。
- 原有 92 行/29 个药物及其靶点表保留为 `tsc2_loss_plastic` 单环境历史结果；合并 plastic/hydrogel 后的 66 个药物已用新脚本完成 ChEMBL、PubChem BioAssay 和 BindingDB 靶点/外部证据整理。当前新表包含 30 个药物的 ChEMBL curated mechanism、12 个药物的文献 fallback，以及其余无 curated mechanism 的明确记录；不同证据层级不混同。
- 66 个药物的 LINCS 基因程序分析也已完成：两个 release 分别汇总，跨 release 只比较 204 个共同可分析基因；这是候选范围扩展后的机制线索，不是 Tier 1 结论。
- 同一批 66 个药物已进一步用 `tsc2_loss_hydrogel` 签名独立分析，并与 plastic 结果比较：两套 top150+top150 面板同方向重叠 165/300，无 top-panel 方向相反重叠；主要比较改为逐药物、逐 LINCS release 的 reversal gene set：median 交集为 33.5 个基因、median Jaccard 为 0.374，且每个 panel 都有约 135 个面板特异基因。共同基因上的 LINCS drug effect 本来来自同一 perturbation，不作为独立生物学证据。
- 槲皮素仍未进入合并候选表，因为它在两个 release 中仍属于 `replication_available_but_weak`，而不是 `replicated_concordant`。

## 关键输出

### 可随项目提供的结果与文档

- [GSE179044 因子 contrast 表](results/tables/GSE179044_factorial_contrasts.csv)
- [GSE179044 功能模块总结](results/tables/GSE179044_functional_module_summary.csv)
- [GSE179044 方向反转诊断](results/tables/GSE179044_direction_reversal_diagnostics.csv)
- [GSE27982 外部验证](results/tables/GSE27982_external_response.csv)
- [GSE16944 历史支持](results/tables/GSE16944_historical_support.csv)
- [CMap 查询 signatures](results/signatures/GSE179044_cmap_query_signatures.csv)
- [GSE135851 人体状态映射](results/human_mapping/GSE135851_state_summary.csv)
- [GSE302356 原始样本程序评分](results/human_mapping/GSE302356_raw_sample_scores.csv)
- [GSE84476 gene-level 机制比较](results/mechanisms/GSE84476_gene_level_log2_tpm_contrasts.csv)
- [GSE104335 gene-level 机制比较](results/mechanisms/GSE104335_gene_level_contrasts.csv)
- [GSE104335 × GSE179044 对齐表](results/mechanisms/GSE104335_cross_dataset_summary.csv)
- [GSE104335 hydrogel-specific overlap](results/mechanisms/GSE104335_hydrogel_specific_overlap.csv)
- [当前缺口与最小补充方案](MISSING_INPUTS.md)
- [候选过滤表](results/candidates/candidate_filtering_table.csv)
- [LINCS perturbation-level 汇总](results/candidates/LINCS_perturbation_summary.csv.gz)
- [LINCS positive-control sanity check](results/candidates/LINCS_positive_control_validation.csv)
- [plastic/hydrogel concordant compound 候选表](results/candidates/tsc2_loss_plastic_or_hydrogel_replicated_concordant_compounds.csv)
- [plastic/hydrogel concordant compound 去重表](results/candidates/tsc2_loss_plastic_or_hydrogel_replicated_concordant_compounds_unique.csv)
- [plastic/hydrogel 候选筛选 manifest](manifests/tsc2_loss_plastic_or_hydrogel_replicated_concordant_analysis.json)
- [TSC2-loss plastic concordant compound 原始 92 行](results/candidates/tsc2_loss_plastic_replicated_concordant_compounds_92_rows.csv)
- [TSC2-loss plastic concordant compound 去重 29 药物](results/candidates/tsc2_loss_plastic_replicated_concordant_compounds_29_unique.csv)
- [候选生成结果说明](results/candidates/README.md)
- [GSE179044 发现报告](reports/01_discovery/GSE179044_discovery.md)
- [GSE27982/GSE16944 外部验证报告](reports/02_external_validation/GSE27982_GSE16944_validation.md)
- [STAT3/SRPK2 机制报告](reports/03_mechanisms/GSE84476_GSE104335_mechanisms.md)
- [人体 LAM 映射报告](reports/04_human_mapping/GSE135851_GSE302356_mapping.md)
- [LINCS 候选生成报告](reports/05_candidate_generation/LINCS_candidate_generation.md)
- [LINCS 候选后续分析报告](reports/06_candidate_analysis/LINCS_candidate_analysis.md)
- [合并 plastic/hydrogel 候选后续分析报告](reports/06_candidate_analysis/tsc2_loss_plastic_or_hydrogel_replicated_concordant_LINCS_gene_program_analysis.md)
- [合并候选 hydrogel panel 分析报告](reports/06_candidate_analysis/tsc2_loss_plastic_or_hydrogel_replicated_concordant_hydrogel_panel_LINCS_gene_program_analysis.md)
- [plastic 与 hydrogel 签名比较报告](reports/06_candidate_analysis/tsc2_loss_plastic_or_hydrogel_replicated_concordant_plastic_vs_hydrogel_comparison.md)
- [阶段报告总览](reports/README.md)
- [候选后续分析说明](candidate_analysis/README.md)
- [候选后续分析报告](reports/06_candidate_analysis/LINCS_candidate_analysis.md)
- [LINCS 分析 manifest](manifests/lincs_analysis.json)
- [GSE277844 翻译程序分析说明](translation_analysis/README.md)
- [GSE277844 翻译程序阶段报告](reports/07_translation_analysis/GSE277844_translation_residual_analysis.md)
- GSE277844 translation effects（`data/processed/translation_analysis/GSE277844_tsc2_loss_translation_effects.csv`）
- GSE277844 与 residual 重叠汇总（`data/processed/translation_analysis/GSE277844_translation_residual_overlap_summary.csv`）
- GSE277844 translation-targeting drug summary（`data/processed/translation_analysis/GSE277844_translation_residual_overlap_drug_summary.csv`）
- [GSE277844 分析 manifest](manifests/GSE277844_translation_residual_analysis.json)
- [GSE277844 背景恢复率与药物一致性研究日志](research_log/2026-08-27_translation_background_and_drug_concordance.md)
- [GSE277844 固定 13 基因功能注释报告](reports/07_translation_analysis/GSE277844_hydrogel_translation_core_functional_annotation.md)
- [GSE277844 固定 13 基因功能注释 manifest](manifests/GSE277844_hydrogel_translation_core_functional_annotation.json)
- [GSE277844 固定 13 基因功能注释研究日志](research_log/2026-08-27_hydrogel_translation_core_functional_annotation.md)
- [研究灵感记录](research_log/2026-08-22_initial_findings.md)

### 运行后生成的未版本化结果

以下文件由分析运行后生成，不随版本库提供；远端重新获取项目后，需要重新运行相应脚本生成。

- GSE302356 paper-derived 状态面板（`data/processed/GSE302356/paper_state_marker_panels.csv`）
- LINCS signature-level WTCS/NCS（GSE92742）（`results/candidates/GSE92742_LINCS_signature_WTCS.parquet`）
- LINCS signature-level WTCS/NCS（GSE70138）（`results/candidates/GSE70138_LINCS_signature_WTCS.parquet`）
- LINCS context-level 汇总（`results/candidates/LINCS_context_summary.csv.gz`）
- LINCS cross-phase/cross-release recurrence（`results/candidates/LINCS_cross_dataset_recurrence.csv`）
- 29 药物主靶点表（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_compound_targets.csv`）
- 靶点汇总（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_target_summary.csv`）
- 靶点家族汇总（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_target_family_summary.csv`）
- 药物–靶点分析（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_drug_target_analysis.csv`）
- PubChem 化合物身份与描述（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_pubchem_identity.csv`）
- PubChem BioAssay target evidence（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_pubchem_assay_target_evidence.csv.gz`）
- BindingDB affinity evidence（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_bindingdb_affinity_evidence.csv.gz`）
- BindingDB target summary（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_bindingdb_target_summary.csv`）
- 合并 plastic/hydrogel 的 66 药物靶点与外部证据表（`data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_or_hydrogel_replicated_concordant_*`）
- 合并候选的 LINCS drug×gene、跨 release、聚类和模块结果（`data/processed/candidate_analysis/programs/tsc2_loss_plastic_or_hydrogel_replicated_concordant_*`）
- 合并候选 hydrogel panel 的 LINCS drug×gene、跨 release、聚类和模块结果（`data/processed/candidate_analysis/programs/tsc2_loss_plastic_or_hydrogel_replicated_concordant_hydrogel_panel_*`）
- plastic/hydrogel 逐基因、逐药物和签名面板比较结果（`data/processed/candidate_analysis/programs/tsc2_loss_plastic_or_hydrogel_replicated_concordant_plastic_vs_hydrogel_*`）
- plastic/hydrogel reversal gene set 明细及按药物跨 release 汇总（`data/processed/candidate_analysis/programs/tsc2_loss_plastic_or_hydrogel_replicated_concordant_plastic_vs_hydrogel_drug_reversal_*`）
- 合并候选的 gene-panel/signature audit 与 target-axis validation（`data/processed/candidate_analysis/audit/tsc2_loss_plastic_or_hydrogel_replicated_concordant_*`、`data/processed/candidate_analysis/validation/tsc2_loss_plastic_or_hydrogel_replicated_concordant_*`）
- 合并候选分析 manifest（`manifests/tsc2_loss_plastic_or_hydrogel_replicated_concordant_LINCS_gene_program_analysis_manifest.json`）
- GSE277844 translation analysis 的运行后结构化输出（`data/processed/translation_analysis/`）和报告中的结果表；远端重新获取项目后，需要重新运行 `scripts/analyze_gse277844_translation_residuals.py` 生成。
- GSE277844 residual 恢复率背景比较（`data/processed/translation_analysis/GSE277844_translation_residual_recovery_background_comparison.csv`）
- GSE277844 逐基因药物一致性（`data/processed/translation_analysis/GSE277844_translation_residual_drug_gene_concordance.csv`、`data/processed/translation_analysis/GSE277844_translation_residual_drug_concordance_summary.csv`）
- GSE277844 固定 13 基因的逐基因功能注释、主题汇总和功能富集（`data/processed/translation_analysis/GSE277844_hydrogel_translation_core_*`）；远端重新获取项目后，需要重新运行 `scripts/annotate_gse277844_hydrogel_translation_core.py` 生成。

## 当前主要限制

1. GSE179044 每格 n=2，moderated variance 只能稳定不确定性，不能替代更多生物学重复。
2. GSE27982 的低血清环境改变 WT 的 mTORC1 基线，不能把所有 G×R 解释成 escape。
3. GSE84476 gene-level 已完成，但样本数和原始设计有限，机制结果主要作支持性证据。
4. GSE104335 已完成 CHP gene-level 解析，但其设计没有 `shSRPK2 + rapamycin`，因此 SRPK2 与 rapamycin 只能做机制方向比较，不能估计联用 interaction。
5. 正式 LAMCORE1/2/3、LAF-seed/LAF-niche 的 cell-level 标签及其空间/ATAC 对应关系仍缺失；当前 paper-derived panel 只能提高可解释性，不能替代作者的 processed metadata。
6. 两个 LINCS 数据集属于同一 L1000 体系的不同阶段/release，因此 cross-release recurrence 的证据强度低于 GSE27982、人体状态和换平台验证。
7. 当前 LINCS 阶段只完成连接性、上下文和 release 稳健性分析；generic cytotoxicity、target KD/KO concordance、人体暴露和正式 LAMCORE/LAF 状态过滤尚未完成，因此不提出 Tier 1 联合治疗药物。
