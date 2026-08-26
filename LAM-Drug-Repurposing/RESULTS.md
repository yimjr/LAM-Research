# 首轮执行结果

项目根目录：`.`

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

### 候选药物

- 已生成 8 类默认 CMap 查询 signature；direction reversal 只进入诊断表，不进入默认查询。
- generic cytotoxicity、target KD/KO concordance、相关状态表达、无关疾病 promiscuity、人体暴露和 sirolimus 互补性均已纳入候选过滤接口。
- 曾加入 CLUE/CMap API 连接器 `scripts/clue_api.py` 和本地查询计划 `results/cmap/clue_query_plan.json`；随着 CLUE 退役，它们仅作为历史记录保留，不再执行在线提交。
- CLUE 已于 2026-01-31 退役；GSE92742/GSE70138 的 Level 5 GCTX 与配套 metadata 已从 GEO 下载、SHA512 校验并移入 `data/raw/LINCS/`。下一步改为本地 connectivity 计算，不再依赖在线 CLUE API。
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

### `tsc2_loss_plastic` concordant compound 子集

- 按 `contrast=tsc2_loss_plastic`、`perturbation_class=compound`、`cross_phase_status=replicated_concordant` 精确筛得 92 行；去除 dataset/query-size 重复后为 29 个唯一药物。
- 92 行中 60 行为正向 `reversal_direction`，32 行为方向一致但属于 `mimic_direction`。因此 `replicated_concordant` 表示两个 LINCS release 的方向类别一致，不表示全部是候选治疗药。
- 29 个唯一药物中，20 个只呈 reversal 方向，9 个只呈 mimic 方向；当前不据此生成 Tier 1。
- 靶点主表覆盖全部 29 个药物、188 条记录：161 条来自 ChEMBL curated mechanism，27 条来自明确标记的 primary-literature fallback。ChEMBL 缺失机制不再被误读为“无靶点”。
- 主要聚集轴是 PI3K/AKT（7 个药物）与 mTOR（7 个药物），其次是 proteasome（3 个）、microtubule（2 个）、受体/转运体/通道（6 个）和其他机制。PI3K/mTOR 轴的集中出现更像 TSC2-loss 相关状态的生物学 sanity signal，但不能排除共同的强转录扰动或 cytotoxicity。
- PubChem 已完成 29/29 化合物身份解析，并保留 21,292 条带 target accession/GeneID 的 BioAssay evidence；这些记录单独保存，不直接当作主要靶点。
- BindingDB 高相似度检索保留 1,177 条 ligand-target affinity evidence。由于其中包含弱结合、面板筛选和物种混杂，BindingDB 只用于直接结合/脱靶背景，不直接改变主靶点等级。

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
- [TSC2-loss plastic concordant compound 原始 92 行](results/candidates/tsc2_loss_plastic_replicated_concordant_compounds_92_rows.csv)
- [TSC2-loss plastic concordant compound 去重 29 药物](results/candidates/tsc2_loss_plastic_replicated_concordant_compounds_29_unique.csv)
- [候选生成结果说明](results/candidates/README.md)
- [GSE179044 发现报告](reports/01_discovery/GSE179044_discovery.md)
- [GSE27982/GSE16944 外部验证报告](reports/02_external_validation/GSE27982_GSE16944_validation.md)
- [STAT3/SRPK2 机制报告](reports/03_mechanisms/GSE84476_GSE104335_mechanisms.md)
- [人体 LAM 映射报告](reports/04_human_mapping/GSE135851_GSE302356_mapping.md)
- [LINCS 候选生成报告](reports/05_candidate_generation/LINCS_candidate_generation.md)
- [LINCS 候选后续分析报告](reports/06_candidate_analysis/LINCS_candidate_analysis.md)
- [阶段报告总览](reports/README.md)
- [候选后续分析说明](candidate_analysis/README.md)
- [候选后续分析报告](reports/06_candidate_analysis/LINCS_candidate_analysis.md)
- [LINCS 分析 manifest](manifests/lincs_analysis.json)
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

## 当前主要限制

1. GSE179044 每格 n=2，moderated variance 只能稳定不确定性，不能替代更多生物学重复。
2. GSE27982 的低血清环境改变 WT 的 mTORC1 基线，不能把所有 G×R 解释成 escape。
3. GSE84476 gene-level 已完成，但样本数和原始设计有限，机制结果主要作支持性证据。
4. GSE104335 已完成 CHP gene-level 解析，但其设计没有 `shSRPK2 + rapamycin`，因此 SRPK2 与 rapamycin 只能做机制方向比较，不能估计联用 interaction。
5. 正式 LAMCORE1/2/3、LAF-seed/LAF-niche 的 cell-level 标签及其空间/ATAC 对应关系仍缺失；当前 paper-derived panel 只能提高可解释性，不能替代作者的 processed metadata。
6. 两个 LINCS 数据集属于同一 L1000 体系的不同阶段/release，因此 cross-release recurrence 的证据强度低于 GSE27982、人体状态和换平台验证。
7. 当前 LINCS 阶段只完成连接性、上下文和 release 稳健性分析；generic cytotoxicity、target KD/KO concordance、人体暴露和正式 LAMCORE/LAF 状态过滤尚未完成，因此不提出 Tier 1 联合治疗药物。
