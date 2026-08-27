# LAMCORE 复现基础与新生物学探索计划

## 1. 目标与结果定位

本项目的核心问题是：

> 在不预先排除弱 marker、低活性或非典型 LAMCORE-like 细胞的前提下，寻找跨患者重复、且不能被现有 LAM 状态框架充分解释的表达程序。

项目结构分为三个层次：

1. 论文核心结果复现：建立与原作者分析相对应的可信基线；
2. 必要的稳健性检验：确认基线和重要结果不是单一参数、donor 或技术因素造成；
3. 新生物学探索：以新的生物学问题为主线，寻找值得进一步实验验证的现象、关联和机制线索。

因此整体顺序是：

> 核心复现 → 必要稳健性检验 → 以探索为主线的后续研究

进入探索阶段后，不再把大规模方法比较作为项目终点；只有当某个具体候选需要确认时，才补充相应的稳健性验证。

此前 `results/independent_reanalysis/` 和旧版 `results/discovery/` 的结果继续保留为同一数据集上的独立再分析，不覆盖，也不直接称为原论文复现。未知程序发现阶段的任何候选，在获得独立 PatientID 和正交证据前都只称为探索性结果或研究线索。

## 当前执行状态（更新至 2026-08-27）

论文核心复现和必要稳健性检验已经完成到可以支撑后续研究的程度：核心基线已冻结，已知程序 benchmark 和 baseline/strict-QC 检查已经建立，GSE179044/GSE84476 扰动数据、GSE302356 三模态空间数据和肺适应参考比较已经完成首轮处理。

当前重点已经转向新生物学探索。首轮结果主要集中在三个相互关联的问题：肺部 protease–antiprotease 空间生态位、rapamycin 后仍保留的 ECM/protease 程序，以及 LAM 细胞进入肺后获得的适应程序。另有未知程序发现、LAMCORE 状态异质性和血浆/EV 蛋白连接两条扩展线。所有结果目前仍按探索性假说或高价值候选机制报告，不直接称为已证实的新机制。

HDBSCAN 分支对原有 KMeans 状态异质性探索进行了独立补充：在 140 个候选细胞中，主设定得到 2 个局部密度簇、28 个非噪声细胞，且只涉及 LAM3/LAM4；目前不支持覆盖四个 donor 的稳定离散亚型，更适合继续研究连续状态和局部状态结构。

已知程序 benchmark 仍提示 top-gene matching 的灵敏度不均匀，严格 QC 将候选数从 140 降至 85，因此每张 Hypothesis Card 仍需报告 baseline/strict-QC 结果。但这些检查现在作为探索质量控制，而不是项目的主要终点。详细状态见 `results/report/parallel_research_status_zh.md` 和英文版。

## 2. 数据与独立性规则

核心肺部基线使用 GSE135851 的 LAM1–LAM4 和 Donor1；外部优先使用 GSE190260、GSE217108 和 GSE302356。LAM uterus、正常/小鼠 uterus 用作谱系参考，正常肺间质用作疾病对照，不作为统一阴性淘汰条件。

独立性以实际 `patient_id` 为主键，而不是 accession、SampleID 或测序条目：

- GSE217108 只有 LAM32、LAM44 两个 LAM donor；同一患者的 RNA 与 ATAC 只增加正交证据，不增加 donor 数。
- GSE302356 按论文给出的 PatientID 映射：SampleID LAM3→PatientID LAM32，LAM4→LAM18，LAM10/LAM14/LAM18→LAM3，LAM13→LAM50，LAM15/LAM19/LAM20→LAM4；LAM16 仍待元数据核实。
- GSE190260 分开记录 donor、specimen 和 sample；LAM1164-1/2/3 在确认前视为一个潜在 donor，不能把 sample 条目数当作患者数。
- GSE139819 与 GSE139534 按同一研究体系记录，不能算两个独立验证集。

完整映射见 `manifests/donor_registry.yaml`。每条证据同时记录 `patient_id`、`donor_id`、`specimen_id`、`sample_id`、`assay`、`tissue`、`treatment` 和 `independence_group`。

只使用公开处理后数据，不下载 FASTQ，不把受控访问的原始数据作为主流程依赖。由于从处理后矩阵开始，QC 的准确表述是：

> 在处理后矩阵允许范围内恢复下游 QC。

不能声称恢复 FASTQ、Cell Ranger、初始 barcode/cell calling 或 empty-droplet 判断。

## 3. 环境和运行约束

- Python 3.12.13：项目 `.venv`；所有包通过 `.venv/bin/python -m pip` 管理。
- 不使用 conda，不向系统 Python 安装包。
- 后续独立分析优先使用 Python；若核心复现需要作者 R/Seurat，则使用隔离的 R 环境。
- 固定随机种子，所有运行参数写入 `results/program_discovery/program_discovery_run_manifest.json` 和项目 run manifest。

## 4. 两层候选池与身份保护

### 4.1 高置信 LAMCORE-like 集

用于可靠的 marker 描述、donor-level pseudobulk 以及与 CORE1/2/3 和原论文的比较。证据可来自 PMEL、MLANA、MITF、CTSK、FIGF/VEGFD、ACTA2、ESR1 等组合、聚类位置、外部参考映射、空间或 ATAC 支持。

### 4.2 宽松 LAM-like 集

用于发现非典型状态，不能要求 PMEL 或 ACTA2 高表达。候选可由部分 LAM marker、LAM-associated mesenchymal/lineage 特征、TSC1/TSC2 相关遗传证据（如数据可用）、跨数据集邻近关系或 LAM 相关空间定位产生。

TSC2 低表达本身不能证明 TSC2 缺失。宽松候选在程序发现后再通过 marker 组合、遗传证据、空间位置、外部 donor 和正常肺/子宫/LAF 对照验证。身份无法确认但程序有意义的细胞保留为 `uncertain_identity_hypothesis`，不强行归入 LAMCORE，也不静默删除。

为了检验“候选池本身把新状态筛掉”的风险，脚本另保留一个不依赖 marker 的 `unrestricted_lam` guardrail pool，即核心数据中全部 LAM 条件细胞；它用于诊断，不直接当作 LAMCORE 结论。

doublet 在第一阶段只计算并记录，不默认删除；doublet 去除与不去除属于稳健性验证。

## 5. 两种程序发现方法

### 5.1 合并数据发现

高置信、宽松和 unrestricted pool 分别进行 consensus NMF，多随机种子、多模块数、donor 平衡抽样和 assay 分层。合并 NMF 只用于提出候选，不能单独证明跨患者稳定。

### 5.2 donor-wise 独立发现

每个 donor 单独进行 NMF/模块提取、donor 内差异表达、程序评分和子抽样稳定性分析，然后用核心基因重叠、gene-set similarity、方向一致性和调控网络匹配构建 meta-program。

优先级最高的候选必须在不同患者中被独立重新发现，而不是一个 pooled 程序在多个 donor 中被动打分。细胞数不足的 donor 使用 donor-level pseudobulk，不能强行独立 NMF。

当前实现：`scripts/discover_lam_programs.py` 已生成 pooled program、donor-wise program、meta-program matching、donor-level score 和 residual sensitivity 表。`meta_program_summary.csv` 只统计 donor 独立发现后达到当前匹配标准的候选。

## 6. 已知程序比较层

已知程序保存在 `config/known_lam_programs.yaml`，每条记录至少包含：

- `program_name`；
- `evidence_scope`；
- `evidence_level`；
- `source_study`；
- `gene_set_or_model`；
- `known_lamcore_relation`；
- `treatment_relation`；
- `microenvironment_relation`。

证据范围区分：

- `LAMCORE-specific`：CORE1、CORE2、CORE3、777-gene LAMCORE；
- `LAM-shared`：不同 LAM 细胞或组织共有程序；
- `TSC-tumor-shared`：LAM、AML 或其他 TSC/mTOR 肿瘤共有；
- `microenvironment`：LAF、巨噬细胞、淋巴管、AT2 等；
- `treatment-associated`：SLS、IS、MDK、dormancy、rapamycin persistence；
- `biomarker`：VEGF-D、PMEL、CCL14、MMP8 等；
- `lineage/reference`：子宫、平滑肌、肺间质参考。

候选与 SLS 等程序部分重叠时不直接淘汰，而是比较：已知程序解释比例、LAM-specific 部分、独有基因、TF/regulon、空间/ATAC/蛋白支持和 donor 复现。主分析不预先回归 ECM、mTOR、SLS、激素或其他生物学程序；回归只作为敏感性分析。

## 7. CORE3 特殊模型

CORE3 不作为普通 gene set，而作为三部分结构化模型：

1. CORE1-like identity；
2. 经 `n_genes_by_counts`、UMI、assay 校正后的低转录活性；
3. protein translation enrichment。

单纯低 UMI、低检测基因数或低 complexity 不能定义 CORE3-like。若低活性在深度校正后消失，优先解释为技术因素；若校正后保留并在外部 donor 重复，才作为低活性状态扩展继续研究。

## 8. 验证与证据分级

综合使用 donor 重复性、独立研究体系、方法敏感性和正交证据。建议的高可信候选条件是：至少 3 个实际 PatientID、至少 2 个研究体系、至少 2 个 donor 独立发现相似 meta-program、至少一个外部体系复现，且不能由单一 donor、assay、批次、doublet 或低质量驱动，并有 ATAC、空间或蛋白中的至少一种支持。

这些是优先级标准，不是机械淘汰规则。结果分为：

1. 已知程序；
2. 已知程序的新扩展；
3. 跨 donor 的高可信新程序；
4. 探索性假说；
5. 技术或身份不确定结果。

## 9. 当前主线：新生物学探索

完成核心复现和必要稳健性检验后，项目当前的主要工作是从公开数据中寻找值得进一步验证的生物学问题，而不是继续堆叠方法。五条探索线可以并行推进；某条线出现具体候选后，再针对该候选补充必要验证。

共同分析策略是：同时保留高置信 LAMCORE-like 集、宽松 LAM-like 集和 unrestricted guardrail 集；用连续表达程序、donor-level 汇总、空间位置、扰动反应和外部数据共同判断。聚类或 NMF 产生的是线索，不直接等同于新的细胞亚型。

### 9.1 多细胞 protease–antiprotease 空间生态位

核心问题是：肺囊性破坏附近的 protease 信号，究竟来自一种细胞，还是来自 LAMCORE、LAF/成纤维细胞和免疫/内皮细胞的共同贡献？

首轮分析已经在 GSE302356 的 Visium、Visium HD 和 Xenium 中分别计算空间 source activity、protease、antiprotease 和 `proteolytic_balance_z`。正式分析不预先规定某个 protease 属于哪种细胞，而是先用实际单细胞中的非 protease identity markers 推断 source，再测量各 protease 的表达贡献。

下一步重点是逐个基因回答：哪些状态表达 CTSK、MMP、ELANE、PRTN3、CTSS 等；各状态贡献多少；protease 增强是否伴随 TIMP/SERPIN 等抑制系统不足或空间分离；不同空间技术是否在各自分辨率下指向同一局部 niche。没有可复核的 cyst wall/lesion mask 时，只报告空间邻域关系，不声称囊壁定位。

### 9.2 rapamycin 后仍保留的 ECM/protease 程序

核心问题是：mTOR 抑制是否降低生长，但没有完全消除与 ECM 重塑和蛋白水解相关的病理程序？

GSE179044 使用 TSC2-loss 与 WT 在 vehicle/rapamycin、hydrogel/plastic 条件下的效应量比较，计算 suppression fraction，而不是用“rapamycin 后不显著”定义持续程序。首轮结果中，ELANE 在 hydrogel 和 plastic 两种环境都达到重复方向一致的 partial-retention 条件；MMP2 目前只在 hydrogel 中满足条件。GSE84476 和已登记的 GSE104335 可作为 TSC2/rapamycin/SRPK2 相关的补充扰动证据。

后续把 replicate-stable retained genes 与人类 LAMCORE、GSE190260、GSE217108、GSE302356 以及空间 protease niche 对接。只有同时获得扰动、人体 LAM/空间和至少一种正交支持时，才升级为 `sirolimus-persistent mechanism` 候选；在缺少治疗患者配对数据时，不称为临床耐药。

### 9.3 LAM 细胞进入肺后获得的新能力

核心问题是：哪些程序属于 LAM 本身的转化或谱系特征，哪些是进入肺后获得的能力，哪些只是在局部肺部 niche 中被诱导？

使用 normal uterus、LAM uterus、normal lung 和 pulmonary LAM 四组，计算：

```text
LAM transformation = LAM uterus - normal uterus
pulmonary disease effect = pulmonary LAM - normal lung
lung-acquired interaction = pulmonary disease effect - LAM transformation
```

首轮结果提示 lung-adaptation 是 lung-acquired candidate，而 ECM 更接近 LAM transformation/lineage-associated program。下一步重点检查 ECM attachment、淋巴管相互作用、迁移/侵袭、免疫逃逸、缺氧/氧化应激适应、生存和分泌程序，并将真正的 lung-acquired 候选与空间 niche 和 ATAC 证据连接。由于各组并非同一患者配对，interaction 只作群体层面的比较证据，不称为个体内因果效应。

### 9.4 CORE1/2/3 之外的新程序与患者特异状态

核心问题不是简单寻找“CORE4”，而是判断现有 CORE1/2/3、ECM、mTOR、SLS/IS、MDK、HOX-PBX、LAF-niche、TREM2/TYROBP 等框架是否仍不能充分解释某些重复表达程序。

程序发现以 pooled consensus NMF 和 donor-wise 独立 NMF 为主，再通过 gene-set、方向、TF/regulon、pathway 和 loading similarity 构建 meta-program。已知程序 benchmark 是必要阳性对照；如果连 CORE2 等已知程序都无法稳定找回，就不能把未知程序匹配弱解释为患者异质性。

状态空间还保留 KMeans 和 HDBSCAN 两个互补视角。KMeans 强制分组，HDBSCAN 允许无法形成高密度结构的细胞保留为噪声。当前 HDBSCAN 主设定在 140 个候选中得到 2 个局部密度簇、28 个非噪声细胞，且只涉及 LAM3/LAM4；这更支持继续研究连续状态或局部状态结构，而不是宣布跨 donor 新亚型。结果见 `results/discovery_hdbscan/`。

### 9.5 肺内机制与血浆/EV 蛋白

核心问题是：能否把肺内的细胞状态和局部机制连接到循环中的可检测蛋白？

优先处理 MSV000099051 plasma EV proteomics 和公开 SomaScan 差异蛋白候选，但只把它们作为跨队列证据整合。重点不是生成一份孤立蛋白列表，而是寻找这样的链条：Protein X 在 LAM 血浆/EV 中异常，单细胞提示其来源于某种 LAMCORE/LAF/免疫状态，空间数据支持其位于病灶相关区域，ATAC 或扰动数据支持其调控，并且能提出蛋白检测或功能实验预测。没有同一患者配对时，不报告患者级 RNA–蛋白相关性。

### 9.6 共同的证据升级路径

每条探索线都沿着以下路径推进：

```text
新现象
→ 多个 PatientID 或跨队列重复
→ 已知文献不能完全解释
→ 空间、ATAC、蛋白或扰动等第二类证据
→ 明确可证伪、可实验验证的预测
```

候选结果分为已知程序、已知程序的新扩展、跨 donor 高可信新程序、探索性假说，以及技术或身份不确定结果。细胞通信、TF、regulon 和 pathway 方法只产生候选机制，不把算法预测直接写成真实通信或已证实机制。

## 10. 交付物

- 高置信、宽松和 unrestricted 候选细胞表；
- pooled program、donor-wise program 和 meta-program 对照表；
- known-program evidence scope 表；
- CORE3-like 结构化评分表；
- PatientID/assay 证据矩阵；
- donor-level expression 和 pseudobulk 表；
- 中文和英文研究报告；
- 分级 `LAM Research Hypothesis Cards`；
- data manifest、run manifest、参数和软件版本记录；
- 可从头执行脚本。

当前 GSE190260、GSE217108 和 GSE302356 已转换为可分析的 `.h5ad`，并完成了完整参数的外部 pooled/donor-wise 程序发现。初步跨数据集 top-gene matching 在 Jaccard 0.15 阈值下尚未发现来自不同 PatientID 集合的强 meta-program；达到阈值的匹配主要包含同一 LAM32 donor，因此目前仍不能把任何候选称为已确认的新发现。跨数据集报告见 `results/program_discovery/cross_dataset/`。

## 11. 默认限制

- 不下载 FASTQ；
- 不进行湿实验；
- 公开处理后数据优先；
- 受控访问原始数据不作为主流程依赖；
- 任何尚未获得独立 donor 或实验支持的结果，均只作为研究线索。

## 执行结果摘要（更新至 2026-08-27）

已完成以下实现：

- 空间 source attribution 改为使用实际单细胞中的非 protease identity-marker 表达，protease/antiprotease 在 source labeling 后独立测量；
- proteolytic balance 改为 modality 内标准化 protease activity 减 antiprotease activity，不使用数学比值；
- GSE179044 改用 TSC2-loss effect before/after rapamycin 和 suppression fraction 定义 retained program；
- 肺适应新增 normal uterus、LAM uterus、normal lung、pulmonary LAM 四组 interaction；
- known-program benchmark 新增 donor-level expression 和 NMF loading similarity；
- 新增 retained genes、空间 source、肺适应和蛋白候选的跨模态证据交叉表。

当前结果显示：CORE1/2/3 identity 和 ECM 程序可被 expression/loading benchmark 找回，但 SLS/IS/MDK 恢复不稳定；lung_adaptation 仍是 lung-acquired candidate，ECM 更接近 LAM transformation；ELANE 在 hydrogel 和 plastic 中均达到重复方向一致的 partial-retention 条件，是目前首个跨环境 protease 候选，MMP2 仅在 hydrogel 中满足条件；蛋白 abundance matrix 尚未完成获取。空间 source attribution 已修正为 donor 内跨 source 正确归一化，并完成重跑；三种空间技术仍分别建模。

新增的 HDBSCAN 状态异质性分支复用了核心重实现中的 140 个候选细胞和已有状态评分，没有修改原 KMeans 脚本。主设定得到 2 个局部密度簇、28 个非噪声细胞，非噪声细胞只来自 LAM3/LAM4；当前不支持覆盖四个 donor 的稳定离散亚型，后续应优先把它作为连续状态和局部状态结构线索。
