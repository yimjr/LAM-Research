# LAMCORE 未知状态程序发现计划

## 1. 目标与结果定位

本项目的核心问题是：

> 在不预先排除弱 marker、低活性或非典型 LAMCORE-like 细胞的前提下，寻找跨患者重复、且不能被现有 LAM 状态框架充分解释的表达程序。

工作顺序为：

> 核心复现完成 → 稳健性验证与新发现并行 → 针对具体候选发现按需加强验证

此前 `results/independent_reanalysis/` 和旧版 `results/discovery/` 的结果继续保留为同一数据集上的独立再分析，不覆盖，也不直接称为原论文复现。未知程序发现阶段的任何候选，在获得独立 PatientID 和正交证据前都只称为探索性结果或研究线索。

## 当前执行状态（2026-08-23）

核心复现基线已冻结，五条研究线已并行进入首轮分析。已完成已知程序阳性对照 benchmark、baseline/strict-QC 固定敏感性门、GSE179044/GSE84476 扰动数据获取、GSE302356 三模态空间首轮分析和肺适应跨组织参考比较。当前最值得继续追踪的是 LAMCORE–protease 空间邻域关系，以及 rapamycin 后 ECM/protease 程序保留的候选机制；两者仍是探索性结果，不是已证实的新机制。

阳性 benchmark 显示当前 top-gene matching 不能稳定找回所有已知程序，尤其 CORE2，因此方向一暂不能把未知程序的弱跨 PatientID 匹配解释为患者特异性。严格 QC 将候选数从 140 降至 85，并显著改变 LAM2/LAM4 的候选分布；后续每张 Hypothesis Card 都必须通过 baseline/strict-QC 检查。详细状态见 `results/report/parallel_research_status_zh.md` 和英文版。

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

## 9. 并行研究方向

- 五条研究线在核心基线建立后并行推进，不要求某一条方向先完成后才能启动其他方向。每条线按照当前证据等级报告，避免把“尚无患者级证据”误写成“不能分析”。
- LAMCORE 内部异质性：寻找 CORE1/2/3、ECM、mTOR、激素、缺氧、增殖、SLS/IS、MDK、HOX-PBX、LAF-niche、TREM2/TYROBP 等已知框架无法充分解释的程序。
- sirolimus-persistent program：可先整合人类 LAM 状态、SLS/IS/MDK、GSE179044 的 TSC2-loss × ECM × rapamycin、GSE84476 的 TSC2/STAT3 干预，以及 GSE104335 的 LAM 621-101 细胞系 rapamycin/SRPK2 微阵列数据，寻找候选的治疗后持续机制；在缺少多个治疗 donor 或配对数据时，不称为患者级治疗耐受，只报告为候选机制。
- 肺囊性破坏：结合 LAMCORE、成纤维细胞、淋巴管、免疫细胞和空间数据，检验 protease/ECM niche 假说。
- 肺适应程序：区分 lineage/origin、lung-adaptation 和 niche-induced program。
- 血浆/EV 蛋白来源：可使用非配对的公开血浆/EV 蛋白组开展跨队列证据整合；只有在同一患者配对时，才进一步报告患者级 RNA–蛋白相关性。

五条线的共同升级路径为：新现象 → 多个 PatientID 或跨队列重复 → 已知文献无法完全解释 → ATAC、空间、蛋白或扰动等第二类证据 → 明确可实验验证的预测。方向 1 重点校准 gene-set matching，并补充 regulon、TF、pathway 和 loading similarity；方向 2–5 可在现有数据基础上直接生成候选 Hypothesis Cards。最终由证据链最完整的方向升级为重点主线，而不是预先规定固定的串行顺序。

细胞通信、TF、regulon 和 pathway 方法只产生候选机制，不把算法预测直接写成真实通信或已证实机制。

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

## 修订版执行结果（2026-08-23）

已完成以下实现：

- 空间 source attribution 改为使用实际单细胞中的非 protease identity-marker 表达，protease/antiprotease 在 source labeling 后独立测量；
- proteolytic balance 改为 modality 内标准化 protease activity 减 antiprotease activity，不使用数学比值；
- GSE179044 改用 TSC2-loss effect before/after rapamycin 和 suppression fraction 定义 retained program；
- 肺适应新增 normal uterus、LAM uterus、normal lung、pulmonary LAM 四组 interaction；
- known-program benchmark 新增 donor-level expression 和 NMF loading similarity；
- 新增 retained genes、空间 source、肺适应和蛋白候选的跨模态证据交叉表。

当前结果显示：CORE1/2/3 identity 和 ECM 程序可被 expression/loading benchmark 找回，但 SLS/IS/MDK 恢复不稳定；lung_adaptation 仍是 lung-acquired candidate，ECM 更接近 LAM transformation；ELANE 在 hydrogel 和 plastic 中均达到重复方向一致的 partial-retention 条件，是目前首个跨环境 protease 候选，MMP2 仅在 hydrogel 中满足条件；蛋白 abundance matrix 尚未完成获取。空间 source attribution 已修正为 donor 内跨 source 正确归一化，并完成重跑；三种空间技术仍分别建模。
