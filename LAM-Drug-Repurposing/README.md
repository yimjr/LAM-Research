# LAM Drug Repurposing

本项目围绕 TSC2 loss 后的 LAM-like 异常状态，寻找值得进一步验证的药物作用程序和联合治疗线索。sirolimus/rapamycin residual 是研究的起点，但不是后来候选药物筛选的唯一门槛。

## 研究路线的实际演变

### 起点：研究 sirolimus 之后仍然残留什么

GSE179044 提供了完整的 `TSC2 × rapamycin × environment` 设计，因此最初的主线是：

```text
TSC2 loss
  → sirolimus 修复了什么
  → plastic/hydrogel 中还残留什么
  → 哪些药物能够 reversal 这些 residual
```

这条线仍然用于解释 residual、环境效应和机制，但在比较 plastic reversal 与 hydrogel reversal 时发现，能够在两个环境中以稳定、共同方向得到印证的 reversal 信号很少。两个环境并非完全没有共同生物学，而是缺少足够稳健的共同 reversal signature，不能把“同时被两个环境支持”继续作为候选药物的必要条件。

### 转折：候选生成改为直接匹配 TSC2-loss 状态

因此，候选生成阶段保留 sirolimus 作为背景、参照和机制问题，但不再强行要求药物先 reversal plastic residual 和 hydrogel residual 的共同部分，而是分别使用 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 的 TSC2-loss disease signature 进行本地 LINCS 匹配：

```text
TSC2-loss transcriptional state
  → plastic / hydrogel 分别匹配 LINCS perturbations
  → 两个 LINCS release 的 cross-phase/cross-release recurrence
  → 任一环境中符合条件的 concordant compounds
  → 去重得到 66 个候选药物
  → 药物靶点聚合与机制分析
```

这里的“按 TSC2 变化”不是只看 TSC2 这个单独基因的表达量，而是使用 WT 与 TSC2-null 之间的整体转录差异 signature。当前主候选表包含 258 条候选记录，去除 dataset/query-size 重复后为 66 个唯一药物；原先基于 `tsc2_loss_plastic` 的 92 行/29 药物结果保留为历史子集。

### 当前的后续验证

候选药物先做基因程序和靶点层面的聚合，再开展 translation 方向验证：

```text
66 个 TSC2-loss 匹配药物
  → drug × gene response、reversal/mimic 程序和跨 release 稳定模块
  → ChEMBL / PubChem / BindingDB 靶点证据聚合
  → generic stress/cytotoxicity 与选择性机制区分
  → 人体 LAM 状态和 niche 映射
  → 形成可检验的机制假说

GSE277844 TSC2-loss translation program
  → 分别比较 plastic residual 与 hydrogel residual
  → 检查 RMC-6272/eFT-508 是否使重叠基因的翻译状态向 WT 靠近
  → 对共同恢复的 13 个基因做功能注释
```

GSE277844 这条线是对 TSC2-loss/residual 关系的跨模型机制验证，不是重新生成 66 个候选药物，也不要求 plastic residual 与 hydrogel residual 先取交集。

完整的数据集证据链仍为：

```text
GSE179044
  ├─ residual / environment 机制发现
  └─ TSC2-loss signature 候选生成
       ↓
GSE27982 / GSE16944：外部支持
GSE84476 / GSE104335：STAT3、SRPK2 等机制解释
GSE135851 / GSE302356：人体 LAMCORE/LAF/niche 映射
GSE277844：translation 方向验证
LINCS 本地 Level 5：药物和遗传扰动匹配
```

前期的 GSE179044 因子分解、GSE27982/GSE16944 外部支持、LINCS 候选生成、候选靶点聚合，以及 GSE277844 translation 与 residual 的初步比较已经完成。代码仅保留支持研究复现、效应量计算、模块富集和候选过滤所需的最小功能。

## 研究边界

- 计算结果用于提出实验假说，不等同于治疗建议。
- GSE179044 内部 contrasts 是正交比较，不是独立重复。
- GSE16944 仅作为经典 rapamycin-insensitive program 的外部支持，不计算完整 residual 或 `G×R`。
- GSE27982 的 `G×R` 只称为 genotype-dependent rapamycin response；是否为 escape 必须结合基线、方向、功能和其他模型判断。
- GSE302356 的 RNA、ATAC 和空间样本不自动视作同一患者的配对多组学验证。

## 运行

```bash
./.venv/bin/python scripts/download_geo.py --dataset GSE179044 --dataset GSE27982
./.venv/bin/python scripts/analyze_factorial.py
./.venv/bin/python scripts/analyze_external.py
./.venv/bin/python scripts/analyze_support.py
./.venv/bin/python scripts/analyze_mechanisms.py
Rscript scripts/analyze_gse104335_chp.R
./.venv/bin/python scripts/analyze_gse302356_raw.py
./.venv/bin/python scripts/build_signatures.py
./.venv/bin/python scripts/map_human_states.py
./.venv/bin/python scripts/filter_candidates.py
./.venv/bin/python scripts/analyze_lincs_local.py
./.venv/bin/python scripts/download_resumable.py \\
  https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277844/suppl/GSE277844_raw_counts.txt.gz \\
  data/raw/GSE277844/GSE277844_raw_counts.txt.gz
./.venv/bin/python scripts/analyze_gse277844_translation_residuals.py
./.venv/bin/python scripts/annotate_gse277844_hydrogel_translation_core.py
```

结果写入 `results/`，科学观察和新假说记录在 `research_log/`。翻译程序方向的说明位于 `translation_analysis/`，结构化结果位于 `data/processed/translation_analysis/`，阶段报告位于 `reports/07_translation_analysis/`。

项目结构按研究环节区分：候选排名和筛选结果位于 `results/candidates/`；候选表之后的结构化分析数据位于 `data/processed/candidate_analysis/`；候选分析环节的说明和机制假说位于 `candidate_analysis/`；各环节的稳定阶段报告位于 `reports/`。

## 本地 LINCS/CMap 分析

本项目实际使用已下载的 GSE92742/GSE70138 Level 5 GCTX 与配套 metadata，由 `scripts/analyze_lincs_local.py` 在本地计算 connectivity、WTCS/NCS、上下文汇总和候选排序，不依赖在线查询服务或 API key：

```bash
./.venv/bin/python scripts/analyze_lincs_local.py --datasets GSE92742 GSE70138
```

如果两个 dataset-level Parquet 已经完成，只需重做聚合：

```bash
./.venv/bin/python scripts/analyze_lincs_local.py --aggregate-only --datasets GSE92742 GSE70138
```

两套 LINCS release 的比较标记为 `cross-phase/cross-release recurrence`，不是独立生物学复现；本阶段不生成 Tier 1 候选。

首轮结果见 [RESULTS.md](RESULTS.md)。

各阶段报告入口：[reports/README.md](reports/README.md)。翻译程序方向说明：[translation_analysis/README.md](translation_analysis/README.md)。
