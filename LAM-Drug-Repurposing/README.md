# LAM Drug Repurposing

本项目围绕一个具体科学问题展开：

> TSC2 loss 之后，rapamycin 修复了什么、留下了什么、诱导了什么；这些残留程序是否存在于特定的人体 LAMCORE/LAF 状态，并能否导向与 sirolimus 互补的药物候选？

研究主线是：

```text
GSE179044 → GSE27982/GSE16944 → GSE84476/GSE104335
           → GSE135851/GSE302356 → LINCS/CMap
```

首轮分析优先完成 GSE179044 的正式 `TSC2 × rapamycin × environment` 分解，以及 GSE27982 的外部 2×2 验证。代码仅保留支持研究复现、效应量计算、模块富集和候选过滤所需的最小功能。

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
```

结果写入 `results/`，科学观察和新假说记录在 `research_log/`。

项目结构按研究环节区分：候选排名和筛选结果位于 `results/candidates/`；候选表之后的结构化分析数据位于 `data/processed/candidate_analysis/`；候选分析环节的说明和机制假说位于 `candidate_analysis/`；各环节的稳定阶段报告位于 `reports/`。

## CLUE/CMap API 接入

已加入 `scripts/clue_api.py`。它使用 CLUE 的 `user_key` 请求头，先通过 `/genes` 将当前 signature 的 gene symbol 映射为 Entrez ID，再把默认的 8 类 residual/escape signature 作为 L1000 batch query 提交到 `/api/jobs`。direction reversal 默认不提交。

先只生成本地计划（不联网、不需要密钥）：

```bash
./.venv/bin/python scripts/clue_api.py prepare
```

在 CLUE 账户页面取得 API key 后，仅在当前终端设置环境变量，不要写入项目文件：

```bash
export CLUE_API_KEY='在本地终端粘贴你的 key'
./.venv/bin/python scripts/clue_api.py doctor
./.venv/bin/python scripts/clue_api.py submit
```

成功提交后，按返回的 `job_id` 查询并下载：

```bash
./.venv/bin/python scripts/clue_api.py poll JOB_ID
./.venv/bin/python scripts/clue_api.py download JOB_ID
```

映射表、GMT 输入、job 响应写入 `results/cmap/`；归档写入 `data/raw/CLUE/`。连接器不会把密钥写盘，也不会把 CLUE 密钥附带到结果归档下载地址。

首轮结果见 [RESULTS.md](RESULTS.md)。

本地 LINCS/CMap 分析使用已下载的 GSE92742/GSE70138 Level 5 GCTX，不依赖 CLUE 在线服务：

```bash
./.venv/bin/python scripts/analyze_lincs_local.py --datasets GSE92742 GSE70138
```

如果两个 dataset-level Parquet 已经完成，只需重做聚合：

```bash
./.venv/bin/python scripts/analyze_lincs_local.py --aggregate-only --datasets GSE92742 GSE70138
```

两套 LINCS release 的比较标记为 `cross-phase/cross-release recurrence`，不是独立生物学复现；本阶段不生成 Tier 1 候选。

各阶段报告入口：[reports/README.md](reports/README.md)。
