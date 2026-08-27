# 研究计划：TSC2-loss 状态、残留程序与药物机制

## 研究路线的演变

最初的假设是：sirolimus/rapamycin 只修复 TSC2-loss 的一部分异常，剩余的 persistent residual 可能是联合治疗靶点。因此先在 GSE179044 中比较 plastic 与 hydrogel 的 residual、reversal 和环境交互。

实际分析显示，plastic reversal 与 hydrogel reversal 之间能够稳定、共同印证的方向一致信号很少。这个结果说明“共同 residual reversal”不适合作为候选生成的硬门槛，并不说明两个环境没有生物学联系。后续候选生成因此转向：保留 sirolimus 作为参照和机制背景，直接以 TSC2-loss 的整体 disease signature 进行药物匹配，分别保留 plastic 与 hydrogel 结果。

当前研究分为三条相互衔接的线：

```text
1. residual 机制线：GSE179044 中 sirolimus 修复、残留和环境依赖程序
2. 候选生成线：TSC2-loss plastic/hydrogel signature → 本地 LINCS → 药物靶点聚合
3. translation 验证线：GSE277844 → 分别比较 plastic/hydrogel residual → 翻译靶向处理
```

## 当前核心问题

哪些转录异常是 TSC2 loss 的稳定组成部分，哪些是 sirolimus response 或 environment-dependent residual；直接匹配 TSC2-loss 状态得到的候选药物及其靶点，是否指向可解释的生物学程序；translation abnormality 是否在 residual 中富集，并能否被 RMC-6272/eFT-508 等翻译靶向处理拉回 WT 附近？

## 因子设计

GSE179044 使用 `G=WT/KO`、`R=vehicle/rapamycin`、`E=plastic/hydrogel` 的完整 2×2×2 设计。令 `D(e,r)=KO(e,r)-WT(e,r)`：

- TSC2-loss：`D(e,vehicle)`；
- rapamycin residual：`D(e,rapamycin)`；
- hydrogel-specific residual：`D(hydrogel,rapamycin)-D(plastic,rapamycin)`，即 rapamycin 条件下的 `G×E` simple interaction；
- environment-dependent escape：`[(KO_R-KO_V)-(WT_R-WT_V)]hydrogel - [...]plastic`，即 `G×R×E`；
- rescue/residual 主指标：shrinkage 后 `signed_residual_ratio=d1/d0`，并保留 `|d1|/|d0|` 作为残留幅度。

direction reversal 先作为机制发现类别，不默认进入 CMap 主查询。

## 数据集角色

- GSE179044：核心发现；
- GSE27982：独立 Tsc2×rapamycin 2×2 外部验证，主要验证 TSC2-loss、rescue/persistence 和 genotype-dependent response；
- GSE16944：MMP2、ECM/invasion 等经典 rapamycin-insensitive program 支持；
- GSE84476：STAT3 与 TSC2/rapamycin 机制；
- GSE104335：SRPK2 与 rapamycin 机制比较；
- GSE135851：人体 LAMCORE 验证；
- GSE302356：LAMCORE1/2/3、LAF-seed/LAF-niche、空间和 ECM niche 验证；
- GSE277844：TSC2-loss translation-up/down 程序与 plastic/hydrogel residual 的分别比较，以及 RMC-6272/eFT-508 的探索性恢复验证；
- LINCS/CMap：以 TSC2-loss disease signature 为主生成药物和遗传扰动候选；residual reversal 作为机制解释线索，不作为所有候选的共同硬门槛。

## 人体映射

以 rank-based gene-set enrichment、patient-level pseudobulk enrichment、snATAC gene activity/regulon 和 spatial module localization 为主；基因交集仅作辅助。

## 候选过滤

当前候选生成的基本条件是：药物在 `tsc2_loss_plastic` 或 `tsc2_loss_hydrogel` 中与 TSC2-loss signature 产生可解释的方向信号，并在两个 LINCS release 间具有 concordant 或可解释的 recurrence。之后再检查 drug × gene 程序、靶点证据、target KD/KO 一致性、LAMCORE/LAF 状态表达、generic cytotoxicity、浓度/人体暴露和与 sirolimus 的互补性。

因此，66 个药物是“进入后续研究的候选集合”，不是已经确认的联合治疗药物。最终 Tier 1 仍允许为 0–5 个，不强行凑数。

## 当前阶段状态

TSC2-loss 候选生成、66 个药物的靶点聚合，以及 GSE277844 translation 与 residual 的初步比较已经完成。下一步重点是 generic stress 去卷积、稳定药物模块的人体 LAM 映射、选择性遗传/药理扰动验证，并将结果整理为机制假说；本阶段不直接生成 Tier 1 联合治疗结论。
