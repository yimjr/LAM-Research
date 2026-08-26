# 研究计划：LAM 残留程序、rapamycin response 与联合治疗候选

## 核心问题

TSC2 loss 后，rapamycin 是否只修复 mTORC1 依赖的一部分异常，而留下 persistent residual、worsened residual 或 genotype-dependent response？这些程序是否定位于人体 LAMCORE1/2/3 或 LAF-seed/LAF-niche，并能否导向联合治疗候选？

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
- LINCS/CMap：药物候选生成，不作为最终证据。

## 人体映射

以 rank-based gene-set enrichment、patient-level pseudobulk enrichment、snATAC gene activity/regulon 和 spatial module localization 为主；基因交集仅作辅助。

## 候选过滤

候选需要通过跨 contrast 重复、外部模型支持、target KD/KO 一致性、LAMCORE/LAF 状态表达、generic cytotoxicity 过滤、浓度/人体暴露和 sirolimus 互补性检查。最终允许 0–5 个 Tier 1 候选，不强行凑数。
