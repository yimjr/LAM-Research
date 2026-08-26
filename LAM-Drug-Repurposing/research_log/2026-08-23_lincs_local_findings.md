# Local LINCS/CMap pass — 2026-08-23

## Execution status

- GSE92742 和 GSE70138 均完成全量 Level 5 signature 评分。
- 每套数据使用完整 10,174-gene BING space 排名；未在 query genes 内部重排。
- 21 个唯一 query 被实际计算。相同实际 gene set 的 target size 被去重；例如 `hydrogel_specific_residual` 的 100/150 合并记录为同一 query。
- GSE92742 产生 9,483,222 条 signature-query 记录；GSE70138 产生 2,337,006 条。
- WTCS=0 的记录均被显式写成 NCS=0；NCS 无缺失。
- 算法 self-test 通过：synthetic ES、up/down 交换、同号 ES 的 WTCS=0、非负 weighted-correlation 权重和 compound/genetic 分离。

## Cross-release observation

两套数据属于同一 LINCS/L1000 体系的不同 phase/release，因此结果被标记为 `cross-phase/cross-release recurrence`，不称为独立生物学复现。

当前 normalized perturbation-level 记录中：

- `replicated_concordant`: 283；
- `replicated_discordant`: 29；
- `replication_available_but_weak`: 22,410；
- `replication_not_available`: 599,970。

大多数扰动在另一 release 中没有测量，不能把它们解释成复现失败。29 个 discordant 记录才是值得优先检查的 release-sensitive/context-sensitive 候选，但它们也可能由 cell line、dose/time、perturbation modality 或 query-specific context 造成，不能直接称为生物学反转。

## Sanity checks and interpretation

- core mTOR/rapamycin、MTOR、RPTOR 和 RHEB perturbations 已进入单独 sanity panel；其结果用于发现算法方向是否合理，不作为硬性成功条件。
- PI3K/AKT/LAMTOR3 单独作为扩展机制 panel，不因未排在前列而判定流程失败。
- 当前没有把任一 compound 或 target 升级为 Tier 1；LINCS 连接性本身不替代外部模型、人体状态、generic cytotoxicity、target perturbation concordance 和人体暴露过滤。

## Research prompts

1. 29 个真实两-release 冲突的 perturbations 可作为“release-sensitive”审计集合，重点检查其是否集中在单一 cell line、dose/time 或 perturbation class。
2. `environment_dependent_escape` 原始 signature 只有 4 up + 4 down；映射到 BING 后为 0 up + 2 down，因此当前仅能作为 exploratory query，不能产生强 CMap 结论。这是输入 feature-space 限制，不是生物学证据。
3. 如果一个候选只在某一 release 出现，应优先查看 `replication_not_available` 和 context coverage，而不是直接删除；真正需要降级的是两边都测到且稳定方向冲突的 `replicated_discordant`。
