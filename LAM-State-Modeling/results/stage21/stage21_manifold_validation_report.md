# Stage 21：State 15-centered manifold validation

本阶段将 State 15 的 200 个细胞严格作为 reference anchor；所有主要 gradient 检验均在其余 22,061 个细胞上完成。复用 Stage 20 的 `distance_to_state15` 和既有 `X_scVI`，不重训 scVI、不重新 Leiden/consensus、不修改 candidate gate。

## Frozen input and independence audit

- Anchor: 200 State 15 cells; ID SHA-256 `bb060490484a7e99ff50675727b6111dd387db51c9ea6eaa1d962626ffeaa5ef`.
- Validation object: 22061 non-State15 cells; candidate-only null pool: 5178 cells.
- LAMCORE formal genes: 777; outside-scVI genes: 557; independent genes: 554.
详见 `lamcore_independence_gene_audit.csv`；Stage 21 的主要独立验证分数为 `LAMCORE_outside_scVI` 和 `LAMCORE_independent`。

## Anchor-excluded distance gradient

- Anchor-excluded binned gradient shows near-vs-far independent LAMCORE decrease: `False`.
- Candidate-only patient-adjusted independent slope: `-0.0154663`.
- 四套 score 的 median/IQR 和完整距离模型分别见 `non_state15_distance_gradient.csv` 与 `gradient_models.csv`。

## Dataset and patient validation

- Dataset-level results: 4 dataset rows for independent LAMCORE.
- Patient-level results: 12 patient rows for independent LAMCORE.
- Stage 20 中非负的患者在 Stage 21 中保留了 distance range、cell count 和 slope/rho 字段，不被直接标记为失败。

## Composition-matched fake-anchor null

- Null repetitions: 500; empirical two-sided p=0.001996; null median slope=0.0037587.
- 每次假 anchor 均按 State 15 的 patient×dataset 组成，从非-State15 candidate pool 抽取；真实比较对象为相同 candidate-only scope。

## State 16 and boundary

- `state16_distance_gradient.csv` 分别给出 pooled 和 patient-stratified 的 near/mid/far profile；`boundary_independent_gradient.csv` 只做 evidence ranking，不产生新 candidate 标签。
- `manifold_connectivity.csv` 使用与 Stage 20 相同的 candidate+boundary `X_scVI` k=30 scope 重建邻接汇总，用于识别单轴或多分支连接；不产生 cluster。

## Stage 21 checkpoint

- `state15_lam_rich_gradient_but_not_robust_manifold`
- Some independent LAM evidence remains outside State 15, but patient/dataset or matched-null evidence is insufficient for a robust manifold claim.

## Outputs

- boundary_independent_gradient.csv
- dataset_independent_gradient.csv
- distance_score_smooth.csv
- gradient_models.csv
- independent_lamcore_scores.csv
- lamcore_independence_gene_audit.csv
- lineage_gradient_by_distance.csv
- manifold_connectivity.csv
- matched_anchor_null.csv
- non_state15_distance_gradient.csv
- patient_independent_gradient.csv
- stage21_manifest.json
- stage21_manifold_validation_report.md
- state16_distance_gradient.csv
