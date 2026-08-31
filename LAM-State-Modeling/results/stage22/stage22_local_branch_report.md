# Stage 22：State 15 局部分支分解

本阶段固定使用 Stage 20/21 的 22,261 个细胞、State 15 的 200 个 anchor、既有 `X_scVI` 和 Stage 21 scores；不重训 scVI、不重新聚类、不修改 candidate gate、不修改 State 1–20 标签。

## Frozen scope

- Main cells: 22261 (5378 candidates + 16883 boundary).
- State 15 anchor: 200 cells; ID SHA-256 `bb060490484a7e99ff50675727b6111dd387db51c9ea6eaa1d962626ffeaa5ef`.
- Local graph: undirected form of the same `X_scVI` k=30 neighbor scope; no Leiden.

## State 15 connectivity and branch selection

- Fixed branch candidates selected: 4.
- Selection rule: external existing state, at least 10 one-hop cells, and at least 2 patients; boundary is not promoted to a branch.
详见 `state15_state_connectivity.csv` 和 `branch_candidates.csv`。

## Branch evidence

branch_id source_state  1hop_cells  patient_count  dataset_count  LAMCORE_independent_latent_slope  CORE3_latent_slope  patient_LAMCORE_decrease_fraction dominant_competing_lineage  dominant_competing_lineage_slope  matched_null_empirical_p                      evidence_label
branch_01     State_16          96             10              4                         -0.023659           -0.216772                           0.714286                       T_NK                          0.050436                  0.025948 LAM_to_lineage_transition_candidate
branch_02     State_12          22              9              4                         -0.024933           -0.040025                           1.000000                endothelial                          0.131434                  0.674651          ordinary_lineage_adjacency
branch_03     State_20          22             10              4                          0.010979            0.002151                           0.625000                 fibroblast                          0.042601                  0.315369          ordinary_lineage_adjacency
branch_04      State_7          11             11              4                          0.002490            0.006365                           0.666667                mesothelial                          0.019902                  1.000000          ordinary_lineage_adjacency

## Boundary

Boundary cells within 1–3 hops are assigned only to a local direction when the direct branch-neighbor count is at least 2 and strictly exceeds the second branch; unresolved cells remain unresolved and no new LAM label is produced.

## Matched null

- Null repetitions per available branch: 500.
- Null scope: non-State15, non-branch cells within 1–3 hops, matched on patient×dataset and branch cell count.
- Per-branch null summaries are recorded in `branch_matched_null.csv` and the manifest.

## Stage 22 checkpoint

- `local_branched_lam_manifold_candidate`
- At least one local branch retains independent LAM identity while moving toward a distinct lineage direction; interpret as a local branched candidate, not a new state.

## Outputs

- all_branch_gradients.csv
- boundary_branch_extension.csv
- boundary_local_branch_assignment.csv
- branch_candidates.csv
- branch_evidence_summary.csv
- branch_gradient_models.csv
- branch_matched_null.csv
- branch_patient_consistency.csv
- stage22_local_branch_report.md
- stage22_manifest.json
- state15_local_graph_cells.csv
- state15_state_connectivity.csv
- state16_branch_gradient.csv
- state16_branch_position.csv
- state16_patient_branch_consistency.csv
