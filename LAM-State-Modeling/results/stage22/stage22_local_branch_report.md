# Stage 22：State 15 局部分支分解

本阶段固定使用 Stage 20/21 的 22,261 个细胞、State 15 的 200 个 anchor、既有 `X_scVI` 和 Stage 21 scores；不重训 scVI、不重新聚类、不修改 candidate gate、不修改 State 1–20 标签。

## Frozen scope

- Main cells: 22261 (5378 candidates + 16883 boundary).
- State 15 anchor: 200 cells; ID SHA-256 `bb060490484a7e99ff50675727b6111dd387db51c9ea6eaa1d962626ffeaa5ef`.
- Local graph: undirected form of the same `X_scVI` k=30 neighbor scope; no Leiden.

## State 15 connectivity and branch selection

- Fixed branch candidates selected: 4.
- Selection rule: external existing state, at least 10 one-hop cells, and at least 2 patients among those one-hop cells; boundary is not promoted to a branch.
详见 `state15_state_connectivity.csv` 和 `branch_candidates.csv`。

## Branch evidence

branch_id source_state  1hop_cells  patient_count  dataset_count  LAMCORE_independent_latent_slope  CORE3_latent_slope  patient_LAMCORE_decrease_fraction  patient_LAMCORE_slope_median  patient_LAMCORE_slope_negative_fraction  patient_LAMCORE_slope_n  LOPO_LAMCORE_slope_negative_fraction  LOPO_LAMCORE_slope_n dominant_competing_lineage  dominant_competing_lineage_slope  matched_null_empirical_left_p  matched_null_empirical_right_p  matched_null_empirical_p  matched_null_q_value  matched_null_two_sided_q_value  matched_null_left_q_value  matched_null_right_q_value           branch_label_pvalue       branch_label_qvalue             evidence_label       raw_evidence_label_before_fdr
branch_01     State_16          96              8              4                         -0.023931           -0.217391                           0.714286                     -0.028166                                 1.000000                        7                              1.000000                    10                       T_NK                          0.056351                       0.439122                        0.562874                  0.878244              0.878244                        0.878244                        1.0                    0.562874 matched_null_empirical_left_p matched_null_left_q_value ordinary_lineage_adjacency          ordinary_lineage_adjacency
branch_02     State_12          22              5              3                         -0.022764           -0.042583                           0.666667                     -0.022906                                 1.000000                        3                              1.000000                     9                endothelial                          0.122994                       0.988024                        0.013972                  0.027944              0.055888                        0.055888                        1.0                    0.027944 matched_null_empirical_left_p matched_null_left_q_value ordinary_lineage_adjacency LAM_to_lineage_transition_candidate
branch_03     State_20          22              6              3                          0.010504            0.004125                           0.625000                      0.000616                                 0.375000                        8                              0.000000                    10                 fibroblast                          0.051935                       1.000000                        0.001996                  0.003992              0.015968                        0.015968                        1.0                    0.007984 matched_null_empirical_left_p matched_null_left_q_value ordinary_lineage_adjacency          ordinary_lineage_adjacency
branch_04      State_7          11              3              3                          0.001486           -0.011434                           0.500000                     -0.003883                                 0.666667                        6                              0.090909                    11                mesothelial                          0.020780                       0.896208                        0.105788                  0.211577              0.282102                        0.282102                        1.0                    0.141051 matched_null_empirical_left_p matched_null_left_q_value ordinary_lineage_adjacency          ordinary_lineage_adjacency

## Boundary

Boundary cells within 1–3 hops are assigned only to a local direction when the direct branch-neighbor count is at least 2 and strictly exceeds the second branch; unresolved cells remain unresolved and no new LAM label is produced.

## Matched null

- Null repetitions per available branch: 500.
- Real and null scopes are identical: non-State15 local 1–3-hop cells; null cells match patient×dataset, cell count and five-bin local distance structure.
- Empirical p-values use direct left/right tails of the observed null distribution; the two-sided p is `2*min(left,right)`, with no zero-centered or symmetry assumption.
- Benjamini–Hochberg q-values are computed separately for left-tail, right-tail and two-sided empirical p-values across all selected branches; LAM/transition labels use the direction-matched left-tail q, while two-sided q remains a general difference statistic.
- Per-patient slopes and leave-one-patient-out fits are recorded in `branch_patient_consistency.csv` and `branch_patient_lopo.csv`.

## Stage 22 checkpoint

- `ordinary_lineage_adjacency_dominates`
- Directly connected branches are better described as ordinary lineage adjacency than LAM-preserving branches.

## Outputs

- all_branch_gradients.csv
- boundary_branch_extension.csv
- boundary_local_branch_assignment.csv
- branch_candidates.csv
- branch_evidence_summary.csv
- branch_gradient_models.csv
- branch_matched_null.csv
- branch_patient_consistency.csv
- branch_patient_lopo.csv
- stage22_local_branch_report.md
- stage22_manifest.json
- state15_local_graph_cells.csv
- state15_state_connectivity.csv
- state16_branch_gradient.csv
- state16_branch_position.csv
- state16_patient_branch_consistency.csv
