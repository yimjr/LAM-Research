# Stage 17 — Cross-dataset identity calibration audit

This is a read-only audit of Stage 16. It does not change the formal candidate gate, retrain scVI, or recluster cells.

Positive-reference cells audited: 2,257
Focus dataset: GSE190260
Formal 777-gene signature at audit time: available

## GSE190260 root-cause summary

- Positive references: 2,117
- Stage 16 core recovered: 0 (0.000)
- Stage 16 core+boundary recovered: 1,499 (0.708)
- Missed: 618
- Primary failure category: competing_penalty_only
- Identity-marker dropout (<2 detected): 1.000
- Median final score: -1.007; shift vs other positive references: -6.611
- Median competing-lineage penalty: 4.950
- LODO core+boundary recovery: 0.000

## Interpretation boundary

Failure categories are computed from the Stage 16 component thresholds and exclusion flags, not from manual state inspection. The core3_like label remains an upstream reference label; it is not treated as a new formal candidate assignment. A newly available formal signature is supplemental in this audit and is not inserted into the historical Stage 16 score.

## Outputs

- positive_reference_failures.csv: cell-level reference decomposition and failure reason.
- component_scores_by_dataset.csv: positive-reference score distributions.
- marker_detection_by_dataset.csv: raw-count dropout and depth audit.
- competing_lineage_penalty_audit.csv: lineage-specific penalty summary.
- counterfactual_calibration.csv: raw, within-dataset z-score and percentile score audit.
- lodo_recovery.csv: explicit held-out positive-reference recovery counts.
- root_cause_by_dataset.csv: dataset-level attribution table.
