# Statistical evidence appendix

| Evidence item | Value | Source |
|---|---:|---|
| Stage6 high-confidence candidate cells | 5,378 | `results/stage1_6/stage6_checkpoint.json` and grid tables |
| Stage6 LAM-only reference clusters | 12 | `results/stage1_6/stage6_checkpoint.json` |
| Stage6 patient/dataset/assay ARI | 0.152554 / 0.086786 / 0.009094 | `results/stage1_6/stage6_driver_ari.csv` |
| Stage7 consensus states | 20 | `results/stage7/state_consensus_state_summary.csv` |
| Stage10 DE rows | 159,059 | `results/stage10/state_de_results.csv` |
| Stage15 marker-combo candidates | 5,238 | `results/stage15/root_cause_evidence.csv` |
| Stage15 1-UMI-only marker-combo cells | 1,443 | `results/stage15/root_cause_evidence.csv` |
| Stage18 State15 LAMCORE median | 0.5125 | `results/stage18/state15_anchor_report.md` |
| Stage19 LAM1163 composition enrichment | 6.8301 | `results/stage19/state15_patient_composition.csv` |
| Stage21 candidate-only independent slope | -0.015466 | `results/stage21/gradient_models.csv` |
| Stage21 matched-null empirical two-sided p | 0.001996 | `results/stage21/matched_anchor_null.csv` |
| Stage22 State16 independent branch slope | -0.023931 | `results/stage22/branch_evidence_summary.csv` |
| Stage22 State16 corrected matched-null left-tail p | 0.439122 | `results/stage22/branch_evidence_summary.csv` |
| Stage22 State16 corrected left-tail BH q | 1.000000 | `results/stage22/branch_evidence_summary.csv` |
| Stage22 State16 corrected matched-null two-sided p | 0.878244 | `results/stage22/branch_evidence_summary.csv` |
| Stage22 State16 corrected two-sided BH q | 0.878244 | `results/stage22/branch_evidence_summary.csv` |

The Stage21 slope, Spearman rho and binned medians are different estimands/scopes. Stage22 uses the corrected local 1–3-hop scope, distance-structure matched null, direct empirical tails and BH correction; its branch labels are not comparable to the superseded raw-p-value-only version. Their signs and shapes are intentionally reported without forced reconciliation.
