# Cross-Dataset Program Comparison Report

> This report compares top-50 gene overlap across discovery runs; overlap is not treated as biological replication by itself.

## Current result

- Program pairs compared: 73;
- Matches with Jaccard ≥ 0.15: 2;
- Strong matches from different PatientID sets: 0;
- Matches involving the same PatientID: 2, mainly LAM32 across GSE217108 and GSE302356.

No stable meta-program from different PatientID sets reached this threshold. This means that the current candidate definition, feature selection, and top-gene matching did not produce a strong signal; it does not prove that no cross-patient program exists.

## Highest-overlap matches

| pool | dataset | program | dataset | program | Jaccard | PatientID relation |
|---|---|---|---|---|---:|---|
| broad_lam_like | GSE217108 | program_1 | GSE302356 | program_5 | 0.220 | same_patient_overlap_present |
| broad_lam_like | GSE217108 | program_4 | GSE302356 | program_4 | 0.220 | same_patient_overlap_present |
| high_confidence | GSE190260 | program_1 | GSE302356 | program_2 | 0.136 | different_patient_sets |
| high_confidence | GSE217108 | program_3 | GSE302356 | program_1 | 0.111 | same_patient_overlap_present |
| unrestricted_lam | GSE190260 | program_4 | GSE302356 | program_5 | 0.099 | different_patient_sets |
| broad_lam_like | GSE190260 | program_1 | GSE302356 | program_2 | 0.087 | different_patient_sets |
| unrestricted_lam | GSE217108 | program_6 | GSE302356 | program_4 | 0.075 | same_patient_overlap_present |
| broad_lam_like | GSE217108 | program_2 | GSE302356 | program_1 | 0.075 | same_patient_overlap_present |
| high_confidence | GSE217108 | program_5 | GSE302356 | program_5 | 0.064 | same_patient_overlap_present |
| broad_lam_like | GSE135851_core | program_2 | GSE190260 | program_3 | 0.053 | different_patient_sets |

## Next steps

1. Match independently discovered donor programs, not only pooled NMF programs;
2. Test rank-based scores, known-state explained variance, and leave-one-donor-out stability;
3. Use GSE217108 ATAC and GSE302356 ATAC/spatial data for orthogonal support;
4. Report same-patient cross-assay repetition separately from true different-patient replication.
