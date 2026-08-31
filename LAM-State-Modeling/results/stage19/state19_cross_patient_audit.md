# Stage 19：State 15 跨患者 identity calibration audit

本阶段固定使用现有 5,378 个 high-confidence candidate 的 consensus 标签，不重新聚类、不重训 scVI、不修改 candidate gate。

## 1. Patient composition baseline

- State 15: 200 cells; 7 patients.
- LAM1163 在全部 candidate pool 中占 0.0930，在 State 15 中占 0.6350，composition enrichment = 6.8301.

```text
patient_id  candidate_pool_total_cells  candidate_pool_fraction  state15_cells  state15_fraction  enrichment  candidate_pool_dataset_count  state15_dataset_count
      LAM1                          75                 0.013946             25             0.125    8.963333                             1                      1
   LAM1158                         182                 0.033842             10             0.050    1.477473                             1                      1
   LAM1163                         500                 0.092971            127             0.635    6.830060                             1                      1
   LAM1164                         745                 0.138527              1             0.005    0.036094                             1                      1
      LAM3                        1678                 0.312012             26             0.130    0.416651                             2                      2
     LAM32                         812                 0.150985             10             0.050    0.331158                             2                      2
     LAM50                         725                 0.134808              1             0.005    0.037090                             1                      1
```

## 2. Author-style annotation availability

只有存在真实逐细胞作者阳性标签的数据集才进入 author enrichment；其他数据集的零值被标记为 `not_assayed`，不解释为 author-negative。

```text
  dataset                                                                                                                     candidate_file  author_style_field_present author_style_annotation_status                               availability_basis  author_style_positive_total  state15_cells  state15_author_style_positive_observed
GSE135851                    /mnt/e/lam-research/LAM-State-Modeling/../LAM-Cell-Research/results/program_discovery/candidate_pool_labels.csv                        True                      available                   positive_author_labels_present                        140.0             50                                    49.0
GSE190260 /mnt/e/lam-research/LAM-State-Modeling/../LAM-Cell-Research/results/program_discovery/external_GSE190260/candidate_pool_labels.csv                        True                    not_assayed external_field_initialized_false_or_file_missing                          NaN            138                                     NaN
GSE217108 /mnt/e/lam-research/LAM-State-Modeling/../LAM-Cell-Research/results/program_discovery/external_GSE217108/candidate_pool_labels.csv                        True                    not_assayed external_field_initialized_false_or_file_missing                          NaN              1                                     NaN
GSE302356 /mnt/e/lam-research/LAM-State-Modeling/../LAM-Cell-Research/results/program_discovery/external_GSE302356/candidate_pool_labels.csv                        True                    not_assayed external_field_initialized_false_or_file_missing                          NaN             11                                     NaN
```

```text
  dataset annotation_status  state15_cells  state15_author_style  other_candidate_cells  other_author_style  state15_author_fraction  other_author_fraction  enrichment_fold  fisher_odds_ratio  fisher_pvalue_greater
GSE135851         available             50                    49                    485                  91                     0.98               0.187629         5.223077         212.153846           8.027702e-31
```

## 3. Patient-matched and leave-one-patient-out evidence

- Patient-matched LAMCORE positive delta fraction: 1.0000.
- LOPO profile reference closer fraction: 0.7714.
- LOPO latent reference closer fraction: 0.9714.

## 4. Removing LAM1163

- Remaining State 15 cells: 73.
- Remaining State 15 LAMCORE median: 0.4683.
- LAMCORE remains above every requested comparator: True.
- Interpretation: 去除 LAM1163 后 State 15 的 LAMCORE 仍高于全部指定 comparator，且 patient-matched LAMCORE 方向保留。

## 5. Stage 19 conclusion

- Classification: `B_patient_enriched_but_biological_profile_persists_after_removal`.
- The patient-composition baseline rules out interpreting the 63.5% LAM1163 fraction as ordinary sampling alone: its State 15 composition enrichment is 6.8301-fold.
- The LAM1163 removal sensitivity, patient-matched direction, and LOPO latent comparison support a real but patient-heterogeneous State 15 biological profile.
- Author-style evidence is available only for GSE135851 in the inherited upstream files; GSE190260, GSE217108 and GSE302356 are `not_assayed`, not author-negative.

## 6. Scope

本阶段仅判断当前冻结 State 15 的跨患者组成基线和生物学复现，不把结果自动写回 candidate gate、consensus clustering 或 atlas。

## Outputs

- author_annotation_availability.csv
- stage19_manifest.json
- state15_author_enrichment_assayed.csv
- state15_lopo_validation.csv
- state15_patient_composition.csv
- state15_patient_matched_comparison.csv
- state15_patient_profiles.csv
- state15_patient_pseudobulk_profiles.csv
- state15_without_LAM1163.csv
- state19_cross_patient_audit.md
