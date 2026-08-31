# Limitations

- State15 仅覆盖 7/12 patients，且 LAM1163 composition enrichment=6.8301。
- External author-style annotation is `not_assayed` for GSE190260/GSE217108/GSE302356.
- Upstream `cell_type` is unknown for all consensus cells; Stage24 analogues are inferred.
- Candidate gate is high-recall and includes ordinary lineage states.
- Dataset heterogeneity, marker dropout and GSE190260 score shift limit cross-dataset calibration.
- No time, spatial, prospective cohort or experimental validation was included.
- Pathway enrichment and regulon outputs are unavailable placeholders.
- Stage21 pooled and patient-adjusted gradient estimands differ; corrected Stage22 branch p-values remain exploratory and State16 no longer passes the corrected matched-null/FDR evidence label.
- Stage22 branch selection and null analysis are local 1–3-hop analyses; 11648 local boundary cells were projected, while farther boundary cells were intentionally not assigned.
