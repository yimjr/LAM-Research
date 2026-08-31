# State 15 LAM-core reference anchor validation

State 15 was frozen from the existing consensus annotation. This stage did not recluster, retrain scVI, or modify any candidate gate.

- Frozen State 15 cells: 200
- Frozen State 15 ID SHA-256: bb060490484a7e99ff50675727b6111dd387db51c9ea6eaa1d962626ffeaa5ef
- State 15 dataset coverage: 4 datasets
- State 15 patient coverage: 7 patients
- Formal signature: available (777 genes)

## Formal LAMCORE and author-label evidence

- State 15 LAMCORE median: 0.5125
- Boundary LAMCORE median: 0.0936
- Normal/control LAMCORE median: 0.0513
- Overall author-style enrichment fold: 13.9408
- Overall Fisher exact one-sided p-value: 2.229e-36

## Anchor decision

- Decision: `provisional_reference_candidate_not_formally_upgraded`
- State 15 is present in 7/12 consensus patients; the largest contribution is LAM1163 (127/200 cells).
- Author-style support is present in datasets: GSE135851.
- Formal LAMCORE/comparator separation is supportive, but the current patient and author-label concentration is not sufficient to promote State 15 to a formally cross-patient reference anchor.
- State 15 remains frozen as a provisional reference candidate for a later, explicitly designed expansion analysis; no gate or existing state artifact is changed here.

## Interpretation

The tables retain State 15 as a fixed validation object. Formal LAMCORE elevation, author-label enrichment, patient-level consistency, comparator separation, and latent-neighborhood continuity must be considered together; no single score or cell count is promoted to an automatic gate.

## Outputs

- state15_lamcore_summary.csv
- state15_marker_profile.csv
- state15_author_enrichment.csv
- state15_vs_comparators.csv
- state15_patient_pseudobulk.csv
- state15_patient_consistency.csv
- state15_latent_neighbors.csv
- state15_latent_neighbor_edges.csv
- state15_latent_distance_by_cell.csv
- state15_latent_distance_gradient.csv
