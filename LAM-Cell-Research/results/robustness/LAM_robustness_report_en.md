# Targeted Robustness Validation for LAMCORE

## Purpose

Phase 2 is a quality filter for concrete candidate findings rather than the endpoint of a methods-comparison project. With the Phase 1 candidate labels fixed, we compare doublet removal, loose/strict QC, clustering seeds/resolutions, the 777-gene module score, a rank-based score, assay strata and leave-one-donor-out analyses.

- Phase 1 without doublet removal: 30708 cells and 140 candidates.
- After removing predicted doublets: 30077 cells and 140 candidates.
- Predicted doublets: 631.
- Loose/strict QC candidate counts: 140 / 85.

These checks assess technical sensitivity; they are not independent donor validation and do not use a mechanical 3/4-donor pass rule.
