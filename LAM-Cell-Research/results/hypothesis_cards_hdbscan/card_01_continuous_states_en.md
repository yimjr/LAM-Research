# LAMCORE May Contain Continuous Expression States

## Class: exploratory methodological result (same dataset; no independent LAM-donor validation yet)

## Observation
Continuous state programs in the author-style marker candidates showed directional signals across multiple LAM donors: contractile, ecm_remodeling, hormone_related, stress_hypoxia, proliferative. See `lamcore_state_programs_by_donor.csv` for the donor-level differences.

## HDBSCAN result
Under the primary setting, HDBSCAN identified 2 density clusters among 140 candidates, assigned 28 cells, and marked 112 as noise across 2 donors. The primary non-noise cells came only from LAM3 and LAM4; no LAM1 or LAM2 candidate was assigned to a density cluster. Thus HDBSCAN did not support a stable discrete state spanning all four donors. A noise label means that a point was not assigned to a sufficiently dense region in this state space, not that the cell lacks biological meaning.

## Donors, cells and pathways
The unit is donor across LAM1–LAM4; programs include contractile, ECM remodeling, stress/hypoxia, inflammatory, hormone-related, metabolic and mTOR-related states.

## Robustness
Candidates were defined by known markers plus the author-style graph. Phase 2 compared doublets, QC, clustering seeds/resolutions, the 777-gene module score, a rank-based score, assay strata and leave-one-donor-out analyses. Because discovery and assessment use the same cohort, this is not independent validation.

## Alternatives
The candidate definition overlaps some contractile/ECM markers; tissue processing, assay, cell cycle, stress and donor biology may contribute.

## Next validation
Test the programs in independent LAM donors, spatial transcriptomics, protein or snATAC data within the same marker-defined LAMCORE-like cells, then assess TSC2/mTOR, ECM and lymphatic mechanisms experimentally.

## Novelty / confidence / priority
Novelty: medium; confidence: low-to-medium; priority: medium. The result supports further study of continuous states and local dense structure, not a new cross-donor subtype claim.
