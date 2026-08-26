# LAM Research Hypothesis Card: unrestricted_lam / program_2

## Current tier

Exploratory hypothesis only. It comes from the GSE135851 cohort and has not met the independent-validation standard.

## Observation

`program_2` emerged in pooled NMF of the `unrestricted_lam` candidate pool. Leading genes are: MALAT1, CLDN5, IFI27, PECAM1, MT1M, SPTBN1, IFITM3, CAVIN2, RGCC, MT1E, HYAL2, TCF4.

## Donor reproducibility

Donors meeting the current top-gene matching rule after donor-wise discovery: . This is independent donor discovery evidence, not passive scoring of a pooled model; an empty value means that cross-donor support is currently insufficient.

## Relation to known frameworks

Strongest known-program match: CORE1 (LAMCORE-specific; overlap=1). Partial overlap does not automatically invalidate a candidate; test whether it is a new LAMCORE implementation of a known state, has a LAM-specific component, or adds a new TF/regulon.

## External evidence

External AnnData currently available: GSE190260, GSE217108, GSE302356. ATAC, spatial, and protein evidence have not yet been used for the claim.

## Alternative explanations

- donor, assay, or batch-specific signal;
- cell cycle, doublet, sequencing depth, or low-quality effects;
- a projection of an existing CORE/SLS/IS/ECM state;
- unresolved identity in the broad candidate pool.

## Next validation

Re-discover the program by PatientID in GSE190260, GSE217108, and GSE302356; compare known-program explained variance; then test ATAC, spatial, or protein support. Do not name a new LAMCORE subtype or mechanism before these checks.

## Novelty / confidence / priority

- Novelty: not assessed;
- Current confidence: low to moderate, same-cohort candidate only;
- Priority: medium, conditional on independent donor re-discovery.
