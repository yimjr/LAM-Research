# LAMCORE Unknown-State Program Discovery Report: GSE217108

> This report documents program-discovery candidates in GSE217108. It is not a completed independent-validation claim.

Run mode: `configured_full`. The configured full rank, seed, and donor-wise analysis parameters were used.

## Completed in this stage

- Input: 12396 cells and 36601 genes;
- Pools: high-confidence 1075, broad 3855, unrestricted guardrail 12396;
- Pooled NMF and donor-wise independent NMF were both run, followed by preliminary meta-program matching;
- Known programs were compared post hoc and were not regressed out before primary discovery;
- CORE3 was modeled as identity, depth-adjusted low activity, and translation enrichment.

## Interpretation

All programs remain candidates. External GSE190260, GSE217108, and GSE302356 are now available as analysis-ready AnnData and have completed the full configured run; orthogonal modality validation remains. Even after the full RNA analysis, orthogonal modalities and PatientID independence checks remain necessary.

## Run summary

- `high_confidence`: selected rank=6, seed=0, factor stability=1.000.
- `broad_lam_like`: selected rank=5, seed=2, factor stability=1.000.
- `unrestricted_lam`: selected rank=6, seed=2, factor stability=1.000.

## Next steps

1. Download and inspect the public processed matrices;
2. Perform PatientID-aware donor-wise discovery and meta-program matching;
3. Test cross-donor candidates against doublet, assay, depth, known-program, and leave-one-donor-out sensitivities;
4. Use GSE217108 ATAC, GSE302356 ATAC/spatial, and protein data for orthogonal evidence;
5. Upgrade a candidate only when identity, independent donors, and orthogonal evidence support it together.
