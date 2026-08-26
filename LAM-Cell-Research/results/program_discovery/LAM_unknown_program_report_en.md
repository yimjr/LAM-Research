# LAMCORE Unknown-State Program Discovery Report: GSE135851_core_reproduction

> This report documents program-discovery candidates in GSE135851_core_reproduction. It is not a completed independent-validation claim.

Run mode: `fast_smoke`. Parameters are for pipeline acceptance and not final rank or stability claims.

## Completed in this stage

- Input: 30708 cells and 63677 genes;
- Pools: high-confidence 535, broad 3940, unrestricted guardrail 23228;
- Pooled NMF and donor-wise independent NMF were both run, followed by preliminary meta-program matching;
- Known programs were compared post hoc and were not regressed out before primary discovery;
- CORE3 was modeled as identity, depth-adjusted low activity, and translation enrichment.

## Interpretation

All programs remain candidates. External GSE190260, GSE217108, and GSE302356 are now available as analysis-ready AnnData; this report is a fast smoke run and does not replace the full parameterized analysis. Even with external data available, a smoke run alone cannot establish the final independent-validation tier.

## Run summary

- `broad_lam_like`: selected rank=3, seed=1, factor stability=1.000.
- `high_confidence`: selected rank=3, seed=0, factor stability=1.000.
- `unrestricted_lam`: selected rank=3, seed=0, factor stability=1.000.

## Next steps

1. Download and inspect the public processed matrices;
2. Perform PatientID-aware donor-wise discovery and meta-program matching;
3. Test cross-donor candidates against doublet, assay, depth, known-program, and leave-one-donor-out sensitivities;
4. Use GSE217108 ATAC, GSE302356 ATAC/spatial, and protein data for orthogonal evidence;
5. Upgrade a candidate only when identity, independent donors, and orthogonal evidence support it together.
