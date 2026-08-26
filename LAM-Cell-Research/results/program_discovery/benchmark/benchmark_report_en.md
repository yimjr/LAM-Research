# Known-program positive-control benchmark (first pass; superseded by v2)

> This file is retained for historical traceability. Use [`benchmark_report_v2_en.md`](benchmark_report_v2_en.md) and `known_program_extended_benchmark_summary.csv` for the current conclusion.

## Purpose

Before interpreting weak cross-PatientID matching of an unknown program, we test whether the current top-gene/Jaccard representation can recover programs that are already expected to be present. If known programs also fail to match, weak matching of an unknown program is not interpretable as patient-specific biology.

## First-pass result

At top-50, using the current approximate gene sets and the union of output-program genes and known-program genes as the null universe, recovery was uneven:

- ECM_remodeling: pooled detectable fraction 0.667;
- SLS_stem_like: 0.444;
- CORE1: 0.333;
- IS_inflammatory: 0.333;
- LAM_myogenic_contractile / contractile: 0.111;
- MDK_dormancy_persistence: 0.111;
- CORE2: 0.000.

This is not evidence that CORE2 is absent. It shows that the current top-gene matching is sensitive to gene-set definition, candidate pool, module number, and patient/assay variation. Because CORE2 is not reliably recovered, a weak cross-PatientID match for an unknown program cannot currently be called patient-specific novel biology.

## Historical conclusion

In this first-pass report, Direction 1 remained in calibration. The v2 pass added donor-level expression scores and NMF loading similarity; v2 shows that CORE2 can be recovered at the expression and loading levels, so the old “CORE2=0” statement must not be used as the current conclusion. Regulon/TF and pathway layers remain pending; pooled and donor-wise discovery paths are retained.

## Reproducible files

- `known_program_matching_benchmark.csv`
- `known_program_matching_summary.csv`
- `benchmark_manifest.json`

Note: the current null universe is the union of generated program genes and known-program genes, not the complete AnnData gene universe. This is a conservative calibration result, not a formal inferential test.
