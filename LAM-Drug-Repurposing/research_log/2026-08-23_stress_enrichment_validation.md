# 2026-08-23 — stress deconvolution, enrichment, and independent target references

## Completed

- Applied operational local generic apoptosis, UPR/ER-stress, heat-shock,
  cell-cycle-arrest and proteostasis/translation-stress templates to the 204
  common cross-release genes.
- Produced raw and stress-adjusted gene effects, generic-stress burden,
  module retention, and leading-gene enrichment tables.
- Ran GO Biological Process, Reactome and MSigDB Hallmark ORA with the actual
  204-gene common analyzable background and BH correction.
- Added an external RET-oriented profile from GSE49414 and an external
  BTK/ibrutinib resistance-contrast profile from GSE207322.
- Added explicit `not_available` status for pharmacological references that
  were not measured in a LINCS release.

## Main observations

1. The mimic module is not wholly generic stress: mean directional retention
   remains about 0.95 in GSE92742 and 0.98 in GSE70138 after local adjustment.
   However, clofarabine, proteasome, microtubule and related mimic
   perturbations carry higher generic-stress burden and must remain down-weighted
   until curated stress references are added.
2. `reversal_only__stable_2` is the strongest module: Jaccard 0.717, high
   cross-release direction fractions, and retention after stress adjustment.
   Its leading biology is glycolysis/hypoxia plus ECM organization/remodeling.
3. `reversal_only__stable_3` is small and unstable; it remains an exploratory
   transcription/EMT/state-switch lead.
4. No GO/Reactome/MSigDB term is BH-FDR significant at 0.05 for the current
   small modules. Leading genes, direction and cross-release recurrence are
   retained rather than threshold-chasing pathway labels.
5. The RET reference reproduces the expected direction in 46/83 genes of
   stable_2. The BTK reference gives 51/86 in WT TMD8 cells, with less
   consistent responses in C481 mutants. These results support pharmacological
   relevance but not single-target causality.

## Data limitations

- GSE49414 is a TPC1 thyroid-cell RPI-1 versus DMSO experiment; RPI-1 is not a
  RET-only perturbation and is not a LAM model.
- GSE207322 is TMD8 lymphoma with ibrutinib and BTK C481 mutants; it is an
  independent BTK-oriented reference, not a LAM model.
- GSE135851 mapping uses preliminary candidate/other labels; GSE302356 output
  is sample/modality-level state-marker scoring rather than formal paired
  LAMCORE/LAF annotation.
- NTRK3 and MKNK1/2 still lack comparably strong independent perturbation
  datasets in the current local analysis; LINCS genetic/pharmacological
  references remain hypothesis-generating only.

## Reproducibility correction

The downstream loader was corrected to derive its background from the 300-gene
panel mapping audit. It now uses the 204 genes available in both LINCS
releases, rather than all gene-level rows present in the broader comparison
table. Stress, enrichment and human-mapping outputs were rerun after this
correction.

## Next decisive analysis

Add curated stress signatures, then repeat H2 gene-module testing with
selective RET/BTK perturbations and independent genetic perturbation. Map the
resulting module to formal LAMCORE1/2/3 and LAF-seed/LAF-niche labels when the
annotated human matrices are staged.
