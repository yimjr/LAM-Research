# GSE277844 translation background and drug concordance

## New observation

The recovery rate of residual-overlap genes was compared with the recovery rate across all selected GSE277844 TSC2-loss translation-abnormal genes. Plastic and hydrogel residual categories were kept separate.

- All 196 selected translation-abnormal genes: RMC-6272 reduced the KO-vs-WT translation distance for 133/196 (67.9%); eFT-508 did so for 147/196 (75.0%).
- Using the same baseline-effect gate used for rescue summaries (absolute baseline translation effect ≥ 0.5), the background was 176 genes: 128/176 (72.7%) for RMC-6272 and 137/176 (77.8%) for eFT-508.
- `hydrogel_residual_q10`: 15/18 eligible genes were reduced by each drug (83.3%); relative to the gated background, this is +10.6 percentage points for RMC-6272 and +5.5 points for eFT-508. The overlap is small and two-sided Fisher tests were not significant, so this is a prioritization signal rather than proof that residual genes are intrinsically more drug-responsive.
- `persistent_residual_hydrogel`: 8/9 and 9/9 eligible genes were reduced, but the denominator is especially small.
- `persistent_residual_plastic`: RMC-6272 reduced 4/5, whereas eFT-508 reduced 2/5. This does not support a uniform residual-specific advantage and may reflect the limited overlap and cross-model noise.

## Gene-level agreement

Within `hydrogel_residual_q10`, 13/18 eligible genes were reduced by both drugs. RMC-6272-only support was FIBIN and HOXC6; eFT-508-only support was APOL1 and TYMS; ZWINT was reduced by neither. The Jaccard overlap of the two reduced-gene sets was 0.765.

Within `persistent_residual_plastic`, GPR27 and RNF182 were supported by both drugs, HOXC6 and NETO1 only by RMC-6272, and ALCAM by neither. This category is too small for a mechanistic comparison.

## Interpretation

The main result is not simply that both translation-targeting perturbations recover many genes: they do so across much of the full translation-abnormal background. Hydrogel residual overlap shows a modest additional enrichment for recovery under the effect-size-gated comparison, while the drug-level agreement suggests a shared core plus a small drug-specific component. These results do not yet establish LAM-specific mechanism because GSE277844 is a human NPC model and the treatment data are not an independent LAM experiment.

## Follow-up questions

- Do the 13 shared hydrogel genes form a coherent translation or proteostasis module rather than a collection of correlated genes?
- Are FIBIN/HOXC6 versus APOL1/TYMS differences supported by treatment effect direction, not only distance reduction?
- Can the shared and drug-specific modules be tested with selective MNK1/2, mTORC1 or genetic perturbations in a LAM-relevant model?
