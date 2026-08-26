# Hypothesis Card 04 (revised): multicellular protease–antiprotease spatial niche

## Observation

Visium, Visium HD and Xenium show directionally concordant spatial associations between LAMCORE-like signal and protease signal. Source attribution now starts from actual single-cell expression and measures protease/antiprotease genes independently.

## Important correction

CTSK, MMP, ELANE and CTSS are not assigned to LAMCORE, LAF or immune cells in advance. The first single-cell attribution pass shows gene-specific source patterns and a substantial unclassified component, so no cell population is yet claimed as a definitive source.

## Balance definition

```text
proteolytic_balance_z = standardized protease activity - standardized antiprotease activity
```

No protease/antiprotease ratio is used. The three spatial technologies remain separate; raw units, scores and p-values are not pooled.

## Current level

High-value exploratory hypothesis. A reproducible cyst-wall/lesion-edge mask is not yet available, and it has not been shown that multicellular contribution predicts lesion location better than a single source.

## QC robustness

The fixed global sensitivity check contains 140 baseline candidates and 85 strict-QC candidates; LAM1/2/3/4 are 31/4/84/21 versus 28/1/52/4. This spatial source attribution has not yet been rerun separately on both candidate pools, so the card is not upgraded to a high-confidence mechanism.

## Next validation

1. improve source-state annotation;
2. report gene-specific and donor-specific protease contributions;
3. compare protease increase, antiprotease decrease and spatial separation independently;
4. test lesion-edge enrichment only after a reproducible lesion mask exists.
