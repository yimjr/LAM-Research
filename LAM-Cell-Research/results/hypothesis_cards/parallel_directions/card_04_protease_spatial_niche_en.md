# Hypothesis Card 04: LAMCORE–protease spatial niche

## Observation

Across Visium, Visium HD and Xenium in GSE302356, spatial units with higher LAMCORE-like module scores showed positive neighborhood enrichment for the protease module. LAMCORE–protease unit-level Spearman correlations were approximately 0.818, 0.402 and 0.895, with top-decile neighborhood enrichment ratios of 6.66, 3.96 and 7.44 for Visium, Visium HD and Xenium, respectively.

## Interpretation boundary

This is a directionally concordant spatial candidate relationship across three technologies. It is not proof of cellular communication and does not establish that LAMCORE cells themselves produce the proteases. Visium is spot-level, Visium HD is bin/segment-level, and Xenium is targeted cell-level; raw units, scores and p-values were not pooled. A missing gene in Xenium is panel-unobserved, not a biological negative.

## Candidate mechanism

LAMCORE, immune cells, fibroblast/LAF-like cells and ECM may form a local proteolytic/ECM-remodeling niche that contributes to cystic tissue destruction. The next step is to resolve source contributions from LAMCORE-associated CTSK/MMP genes, immune proteases, fibroblast ECM programs and antiprotease localization.

## Robustness and evidence level

- PatientIDs: LAM4 (Visium LAM20 plus Xenium LAM19) and LAM3 (Visium HD LAM18); these are not two independent patients.
- Orthogonal support: directionally concordant across three spatial technologies.
- QC gate: spatial units are not identical to the single-cell 140→85 QC gate; any cell-level upgrade still requires baseline/strict-QC sensitivity.
- Current level: exploratory hypothesis, not a high-confidence novel program.

## Testable predictions

1. LAMCORE-rich regions should also show MMP/CTSK activity and ECM breakdown/remodeling markers;
2. protease sources should partition among LAMCORE, immune and/or LAF-like cells rather than a single source;
3. spatial separation or co-localization of antiproteases and proteases may distinguish lesion regions better than global expression means.

## Files

See the three `*_unit_scores.csv`, three `*_co_localization.csv` and modality manifests under `results/spatial/GSE302356/`.
