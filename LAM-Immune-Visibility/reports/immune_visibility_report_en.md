# First-pass LAMCORE Immune Visibility Report

## Positioning

This is a computational analysis of existing public processed matrices. It does not establish clinical efficacy, peptide presentation, or direct cellular communication.

## Expression-state rule

Any nonzero raw count is detected. `detected_low` means detected and located in a predefined low-expression interval; low/high describe the expression interval after detection. Genes absent from a matrix or panel are `not_assayed`, not negative.

## Current outputs

The run produced 6 antigen-associated candidates, 735 patient-level module rows, 40 state-association rows, and 20 immune-context association rows.

## Interpretation

Antigen modules represent antigen-associated expression, while presentation modules represent machinery expression. Without immunopeptidomics, no specific HLA peptide is confirmed. A single-cell zero is reported as not detected in the current assay, not as true biological absence.
Spatial and immune associations are candidate associations only; same-patient modalities do not increase independent donor counts. Rapamycin retention is perturbation-supported persistence evidence, not patient-level treatment causality.