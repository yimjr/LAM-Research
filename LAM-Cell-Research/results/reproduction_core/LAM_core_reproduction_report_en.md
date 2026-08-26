# LAMCORE Core Lung Reproduction: Phase 1 Baseline

## Interpretation

This is an independent reimplementation guided by the authors' public R scripts, not a strict byte-for-byte Seurat reproduction. The public scripts mix Seurat 2 and Seurat 3, whereas the current executable path uses AnnData/Python; deviations are recorded in `method_deviation_table.csv`.

- The paper reported approximately 125 LAMCORE cells from LAM1, LAM3 and LAM4.
- Author-style marker/cluster reimplementation candidates: **140**.
- Candidate counts by donor: `{"LAM1": 31, "LAM2": 4, "LAM3": 84, "LAM4": 21}`.
- LAM2 candidates: **4**; this is an operational marker candidate count, not the paper's original LAMCORE label.

## Method boundary

Candidates were located using known features (PMEL, ACTA2, ESR1, FIGF/VEGFD, CTSK and MLANA) together with an author-style graph. The 777-gene signature was scored only after candidate selection as a consistency check, avoiding circular use of a signature derived from the original LAMCORE cells.

QC means downstream QC recoverable from processed matrices; FASTQ, Cell Ranger, initial barcode/cell calling and empty-droplet decisions were not reconstructed. Doublet scores and predictions were recorded, but Phase 1 did not remove cells on that basis.

## Next

The baseline is ready for targeted Phase 2 robustness checks and Phase 3 biological discovery in parallel. If the separate R/Seurat run narrows the implementation differences, the reproduction status will be updated.
