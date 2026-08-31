# Detailed methods appendix

## Input contract

The project inherited converted AnnData and upstream annotations from `LAM-Cell-Research` where available, with `data-temp` as the supplied fallback. GSE135851 inherited its upstream QC status; the three external AnnData objects received project-specific supplementary QC. The canonical candidate rule remained `pool_high_confidence`; broad-minus-high was boundary; unrestricted was audit-only.

## Matrix separation

NMF used raw counts followed by library-size normalization, log1p, State Modeling HVG selection and NMF. scVI used `layers["counts"]` as raw counts with `batch_key="dataset"` only. `assay` remained metadata. Later stages used the existing `X_scVI` and did not retrain.

## Statistical scope

Stage10 fitted each state independently against same-patient Rest_of_LAM using patient×state pseudobulk and `~ patient_id + group`. Stage8 LOO metrics were computed on retained cells only. Stage21/22 nulls were composition/cell-count matched exploratory controls, not replacements for independent biological replication.

## Stage24 boundary

Stage24 is a read-only synthesis stage. It performs deterministic reading, table joining and document audit. It does not run scanpy clustering, scVI training, DE, or new biological discovery.
