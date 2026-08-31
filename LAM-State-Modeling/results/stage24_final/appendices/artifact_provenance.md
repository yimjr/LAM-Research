# Artifact provenance appendix

The authoritative machine-readable provenance is `artifact_index.csv`. It lists all project scripts, reports, result tables, AnnData/model files and declared upstream paths, with existence, file size, CSV row/column counts where applicable, and SHA-256 for small files. Large binary artifacts are intentionally not hashed during Stage24 to avoid unnecessary I/O and memory pressure.

Core frozen artifacts include:

- `data/processed/state_model_scvi.h5ad` and `data/processed/scvi_model/`
- `results/stage7/state_consensus_assignments.csv`
- `results/stage7/state_consensus_state_summary.csv`
- `results/stage13/state_atlas.csv`
- `results/stage18/state15_anchor_summary.json`
- `results/stage20/state15_centered_manifold.csv`
- `results/stage21/gradient_models.csv`
- `results/stage22/branch_evidence_summary.csv`
- `results/stage23_visualization/visualization_manifest.json`

Stage24-generated files are confined to `results/stage24_final/`.
