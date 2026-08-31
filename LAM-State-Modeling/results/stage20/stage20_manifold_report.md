# Stage 20：State 15-centered transcriptional manifold

本阶段固定使用现有 State 15 的 200 个细胞和 `X_scVI`，分析对象为 5,378 个 high-confidence candidates 加 16,883 个 inherited boundary cells。没有重训 scVI、重新 Leiden/consensus clustering 或修改 candidate gate。

## Frozen input

- State 15 cells: 200
- State 15 ID SHA-256: `bb060490484a7e99ff50675727b6111dd387db51c9ea6eaa1d962626ffeaa5ef`
- Latent artifact: `/mnt/e/lam-research/LAM-State-Modeling/data/processed/state_model_scvi.h5ad` (`X_scVI`)
- Main geometry cells: 22261 (5378 candidates + 16883 boundary)

## Geometry-first distance axis

距离、State15 邻居比例和 30-NN 图均只由 `X_scVI` 计算；LAMCORE、program 和 lineage score 在几何结果生成后才作为描述性映射。

- State 15-nearest distance median: 0.0000
- State 15 neighbor fraction median: 1.0000

## Distance-bin composition

详见 `state15_distance_bins.csv`；该表只报告当前 State 1–20 标签和 boundary 的几何组成，不重新定义 state。

## Identity and lineage gradients

- LAMCORE median, nearest bin → farthest bin: 0.2152 → 0.1556.
- Patient-level distance~LAMCORE Spearman rho < 0 in 8/12 patients with estimable values.
- Dataset-level distance~LAMCORE Spearman rho < 0 in 4/4 datasets with estimable values.

## State 16 audit

- State 16 cells audited: 396.
- Raw-count LAM-marker/immune-marker coexpressing State 16 cells: 80.
- LAM identity threshold: 1.1012; immune threshold: 0.3553.
详见 `state16_lam_immune_coexpression.csv` 和 `state16_doublet_audit.csv`；阈值只用于诊断分类，不改变 State 16 标签。

## Boundary and normal scope

Boundary 仅投射到 State 15-centered geometry，normal 仅作为远端对照，不参与 State 15 邻域图或 state 数量。

       cohort  n_cells  n_datasets  n_patients  nearest_state15_distance_median  nearest_state15_distance_q25  nearest_state15_distance_q75  nearest_state15_distance_mean  state15_centroid_distance_median
normal_remote    32977           1           6                         4.821811                      4.484746                      5.198003                       4.869757                          5.940022

## Stage 20 checkpoint

- Current geometry checkpoint: `supports_lam_centered_transcriptional_manifold`.
- Interpretation: LAM identity signals decline along the State 15 distance axis with majority patient- and dataset-level negative distance correlations.

## Outputs

- boundary_state15_projection.csv
- dataset_gradient_consistency.csv
- normal_remote_summary.csv
- patient_gradient_consistency.csv
- stage20_manifest.json
- stage20_manifold_report.md
- state15_cell_distances.csv
- state15_centered_manifold.csv
- state15_distance_bins.csv
- state15_identity_gradient.csv
- state15_lineage_gradient.csv
- state16_cell_audit.csv
- state16_doublet_audit.csv
- state16_lam_immune_coexpression.csv
- state16_lam_immune_per_cell.csv
