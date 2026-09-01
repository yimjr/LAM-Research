# Stage 6 checkpoint

- Status: **GO**
- High-confidence candidate cells: 5378
- Independent patients: 12
- Reference LAM-only parameters: n_neighbors=30, resolution=0.4
- Reference LAM-only scVI latent clusters: 12
- Full-cohort clusters retained for audit only: 33
- Reference singleton clusters: 0
- Reference clusters with fewer than 5 cells: 0
- Reference median/largest cluster size: 440.5/865
- Reference clusters shared by at least 2 patients: 12
- Reference clusters shared by at least 2 patients with ≥5 cells per patient: 12
- Shared cluster across patients: True
- Dominant drivers: none detected
- scVI training mode: full

## Parameter grid
The nine configurations are reported as a stability analysis. Cluster count alone is not used to choose a configuration.

| grid_id | n_neighbors | resolution | clusters | singletons | <5 cells | largest | median | shared ≥2 patients | shared ≥2 patients and ≥5 cells/patient | patient ARI | dataset ARI | assay ARI |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| nn15_res0.2 | 15 | 0.2 | 7 | 0 | 0 | 1575 | 710 | 7 | 7 | 0.1409531021900585 | 0.13644983108405515 | 0.019676506322040134 |
| nn15_res0.4 | 15 | 0.4 | 11 | 0 | 0 | 1019 | 541 | 11 | 11 | 0.20326056979746596 | 0.1029630954631112 | 0.007797747823365938 |
| nn15_res0.6 | 15 | 0.6 | 19 | 0 | 0 | 614 | 214 | 18 | 18 | 0.17241171816958673 | 0.07833942774388557 | 0.012258519700917288 |
| nn30_res0.2 | 30 | 0.2 | 8 | 0 | 0 | 1522 | 657.5 | 8 | 8 | 0.1339826177967508 | 0.12170446121345098 | 0.011875622855931911 |
| nn30_res0.4 | 30 | 0.4 | 12 | 0 | 0 | 865 | 440.5 | 12 | 12 | 0.1525538138508663 | 0.08678601736572315 | 0.00909384843483104 |
| nn30_res0.6 | 30 | 0.6 | 14 | 0 | 0 | 870 | 342 | 14 | 14 | 0.15048303414648023 | 0.08727382437234399 | 0.007934603349720943 |
| nn50_res0.2 | 50 | 0.2 | 6 | 0 | 0 | 1589 | 796.5 | 6 | 6 | 0.1084171783349124 | 0.0936398536721832 | -0.010103072706577985 |
| nn50_res0.4 | 50 | 0.4 | 9 | 0 | 0 | 952 | 627 | 9 | 9 | 0.12326618114995518 | 0.09911177322479109 | 0.005527927092696521 |
| nn50_res0.6 | 50 | 0.6 | 13 | 0 | 0 | 867 | 403 | 13 | 12 | 0.1442937192368869 | 0.0862294407277758 | 0.009020065512002393 |

Pairwise partition ARI across grid: mean=0.7109048533752511, min=0.4205495031332101, max=0.942057569432587.

## Decision reasons
- high-confidence LAM latent structure is multi-cluster and cross-patient without a dominant measured driver

## Interpretation boundary
The full-cohort Leiden labels are retained for audit only. The Go/No-Go decision uses the reference LAM-only clustering of high-confidence cells; boundary and normal are auxiliary analyses.
Upstream candidate/state/program correspondence is post-hoc interpretation and is not a Go/No-Go criterion. Any unmatched scVI cluster is retained as a `novel_or_unexplained` candidate for later stages.
The parameter grid is intended to identify a stable parameter interval, not to optimize for a preferred number of clusters. Singleton/small-cluster counts and pairwise partition ARI should be reviewed before treating clusters as biological states.

Boundary auxiliary status: available; normal reference status: available.
