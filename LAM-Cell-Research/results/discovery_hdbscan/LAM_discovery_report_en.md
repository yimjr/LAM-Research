# New Biological Exploration of LAMCORE (Phase 3)

## Interpretation

This phase starts from the author-style marker candidates and explores continuous expression programs, the LAM2 difference and microenvironment associations. Results are research leads; algorithmic structures are not called new subtypes or physical cell communication.

Candidate cells: 140; candidate donors: LAM1, LAM2, LAM3, LAM4.

## HDBSCAN result

Under the primary setting, HDBSCAN identified 2 density clusters among 140 candidates, assigned 28 cells, and marked 112 as noise across 2 donors.
The primary non-noise cells came only from LAM3 and LAM4; LAM1 and LAM2 did not enter a density cluster. This does not support a stable discrete state spanning all four donors among the 140 candidates; the state should remain framed as continuous or locally structured until donor-wise validation is completed.

## Main outputs

- `lamcore_state_programs_by_donor.csv`: donor-level state differences;
- `lamcore_state_heterogeneity_hdbscan.csv`: HDBSCAN sensitivity in continuous state space;
- `lamcore_state_hdbscan_cluster_summary.csv`: HDBSCAN cluster summaries;
- `candidate_microenvironment_associations.csv`: candidate expression associations, not direct communication;
- `external_validation_status.json`: processing status for normal-lung GSE122960 and mouse-uterus GSE118180;
- `results/hypothesis_cards_hdbscan/`: three bilingual research-lead cards.

The highest-priority lead is continuous LAMCORE state variation, followed by the weaker/heterogeneous LAM2 signal and candidate lymphatic/ECM associations. Independent LAM-donor validation remains outstanding.
