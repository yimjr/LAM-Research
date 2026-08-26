# Hypothesis Card 06: ECM/protease programs retained after mTOR inhibition

## Observation

In the TSC2-knockout model of GSE179044, rapamycin weakly suppressed the ECM-remodeling and protease/ECM-niche scores, with similar directions in hydrogel and plastic conditions. In the TSC2-null LAM-derived cells of GSE84476, ECM-remodeling and protease/ECM-niche scores were higher under rapamycin than siCtrl.

## Interpretation boundary

This is a cross-perturbation candidate mechanism, not patient-level sirolimus-resistance evidence. GSE179044 has two biological replicates and supports replicate-aware gene-level follow-up; GSE84476 has limited samples per condition and is descriptive here. Absolute scores cannot be compared directly across matrices and scales.

## Candidate mechanism

mTOR inhibition may reduce parts of growth/inflammatory activity without eliminating ECM adaptation or protease-associated programs, which may even become relatively stronger in some TSC2-null models. This is a candidate “growth suppression with matrix adaptation retained” mechanism.

## Evidence level

High-value exploratory hypothesis. It is not yet patient-level treatment persistence or a novel mechanism. It requires gene-level replicate consistency, mapping to human LAM states and ECM-conditioned validation.

## Testable predictions

1. retained ECM/protease genes should show concordant directions across both GSE179044 replicates;
2. retention should be stronger under ECM hydrogel than plastic, or show a clear genotype-by-environment interaction;
3. the program should map to human LAMCORE/spatial protease niches rather than being restricted to the in-vitro models.

## Files

See `GSE179044_program_contrasts.csv`, `GSE84476_program_contrasts.csv` and the corresponding analysis manifests under `results/perturbation/`.
