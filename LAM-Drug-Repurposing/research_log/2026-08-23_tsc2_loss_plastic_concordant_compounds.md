# `tsc2_loss_plastic` replicated-concordant compounds — 2026-08-23

## Filter and counting

Exact filter:

```text
contrast = tsc2_loss_plastic
perturbation_class = compound
cross_phase_status = replicated_concordant
```

The filter returns 92 dataset/query-size/perturbation summary rows, representing 29 unique `pert_iname` values. The 92 rows must not be described as 92 distinct drugs.

Among the 92 rows, 60 are positive `reversal_direction` and 32 are negative `mimic_direction`. Among the 29 unique drugs, 20 are reversal-direction only and 9 are mimic-direction only. `replicated_concordant` means the two LINCS releases agree in direction category; it does not mean that a compound reverses the disease signature.

## Target annotation scope

The primary target table uses ChEMBL exact InChIKey/name matching and ChEMBL curated mechanism records. Twelve compounds had no ChEMBL mechanism record, so a small set of well-established mechanisms was added from explicitly labeled primary-literature fallbacks. The final table covers all 29 drugs with 188 target records:

- 161 ChEMBL curated mechanism records;
- 27 literature-fallback records;
- no drug is marked as “no target”; missing ChEMBL mechanism is not interpreted as biological absence of targets.

The target table is a pharmacological target table, not an exhaustive list of every weak assay hit. PubChem BioAssay and BindingDB are retained separately as evidence layers.

## Drug-level target map

The following is the compact drug-level map. The machine-readable table contains target IDs, UniProt accessions where available, mechanism, action type, evidence source and source URL.

| Drug | Main target interpretation |
|---|---|
| AZD-8055 | MTOR |
| CGS-20625 | GABA-A receptor complex: GABRA1/GABRB2/GABRG1 |
| GDC-0941 | Class I PI3K: PIK3CA/B/C/D and regulatory PIK3R1/2/3/5 |
| GSK-1059615 | Class I PI3K: PIK3CA/B/C/D and PIK3R1/2/3/5 |
| GSK-2126458 | Class I PI3K plus MTOR |
| LY-294002 | Class I PI3K plus reported mTOR activity; tool compound with off-target risk |
| MG-132 | Proteasome, with PSMB5 as the compact primary gene-level representation; calpain/cysteine-protease off-targets are possible |
| NVP-BEZ235 | Class I PI3K plus MTOR |
| OSI-027 | MTOR |
| PI-103 | PIK3CA, MTOR and PRKDC/DNA-PK evidence |
| QL-X-138 | BTK plus MKNK1/MKNK2 |
| ZSTK-474 | Class I PI3K |
| bortezomib | 26S proteasome complex |
| buparlisib | Class I PI3K |
| clofarabine | RRM1, POLA1; DCK is an activation enzyme |
| colchicine | Tubulin/microtubule complex |
| dinoprost | PTGFR |
| gatifloxacin | Bacterial DNA gyrase/topoisomerase IV; not a human LAM target |
| ixazomib | 26S proteasome complex |
| lacidipine | L-type calcium channels: CACNA1C/D/F/S |
| lestaurtinib | FLT3, JAK2, RET and NTRK1/2/3 |
| letrozole | CYP19A1/aromatase |
| milnacipran | SLC6A2 and SLC6A4 |
| niacin | HCAR2 primary; HCAR3 lower-affinity secondary receptor |
| pevonedistat | NAE1/UBA3 NEDD8-activating enzyme |
| phorbol-myristate-acetate | PKC family, represented by PRKCA/PRKCD |
| podophyllotoxin | Tubulin, represented by TUBB |
| roxithromycin | Bacterial 50S ribosome/23S rRNA; not a human gene target |
| sirolimus | Direct FKBP1A mechanism record; pharmacology is the FKBP12–mTORC1 complex |

## Main interpretation

1. The strongest positive reversal cluster is dominated by PI3K/mTOR compounds: ZSTK-474, GSK-2126458, QL-X-138, PI-103, OSI-027, buparlisib, GSK-1059615, GDC-0941, LY-294002, AZD-8055, NVP-BEZ235 and sirolimus are all positive-direction examples. This is compatible with a residual program containing a TSC2–PI3K–mTOR-related component, but does not establish a drug-specific LAM vulnerability.
2. Several positive-direction compounds are pharmacologically non-LAM-specific or non-human-target controls: gatifloxacin, roxithromycin, milnacipran, niacin, letrozole and dinoprost. Their connectivity should be treated as hypothesis-generating until human LAMCORE/LAF expression, target perturbation and generic-cytotoxicity filters are applied.
3. The mimic-direction group is enriched for proteasome, microtubule and stress-inducing compounds: bortezomib, ixazomib, MG-132, colchicine, podophyllotoxin, pevonedistat, PMA, clofarabine and lacidipine. Their cross-release agreement is real at the LINCS direction-category level, but it is not evidence that they are useful combination partners with sirolimus.
4. BindingDB is most useful here for showing that broad kinase inhibitors can have large panels of weak or context-dependent associations. GDC-0941, PI-103 and lestaurtinib return especially broad target panels; these should be filtered by affinity strength, human LAM-state expression and target perturbation concordance rather than counted as 300 independent mechanisms.
5. DrugBank, TTD and STITCH are not included in the primary pipeline: DrugBank access/licensing is restrictive, TTD is better suited to disease-target annotation than compound-level binding completeness, and STITCH associations mix evidence and inference. They can be added later as exploratory cross-references, but would not strengthen the current primary target calls enough to justify the extra complexity.

## Reproducibility outputs

- `results/candidates/tsc2_loss_plastic_replicated_concordant_compounds_92_rows.csv`
- `results/candidates/tsc2_loss_plastic_replicated_concordant_compounds_29_unique.csv`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_compound_targets.csv`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_target_summary.csv`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_target_family_summary.csv`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_drug_target_analysis.csv`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_pubchem_identity.csv`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_pubchem_assay_target_evidence.csv.gz`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_pubchem_assay_target_summary.csv`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_bindingdb_affinity_evidence.csv.gz`
- `data/processed/candidate_analysis/drug_targets/tsc2_loss_plastic_replicated_concordant_bindingdb_target_summary.csv`
- `scripts/analyze_tsc2_loss_plastic_concordant_compounds.py`

## Source notes

- ChEMBL REST API: `https://www.ebi.ac.uk/chembl/api/data`
- PubChem PUG REST: `https://pubchem.ncbi.nlm.nih.gov/rest/pug`
- BindingDB REST `getTargetByCompound`: `https://bindingdb.org/rest/getTargetByCompound`
- Fallback literature URLs are stored row-wise in `target_source_url` and `mechanism_reference`.
