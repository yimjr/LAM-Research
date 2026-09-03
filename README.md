# LAM Research

A computational research project on **LAM (Lymphangioleiomyomatosis)** based primarily on publicly available biomedical datasets.

This project reuses and connects transcriptomic, single-cell, spatial omics, drug perturbation, and other public datasets to explore research questions related to **LAM pathogenesis, lung tissue destruction, therapeutic response, immune recognition, and cellular-state organization**.

The repository currently contains four major research directions. They are relatively independent, but can also provide complementary evidence for each other:

| Direction                                       | Core Question                                                                                                                       |
| ----------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| [LAM Cell Research](LAM-Cell-Research/)         | What cellular states exist in LAM, and how might LAM cells adapt to the lung and contribute to tissue destruction?                  |
| [LAM Drug Repurposing](LAM-Drug-Repurposing/)   | Can existing drugs reverse abnormal LAM/TSC2-loss states, including programs that remain after sirolimus treatment?                 |
| [LAM Immune Visibility](LAM-Immune-Visibility/) | Which LAM-associated features may be visible to the immune system, and how are antigen expression and antigen presentation related? |
| [LAM State Modeling](LAM-State-Modeling/)       | Can integrated single-cell data identify robust LAM-rich latent states, and do nearby states form reproducible biological trajectories or transitions? |

[中文](README_zh.md)

---

## Why This Project?

LAM is a rare disease in which **TSC1/TSC2 dysfunction and subsequent activation of the mTOR pathway** play a central role.

Sirolimus has substantially changed the treatment of LAM, but many biological questions remain open:

* LAM cells are heterogeneous and may exist in multiple cellular states;
* mTOR inhibition controls many abnormal programs, while some cellular programs may persist;
* the mechanisms that allow LAM cells to survive and cause destructive lesions in the lung remain incompletely understood;
* lung destruction may involve interactions among LAM cells, fibroblasts, immune cells, and protease systems;
* existing drugs may affect abnormal programs outside the canonical mTOR pathway;
* LAM cells express several unusual melanocytic and mesenchymal features, raising questions about their potential immune visibility;
* integrated single-cell state models may help distinguish reproducible LAM-rich structure from donor effects, ordinary lineage adjacency, and latent-space geometry.

At the same time, an increasing number of LAM single-cell, spatial, perturbation, and transcriptomic datasets have become publicly available.

The central idea of this project is therefore:

> **Reconnect public datasets originally generated for different research questions, and use them to identify new biological questions, candidate mechanisms, and experimentally testable hypotheses.**

The goal is to progressively reduce large public datasets into a smaller number of clear research directions that may be useful for future experimental studies.

---

# 1. LAM Cell Research

📁 [LAM-Cell-Research](LAM-Cell-Research/)

## Why Study LAM Cell States?

Single-cell studies have identified a characteristic population of LAM cells commonly referred to as **LAMCORE**, together with several cellular states.

However, this still leaves many questions:

* Are there additional LAMCORE states that have not been fully characterized?
* Do LAM cells from different patients share hidden gene-expression programs?
* Are LAM states best described as discrete subtypes, or as a continuous state space?
* How do LAM cells adapt to and survive in the lung?
* Do surrounding stromal and immune cells cooperate with LAM cells in tissue destruction?
* Which extracellular-matrix and protease-related programs remain active after mTOR inhibition?

This project therefore combines public single-cell data with spatial transcriptomics, perturbation datasets, and tissue-reference datasets.

---

## General Research Strategy

The overall workflow is:

```text
Public LAM single-cell datasets
        ↓
Identification of LAMCORE-like cells
        ↓
Reproduction of known LAMCORE states
        ↓
Analysis of within-patient and between-patient variation
        ↓
Discovery of candidate gene programs
        ↓
Interpretation using spatial, perturbation, and external datasets
```

In addition to individual genes, much of the analysis focuses on **gene programs**:

> Groups of genes that change together and may represent a shared cellular function or state.

This makes it possible to study broader cellular behavior rather than relying only on individual markers.

---

## Current Research Directions

### 1. Are LAMCORE States More Continuous Than Discrete?

The current implementation identifies approximately **140 LAMCORE-like candidate cells**, broadly comparable with the approximately 125 cells reported in the original study.

Under stricter quality-control criteria, approximately **85 cells** remain.

Further clustering analyses have not revealed a clearly separated new subtype that is consistently present across multiple donors.

For example, the HDBSCAN branch currently identifies:

* 2 local density clusters;
* 28 non-noise cells among 140 candidates;
* clusters mainly involving LAM3 and LAM4;
* no stable discrete subtype shared across all donors.

This has motivated an alternative possibility:

> **LAMCORE heterogeneity may be better described partly as continuous cellular variation rather than as a small number of sharply separated subtypes.**

The project therefore studies continuous states, local state structures, and donor-specific variation in addition to conventional clustering.

---

### 2. Searching for Under-Characterized LAM Gene Programs

The project uses pooled NMF, donor-wise NMF, meta-program matching, and related approaches to identify groups of genes that repeatedly vary together.

Initial candidate programs involve themes such as:

* ECM / extracellular matrix;
* protease activity;
* migration;
* proliferation;
* hormone-related signaling;
* LAM lineage;
* microenvironmental interactions.

So far, no single program has emerged as an obvious completely new LAMCORE program that is independently rediscovered across multiple patients.

This itself suggests an interesting possibility:

> Important LAM heterogeneity may arise from different combinations and intensities of known biological programs rather than from entirely separate cell types.

Several candidate programs are still being evaluated across additional LAM datasets.

---

### 3. Protease–Antiprotease Spatial Niche

One of the defining pathological features of LAM is progressive cystic destruction of lung tissue.

This raises an important question:

> **Does destructive protease activity come mainly from LAM cells themselves, or from a multicellular proteolytic niche involving LAMCORE cells, fibroblasts, immune cells, and other surrounding populations?**

The project currently analyzes spatial datasets generated using:

* Visium;
* Visium HD;
* Xenium.

Across these spatial technologies, LAMCORE-like spatial signals show directionally consistent associations with protease-related signals.

Current genes of interest include:

* `CTSK`
* `MMP` family members
* `ELANE`
* `PRTN3`
* `CTSS`

Antiprotease activity is analyzed in parallel to derive a broader concept of:

```text
protease activity
        -
antiprotease activity
        ↓
proteolytic balance
```

The long-term question is whether LAM lung destruction may involve a **multicellular local proteolytic environment** rather than a single destructive cell type.

---

### 4. Do LAM Cells Acquire a Lung-Adaptation Program?

LAM raises an unusual biological question:

> How do cells with a likely extrapulmonary origin successfully survive, expand, and form lesions in the lung?

To study this, the project compares:

```text
Normal uterus
LAM uterus
Normal lung
Pulmonary LAM
```

and evaluates contrasts conceptually similar to:

```text
Pulmonary LAM - Normal lung
        vs.
LAM uterus - Normal uterus
```

Current results suggest that:

* a `lung_adaptation` program is more consistent with an acquired pulmonary state;
* ECM-related programs are elevated in both LAM uterus and pulmonary LAM.

This motivates a model in which pulmonary LAM may combine:

```text
LAM transformation / lineage program
        +
lung-acquired adaptation program
        =
pulmonary LAM state
```

Understanding this distinction may help explain why LAM cells are particularly successful in the lung environment.

---

### 5. ECM / Protease Programs Retained After Rapamycin

Sirolimus and rapamycin strongly suppress mTOR-related growth signaling.

The project asks a further question:

> **When cellular growth is controlled, do matrix-remodeling and tissue-destruction programs disappear to the same extent?**

In TSC2-loss perturbation datasets, several ECM- and protease-related genes retain abnormal expression after rapamycin treatment.

Current examples include:

* `ELANE`, showing directionally consistent partial retention in both plastic and hydrogel conditions;
* `MMP2`, showing retention in the hydrogel condition.

This supports further study of the possibility that:

> **mTOR inhibition can strongly control growth while leaving parts of matrix-related pathology relatively active.**

This research line also connects directly with the drug-repurposing project.

---

# 2. LAM Drug Repurposing

📁 [LAM-Drug-Repurposing](LAM-Drug-Repurposing/)

## Why Study Drug Repurposing?

LAM already has an effective therapy in sirolimus.

However, sirolimus mainly targets the mTOR axis.

This raises two related questions:

> **Which parts of the abnormal TSC2-loss state remain after mTOR inhibition?**

and:

> **Can existing drugs shift those abnormal states back toward a WT-like state?**

Such compounds may serve as:

* mechanistic research tools;
* candidates for combination-treatment research;
* potential drug-repurposing directions.

---

## Evolution of the Research Strategy

The project originally started from:

```text
TSC2 loss
    ↓
rapamycin
    ↓
Which abnormalities remain?
```

Analysis later showed that residual programs differ substantially between plastic and hydrogel environments.

The strategy therefore expanded to:

```text
TSC2-loss transcriptional state
        ↓
Separate plastic / hydrogel disease signatures
        ↓
Comparison with large-scale LINCS/CMap perturbations
        ↓
Identification of drugs producing opposite expression patterns
        ↓
Candidate drug generation
        ↓
Target, mechanism, human-state, and perturbation validation
```

The analysis does not focus simply on restoring the expression level of the `TSC2` gene itself.

Instead, it uses:

> **The broader transcriptional state caused by TSC2 loss.**

A drug therefore does not need to directly target TSC2 to potentially reverse downstream consequences of TSC2 loss.

---

## Core Dataset: TSC2 × Rapamycin × Environment

GSE179044 provides a particularly useful factorial design:

```text
WT / TSC2-null
        ×
Vehicle / Rapamycin
        ×
Plastic / Hydrogel
```

The current analysis includes:

* 16 samples;
* 59,055 genes.

This makes it possible to separately ask:

1. What changes after TSC2 loss?
2. What does rapamycin restore?
3. What remains abnormal after rapamycin?
4. How does the extracellular environment modify these effects?

---

## Transcriptional States After Rapamycin

In the hydrogel condition, the current classification includes:

| Category                                 | Number of genes |
| ---------------------------------------- | --------------: |
| Near-complete rescue                     |           1,056 |
| Partial rescue with residual abnormality |           3,099 |
| Persistent residual                      |             668 |
| Worsened residual                        |             302 |
| Direction reversal                       |           1,385 |

One notable observation is that the canonical **mTORC1 program itself is relatively weak among the residual signals**.

The more persistent signals are enriched in areas such as:

* ECM / invasion;
* myogenic programs;
* metabolism;
* autophagy;
* stress-related programs.

This suggests that:

> **Persistent abnormalities after sirolimus may involve biological programs beyond incomplete suppression of mTORC1 itself.**

Hydrogel-related analyses have also highlighted genes such as:

* `NNMT`
* `COL8A1`
* `MIR210HG`
* `SLC40A1`
* `FBLN5`
* `DCN`
* `LUM`

as candidates for further mechanistic investigation.

---

## SRPK2 / ECM-Related Signals

Additional perturbation datasets provide a potentially interesting mechanistic connection.

In a related cellular model, rapamycin is associated with increased expression of genes including:

* `COL8A1`
* `NNMT`
* `DCN`
* `ACTA2`
* `MMP2`

while SRPK2 knockdown reduces several ECM- and metabolism-related genes including:

* `NNMT`
* `COL3A1`
* `LUM`
* `FBLN5`

This suggests a possible research axis:

```text
TSC2 loss
   ↓
mTOR inhibition
   ↓
Partial ECM / metabolic state persists
   ↓
SRPK2-sensitive program?
```

Genes such as `NNMT` and `COL8A1` are therefore being followed as potentially important mechanistic signals.

---

## Local LINCS / CMap Drug Screening

The project currently analyzes two LINCS Level 5 releases locally:

* GSE92742
* GSE70138

The basic idea is straightforward:

```text
LAM / TSC2-loss: A↑ B↑ C↓ D↓

Drug response:   A↓ B↓ C↑ D↑

               ↓

Potential reversal
```

Candidate generation does not require a drug to reverse both the plastic and hydrogel states simultaneously.

A compound can enter downstream analysis if it shows a credible reversal signal in either environment.

The current pipeline has produced:

* **258 candidate records**;
* corresponding to **66 unique candidate compounds** after removing dataset- and query-size duplicates.

Candidates are then evaluated according to:

* known drug targets;
* generic stress or cytotoxicity signatures;
* agreement with genetic perturbations;
* expression of relevant targets or programs in human LAM datasets;
* possible mechanistic complementarity with sirolimus.

The candidate set includes both PI3K/AKT/mTOR-related compounds and mechanisms outside the canonical mTOR pathway.

Accordingly, the project is increasingly interested in **non-mTOR mechanisms** rather than simply searching for additional mTOR inhibitors.

---

## Translation Programs: Another Layer of TSC2-Loss Biology

The project also analyzes GSE277844 to investigate whether TSC2 loss changes not only RNA abundance, but also:

> Which RNAs are preferentially translated into protein.

The current analysis identifies approximately:

* 89 translation-up genes;
* 107 translation-down genes;

for a total of **196 translation-abnormal genes**.

An interesting observation appears among genes overlapping the hydrogel residual program.

For 18 directly comparable hydrogel-residual genes:

* RMC-6272 moves 15/18 closer to WT;
* eFT-508 moves 15/18 closer to WT.

A shared set of **13 genes** is shifted toward WT by both perturbations:

```text
CDC42EP3
RND3
SERPINE2
GPR27
WWTR1
FBN2
REEP2
ZNF354C
PNMA2
CACFD1
NFATC4
SPIN4
GPC4
```

These genes collectively point toward themes including:

* ECM;
* cell adhesion;
* Rho GTPase signaling;
* TGF-β;
* migration.

Translation regulation may therefore represent another important layer of the abnormal TSC2-loss state.

---

# 3. LAM Immune Visibility

📁 [LAM-Immune-Visibility](LAM-Immune-Visibility/)

## Why Study Immune Visibility?

LAM cells express several distinctive lineage-associated genes, including:

* `PMEL`
* `MLANA`
* `MITF`
* `GPNMB`
* `TYRP1`
* `DCT`

Many of these genes are also relevant to melanocytic biology and tumor immunology.

This raises a natural question:

> **Do LAM cells express molecular features that could potentially be recognized by the immune system, while antigen processing, presentation, or the surrounding immune environment remains insufficient for effective recognition?**

The project therefore separates several biological layers:

```text
Antigen-related expression
       ↓
Antigen processing
       ↓
HLA / antigen presentation machinery
       ↓
Immune environment
       ↓
Potential immune visibility
```

---

## General Research Strategy

This project reuses public single-cell and spatial datasets already processed in the other research directions.

It examines several components separately.

### Antigen-Related Expression

Whether LAM cells consistently express candidate antigen-related genes.

### Antigen Presentation

Whether the same cells also express genes involved in:

* HLA;
* antigen processing;
* antigen presentation machinery.

### Immune Context

Whether candidate LAM states are spatially or transcriptionally associated with signals related to:

* T cells;
* NK cells;
* macrophages;
* immunosuppressive environments.

### Treatment Persistence

Whether selected LAM-associated features remain detectable after rapamycin perturbation.

---

## Initial Results

The current analysis has generated:

* **6 candidate antigen / lineage-marker rankings**;
* **735 patient-level module records**;
* **40 state-association records**;
* **20 immune-context associations**.

The current candidate-expression summary includes:

| Gene  | Proportion of patients with at least one detection |
| ----- | -------------------------------------------------: |
| MITF  |                                             100.0% |
| GPNMB |                                             100.0% |
| PMEL  |                                             100.0% |
| MLANA |                                              84.6% |
| TYRP1 |                                              69.2% |
| DCT   |                                              53.8% |

`MITF`, `GPNMB`, and `PMEL` are particularly consistently detected across the currently included patients.

This suggests that melanocytic and lineage-associated signals are recurring features of LAM rather than isolated observations.

A natural next set of questions is therefore:

```text
Are these proteins processed into antigenic peptides?
            ↓
Are those peptides presented by HLA?
            ↓
Can patient HLA genotypes present them?
            ↓
Can corresponding T-cell responses be detected?
```

Existing immunological knowledge around targets such as `PMEL/gp100` provides a useful reference point, while additional LAM-associated antigen candidates can be explored in parallel.

---

# 4. LAM State Modeling

📁 [LAM-State-Modeling](LAM-State-Modeling/)

## Why Model LAM States Separately?

The broader cell-research direction asks which gene programs and biological mechanisms vary across LAM cells. This project focuses more narrowly on a complementary question:

> **Can integrated single-cell data define reproducible LAM-rich latent states, and can the geometry around those states support stronger claims about continuous manifolds or biological transitions?**

The analysis inherits the processed AnnData objects, patient mappings, candidate pools, LAMCORE annotations, and upstream state/program information from `LAM-Cell-Research`, then builds a dedicated latent-state workflow without repeatedly rediscovering those upstream features.

The current implementation covers **Stage 1–24** and includes:

* harmonization and QC;
* PCA/NMF baselines;
* scVI latent-space modeling;
* consensus state construction;
* patient/dataset leave-one-out robustness;
* state hierarchy and biology annotation;
* candidate-identity audits;
* State15 anchor validation;
* global-manifold and local-branch testing;
* matched-null, patient-level, and LOPO robustness analyses.

---

## Current Findings

The current results support a deliberately conservative interpretation.

### State15 Is the Strongest LAM-Rich Consensus State

Within the current high-recall candidate pool, **State15 is the most LAM-rich frozen consensus state**.

Its profile is supported by formal LAMCORE signal, available author-label enrichment, patient-matched comparisons, and sensitivity analysis after removing the strongly represented LAM1163 donor.

However, State15 is not treated as a fully independent or evenly cross-patient reference anchor.

### A Single Global State15-Centered Manifold Is Not Robustly Supported

An initial pooled State15-centered gradient was visible, but stricter validation showed that it does not behave like one robust, unified LAM manifold across the data.

The project therefore does **not** interpret the latent space as a single established developmental or temporal trajectory.

### Local State16/12/20/7 Branches Do Not Support a LAM-to-Lineage Transition Claim

Stage22 identifies four main local neighboring branches around State15: **State16, State12, State20, and State7**.

After correcting branch eligibility, restricting real and null analyses to the same local 1–3-hop scope, matching local distance structure, using empirical tail probabilities, applying FDR correction, and adding patient-level/LOPO checks, none of these branches shows statistical evidence for a specific LAM-to-lineage transition beyond the matched local null.

State16 remains an interesting and reproducible local neighbor: patient-level slopes are directionally consistent, but the effect does not exceed distance-matched null expectations. It is therefore retained as an **ordinary lineage adjacency / mixed neighboring state**, not a transition candidate.

The current conclusion is therefore:

> **A reproducible LAM-rich core state can be identified, but the present data do not establish either a single global LAM trajectory or a specific local LAM-to-lineage transition.**

The project also does not claim temporal conversion, a diagnostic classifier, or that every candidate state represents true LAM identity.

---

# How the Four Directions Connect

The four projects can be studied independently.

At the same time, they examine the same disease from four complementary perspectives:

```text
                         LAM biology
                             │
       ┌─────────────────────┼─────────────────────┐
       │                     │                     │
       ▼                     ▼                     ▼
  Cell Programs         State Geometry        Drug Response
       │                     │                     │
What cells express      Which latent states   What can change
and how they interact   are reproducible      abnormal programs
       │                     │                     │
       └──────────────┬──────┴──────────────┬──────┘
                      │                     │
                      ▼                     ▼
             Immune Visibility      Testable Hypotheses
                      │
             What the immune
             system may see
```

`LAM-Cell-Research` and `LAM-State-Modeling` are especially complementary:

* **LAM Cell Research** emphasizes gene programs, spatial niches, lung adaptation, and biological mechanisms;
* **LAM State Modeling** tests whether integrated latent states and state-to-state geometry are reproducible enough to support stronger structural or transition hypotheses.

For example, a rapamycin-persistent ECM program discovered in **LAM Cell Research** can be further investigated by:

* locating the program within robust states from **LAM State Modeling**;
* searching for compounds that reverse it in **LAM Drug Repurposing**;
* locating it in spatial regions associated with lung destruction;
* examining whether the same cellular state has a distinct immune context in **LAM Immune Visibility**.

Conversely, a candidate drug mechanism can be mapped back into single-cell datasets to ask:

> Which reproducible LAM-rich state or program is most likely to respond to this perturbation?

The four directions therefore gradually converge on a broader framework:

> **Identify reproducible LAM states and programs, determine what can change them, and understand how they interact with lung tissue and the immune system.**

---

# Current Research Leads

Several questions currently appear particularly worth pursuing:

1. **LAMCORE heterogeneity may involve continuous and local variation, but the current state-modeling analysis does not support a single robust global LAM trajectory;**
2. **State15 is currently the strongest LAM-rich frozen consensus state, while its cross-patient balance and use as a formal reference anchor remain limited;**
3. **State16 is a reproducible local neighbor of State15, but current distance-matched null analyses do not support interpreting it as a LAM-to-lineage transition;**
4. **The lung may contain a multicellular protease–antiprotease spatial niche associated with LAM lesions;**
5. **LAM may combine lineage/transformation programs with lung-acquired adaptation programs;**
6. **ECM, protease, and metabolic programs may partly persist after rapamycin treatment;**
7. **ELANE and MMP2 may connect rapamycin persistence with mechanisms of lung tissue destruction;**
8. **NNMT, COL8A1, and an SRPK2-sensitive ECM/metabolic program deserve further investigation;**
9. **TSC2-loss transcriptional signatures have produced a substantial set of drug-repurposing candidates;**
10. **Translation regulation may represent an additional layer of the abnormal TSC2-loss state;**
11. **PMEL, MITF, GPNMB, and related lineage signals are consistently detected across multiple LAM patients;**
12. **The relationship among LAM antigen expression, antigen presentation, and immune context remains a promising research question.**

The common objective is to progressively reduce large public datasets into a smaller number of clearly defined questions that can be experimentally tested.

---

# Data Sources and Project Organization

The project primarily uses publicly available processed datasets, including:

* GEO transcriptomic datasets;
* single-cell and single-nucleus RNA-seq;
* spatial transcriptomics;
* Visium / Visium HD / Xenium;
* TSC2-loss perturbation datasets;
* rapamycin perturbation datasets;
* LINCS / CMap Level 5 perturbation data;
* selected translation / polysome datasets.

The typical analysis workflow is:

```text
Public dataset
    ↓
Standardized processing
    ↓
Gene / program / latent-state analysis
    ↓
Cross-dataset and robustness testing
    ↓
Candidate mechanism or constrained conclusion
    ↓
Testable hypothesis
```

Each subproject maintains its own:

* analysis scripts;
* data manifests;
* intermediate and final results;
* result tables;
* figures;
* research logs or reports.

---

# Repository Structure

```text
LAM-Research/
│
├── LAM-Cell-Research/
│   ├── single-cell / spatial analysis
│   ├── LAMCORE states
│   ├── gene program discovery
│   ├── protease spatial niche
│   ├── lung adaptation
│   └── rapamycin-persistent ECM/protease
│
├── LAM-Drug-Repurposing/
│   ├── TSC2-loss / rapamycin factorial analysis
│   ├── residual programs
│   ├── LINCS / CMap
│   ├── drug candidate analysis
│   ├── mechanism integration
│   └── translation analysis
│
├── LAM-Immune-Visibility/
│   ├── antigen-related expression
│   ├── antigen presentation
│   ├── immune context
│   ├── candidate antigen ranking
│   └── hypothesis cards
│
└── LAM-State-Modeling/
    ├── scVI latent-state modeling
    ├── consensus and robustness analysis
    ├── candidate-identity audits
    ├── State15 anchor validation
    ├── manifold / local-branch testing
    └── final state-modeling reports
```

For detailed methods, execution instructions, and result files, see the README and research documents within each subdirectory.

---

# Project Status

This repository is under active development.

The four current directions are approximately at the following stages:

| Project               | Current Stage                                                                                                          |
| --------------------- | ---------------------------------------------------------------------------------------------------------------------- |
| LAM Cell Research     | Reproduction baseline established; multiple new biological questions under active exploration                          |
| LAM Drug Repurposing  | Major TSC2-loss / LINCS candidate-generation stage completed; mechanism and cross-dataset validation ongoing           |
| LAM Immune Visibility | Initial computational analysis and candidate-antigen ranking completed; follow-up validation questions being developed |
| LAM State Modeling    | Stage 1–24 analysis completed; State15 retained as the strongest LAM-rich consensus state, while global-manifold and local-transition claims are not supported by the final robustness analyses |

The project will continue to evolve as new public LAM datasets become available and existing datasets are reanalyzed from new perspectives.

---

## About the Results

The results in this repository are primarily derived from computational analysis of public datasets and are intended to generate **research hypotheses**, **candidate mechanisms**, and well-defined negative or constraining results.

They are best interpreted as starting points for further experimental and mechanistic research and are not intended as clinical treatment recommendations.

---

## License

The source code in this repository is licensed under the [Apache License 2.0](LICENSE).

Copyright and attribution notices are provided in [NOTICE](NOTICE).

Third-party datasets, software, publications, and other materials included or referenced by this project remain subject to their respective original licenses and terms.
