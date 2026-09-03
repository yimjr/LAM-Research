# LAM Research

A computational research project on **LAM (Lymphangioleiomyomatosis)** based primarily on publicly available biomedical datasets.

This project reuses and connects transcriptomic, single-cell, spatial omics, drug perturbation, and other public datasets to explore research questions related to **LAM pathogenesis, lung tissue destruction, therapeutic response, immune recognition, and cellular heterogeneity**.

The repository currently contains four major research directions. They are relatively independent, but can also provide complementary evidence for each other:

| Direction | Core Question |
| --- | --- |
| [LAM Cell Research](LAM-Cell-Research/) | What cellular programs characterize LAM cells, and how might they adapt to the lung and contribute to tissue destruction? |
| [LAM Drug Repurposing](LAM-Drug-Repurposing/) | Can existing drugs reverse abnormal LAM/TSC2-loss states, including programs that remain after sirolimus treatment? |
| [LAM Immune Visibility](LAM-Immune-Visibility/) | Which LAM-associated features may be visible to the immune system, and how are antigen expression and antigen presentation related? |
| [LAM State Modeling](LAM-State-Modeling/) | Can LAM cells be identified more reproducibly across multiple single-cell datasets, and what cellular states and biological features characterize them? |

LAM State Modeling introduces a methodological component that differs from the other directions: it uses **scVI, a neural-network-based single-cell model**, to integrate high-dimensional gene-expression information across datasets and patients.

Its main research objective is to improve the identification of LAM-rich cell populations and then study their biological characteristics and internal heterogeneity.

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
* LAM cells express several unusual melanocytic and mesenchymal features, raising questions about their potential immune visibility.

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

## Why Build a LAM State Model?

A fundamental challenge in LAM single-cell analysis is the identification of the disease cells themselves.

LAM cells are rare, their expression varies across patients, and individual markers can be affected by sequencing depth, dropout, biological state, and overlap with other lung cell populations.

A broad candidate-selection strategy improves sensitivity, but it can also retain cells with substantial features of normal lung lineages.

This project therefore asks:

> **Can multiple LAM single-cell datasets be integrated to identify reproducible LAM-rich cell populations and then characterize the biological variation within and around those populations?**

The analysis combines information from thousands of genes rather than relying on a small marker panel.

---

## Neural-Network Integration with scVI

The project uses **scVI**, a neural-network-based probabilistic model developed for single-cell RNA-seq analysis.

Each cell contains expression measurements for thousands of genes. scVI learns a lower-dimensional representation that summarizes this high-dimensional information while accounting for major differences among datasets.

Conceptually:

```text
Expression of thousands of genes
        ↓
scVI neural-network integration
        ↓
Low-dimensional representation of each cell
        ↓
Repeated neighborhood and clustering analysis
        ↓
Reproducible consensus states
        ↓
LAM-identity and biological validation
```

This representation provides a common space in which cells from different patients and datasets can be compared using their broader transcriptional profiles.

---

## General Research Strategy

The current workflow begins with a high-confidence candidate pool inherited from the upstream LAM single-cell analyses.

The overall process is:

```text
Candidate LAM cells from multiple datasets
        ↓
Quality control and harmonization
        ↓
scVI integration
        ↓
Repeated clustering across parameter settings
        ↓
Consensus state construction
        ↓
Patient- and dataset-level robustness analysis
        ↓
LAM-identity audit
        ↓
Identification of LAM-rich states
        ↓
Biological characterization
        ↓
Analysis of local and global state relationships
```

Several forms of evidence are considered together:

* stability across clustering parameters;
* reproducibility across patients and datasets;
* formal LAMCORE expression;
* known LAM-associated markers and gene programs;
* original study annotations where available;
* competing normal-lung lineage signals;
* patient-level sensitivity analyses.

This allows the project to distinguish reproducible state structure from patterns driven mainly by a single dataset, patient, or ordinary lung lineage.

---

## Building Reproducible Consensus States

The integrated analysis contains more than five thousand high-confidence candidate cells from multiple patients and datasets.

The scVI representation is analyzed under multiple combinations of neighborhood size and clustering resolution.

Instead of relying on one clustering configuration, the project measures how often pairs of cells remain grouped together across different parameter settings and seeds.

These repeated relationships are summarized into **consensus states**.

Patient- and dataset-level leave-one-out analyses are then used to evaluate whether the same state structure remains detectable when individual sources of data are removed.

This creates a more stable basis for downstream biological interpretation.

---

## Re-Evaluating LAM Identity

A major step in the project was the recognition that the inherited candidate pool contains substantial biological heterogeneity.

Some consensus states show strong LAM-associated expression, while others resemble endothelial, fibroblast, immune, epithelial, or other lung lineages.

The project therefore added a dedicated identity-analysis branch.

Each state was evaluated using:

* the formal LAMCORE program;
* melanocytic and LAM-associated markers;
* supportive LAM-related programs;
* competing lineage programs;
* original annotations from the source datasets;
* patient-level reproducibility.

This analysis shifted the emphasis from simply cataloguing consensus states to identifying which states carry the strongest evidence of LAM-cell identity.

---

## Identification of a LAM-Rich Consensus State

Within **LAM State Modeling**, the consensus state with the strongest combined LAM evidence is referred to as **State15**.

State15 shows:

* strong formal LAMCORE expression;
* enrichment for available source-study LAM annotations;
* elevated LAM-associated marker and program signals;
* consistent patient-matched differences from comparison states.

One patient contributes a relatively large proportion of State15 cells, so the analysis was repeated after excluding that patient.

The LAM-rich profile remained detectable.

These analyses support State15 as the strongest current **LAM-rich consensus population** in this modeling framework.

Its patient distribution remains uneven, so the state is used as a high-confidence LAM-rich reference population for downstream analysis rather than as a complete definition of all LAM cells.

---

## Exploring the Structure Around LAM-Rich Cells

After identifying State15, the project examined whether nearby cells form systematic biological structure.

The first analysis considered the broader geometry around State15 and tested whether LAM-related expression changed continuously with increasing distance from the LAM-rich state.

A pooled gradient was detectable.

Further analyses examined:

* independent LAMCORE measurements;
* patient-specific gradients;
* dataset consistency;
* matched null models;
* local graph connectivity.

These results suggested that part of the LAM-related signal varies continuously in the surrounding latent space, while the complete dataset does not support a single uniform global trajectory.

This motivated a more local analysis of the immediate neighborhood around State15.

---

## Local Neighboring States

Within **LAM State Modeling**, four consensus states form the major local branches around State15:

* State16;
* State12;
* State20;
* State7.

State16 showed the clearest patient-level gradient.

Among patients with sufficient cells, LAMCORE decreased consistently with increasing distance from State15, and the same direction remained stable in leave-one-patient-out analyses.

The project then introduced stricter matched-null comparisons that controlled for:

* patient composition;
* dataset composition;
* local graph distance;
* continuous distance from State15.

After this matching, the State16 gradient fell within the range expected from comparable local structures.

The same framework was applied to the other neighboring states.

The current interpretation is therefore that these states describe reproducible local organization around the LAM-rich population, while their biological meaning remains open for further investigation.

---

## Current Interpretation

The modeling analysis currently supports several main conclusions:

* integrated multi-dataset analysis can identify a reproducible LAM-rich consensus population;
* LAM identity is unevenly distributed across the original broad candidate pool;
* consensus-state stability and biological identity provide complementary information;
* LAM-related expression varies across the local neighborhood of the LAM-rich state;
* local neighboring states show reproducible structure across patients;
* the biological meaning of these neighboring structures requires additional orthogonal evidence.

These results provide a framework for linking cell identification, state heterogeneity, spatial context, and downstream functional analyses.

---

# How the Four Directions Connect

The four projects examine LAM at different biological levels:

```text
                         LAM biology
                             │
        ┌────────────────────┼────────────────────┐
        │                    │                    │
        ▼                    ▼                    ▼
 Cell Programs        Cell Identification      Drug Response
 and Mechanisms       and State Modeling
        │                    │                    │
What LAM cells        Which cells form         What can alter
express and do        reproducible LAM-rich    abnormal programs
                      populations
        │                    │                    │
        └─────────────┬──────┴─────────────┬─────┘
                      │                    │
                      ▼                    ▼
              Immune Visibility    Testable Mechanisms
```

**LAM Cell Research** focuses on the biological programs expressed by LAM cells, including lung adaptation, extracellular-matrix remodeling, protease activity, and microenvironmental interactions.

**LAM State Modeling** focuses on identifying reproducible LAM-rich populations across multiple single-cell datasets and defining the structure of variation within and around those populations.

The resulting state framework can support the other directions:

* gene programs from **LAM Cell Research** can be mapped onto reproducible LAM-rich populations;
* drug-response signatures from **LAM Drug Repurposing** can be evaluated in specific LAM-rich states;
* antigen and presentation programs from **LAM Immune Visibility** can be examined in the same cellular populations.

The four directions therefore converge on a broader objective:

> **Identify LAM cells more reliably, understand their biological programs and heterogeneity, determine which perturbations may alter those programs, and study how these cells interact with lung tissue and the immune system.**

---

# Current Research Leads

Several hypotheses currently appear particularly worth pursuing:

1. **LAM-cell heterogeneity may arise from continuous variation and different combinations of shared biological programs rather than from a small number of sharply separated subtypes;**
2. **A more selective cross-patient LAM-cell identity model may reveal additional disease-specific states that are obscured within broad high-recall candidate pools;**
3. **The reproducible local states surrounding the LAM-rich population may reflect interactions between LAM cells and specific stromal, endothelial, epithelial, or immune environments;**
4. **Spatial transcriptomic data may help determine whether these neighboring states occupy reproducible anatomical positions around LAM lesions;**
5. **The lung may contain a multicellular protease–antiprotease spatial niche associated with LAM lesions;**
6. **LAM may combine lineage/transformation programs with lung-acquired adaptation programs;**
7. **ECM, protease, and metabolic programs may partly persist after rapamycin treatment;**
8. **ELANE and MMP2 may connect rapamycin-persistent programs with mechanisms of lung tissue destruction;**
9. **NNMT, COL8A1, and an SRPK2-sensitive ECM/metabolic program deserve further investigation;**
10. **TSC2-loss transcriptional signatures may identify drug-repurposing candidates with mechanisms complementary to mTOR inhibition;**
11. **Translation regulation may represent an additional layer of the abnormal TSC2-loss state;**
12. **LAM-associated lineage antigens and antigen-presentation programs may define subsets of LAM cells with different immune visibility.**

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
Gene / program-level analysis
    ↓
Cross-dataset comparison
    ↓
Candidate mechanism
    ↓
Hypothesis Card
```

Each subproject maintains its own:

* analysis scripts;
* data manifests;
* intermediate results;
* result tables;
* figures;
* research logs;
* hypothesis cards.

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
    ├── neural-network single-cell integration
    ├── consensus state modeling
    ├── cross-patient robustness analysis
    ├── LAM identity analysis
    └── state characterization
```

For detailed methods, execution instructions, and result files, see the README and research documents within each subdirectory.

---

# Project Status

This repository is under active development.

The four current directions are approximately at the following stages:

| Project | Current Stage |
| --- | --- |
| LAM Cell Research | Reproduction baseline established; multiple biological questions under active exploration |
| LAM Drug Repurposing | Major candidate-generation analysis completed; mechanism and cross-dataset validation ongoing |
| LAM Immune Visibility | Initial computational analysis completed; follow-up validation questions being developed |
| LAM State Modeling | Main modeling workflow completed; biological interpretation and cross-direction integration ongoing |

The project will continue to evolve as new public LAM datasets become available and existing datasets are reanalyzed from new perspectives.

---

## About the Results

The results in this repository are primarily derived from computational analysis of public datasets and are intended to generate **research hypotheses** and **candidate mechanisms**.

They are best interpreted as starting points for further experimental and mechanistic research and are not intended as clinical treatment recommendations.

---

## License

The source code in this repository is licensed under the [Apache License 2.0](LICENSE).

Copyright and attribution notices are provided in [NOTICE](NOTICE).

Third-party datasets, software, publications, and other materials included or referenced by this project remain subject to their respective original licenses and terms.
