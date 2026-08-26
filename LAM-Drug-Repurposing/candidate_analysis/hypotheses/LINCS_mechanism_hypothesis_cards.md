# LINCS stable-module mechanism cards

These are working hypotheses, not Tier 1 treatment recommendations. The
modules were defined by cross-release gene clustering and Jaccard matching;
the two LINCS releases are cross-phase/cross-release recurrence, not
independent biological replication.

## H1 — mimic_only__stable_1

Working hypothesis: this is a mixed proteostasis/UPR, oxidative-stress and
inflammatory-metabolic program. It is partly generic perturbation, but the
module is not completely removed by the operational local stress templates.

Evidence: cross-release Jaccard 0.614 and mean post-adjustment direction
retention about 0.95 in GSE92742 and 0.98 in GSE70138. Leading genes include DDIT3, SQSTM1,
LAMP3, TRIB3, MT3 and VEGFA. The enrichment output points toward stress,
protein catabolism, PI3K–AKT–mTOR and mTORC1-related biology, but no term is
BH-FDR significant at 0.05.

Interpretation: high-stress mimic compounds should not be treated as
mechanistic hits until this module is re-scored with curated generic
apoptosis/UPR/heat-shock/cell-cycle references.

## H2 — reversal_only__stable_2

Working hypothesis: this is the strongest residual/reversal module and may
represent glycolytic/hypoxic adaptation coupled to ECM remodeling at a
LAM-cell/LAF-niche interface, rather than a pure mTORC1 response.

Evidence: strongest cross-release match, Jaccard 0.717; baseline direction
fractions are approximately 0.83 and 0.80, with stress-adjusted retention
approximately 0.95 and 0.97. Leading genes include ENO2, GAPDH, PGK1, PKM,
TPI1, ASPN, DCN, EFEMP1, MFAP2, P4HA2 and PLOD3. In GSE49414, the RET-oriented
RPI-1 profile has the expected reversal direction for 46/83 analyzable module
genes. In GSE207322, ibrutinib has the expected reversal direction for 51/86
genes in WT TMD8 cells; C481 mutants are less consistent, but the attenuation
is not clean enough to assign causality to BTK.

Human mapping is compatible with this idea but remains preliminary: GSE135851
candidate groups trend higher for this module, while GSE302356 sample scores
are highest in LAM20 and LAF-niche marker fractions are visible in LAM3/LAM4.
These are not formal paired LAMCORE/LAF annotations.

Decisive next test: use selective RET and BTK pharmacology plus independent
genetic perturbation, and ask whether the same gene module changes. Only then
test sirolimus combinations.

## H3 — reversal_only__stable_3

Working hypothesis: a smaller, context-sensitive transcription/EMT or
state-switch branch may exist, but it is not stable enough for target
prioritization.

Evidence: Jaccard 0.519 and leading genes including HOXA5, LOXL2, RBBP4,
SDCBP and CSRP2. RPI-1 gives the expected reversal direction for 10/26
analyzable genes.

Limitations: direction fractions are low, retention ratios are unstable when
the baseline module effect is small, and the BTK reference does not show a
clean WT-to-mutant attenuation. Enrichment is not BH-FDR significant.

Decision rule: keep as an exploratory state-switch module; do not make it a
drug-selection gate without a third independent model or selective genetic
perturbation.

## Overall working conclusion

The most useful current hypothesis is H2: a reproducible residual program
linking metabolic adaptation and ECM/niche remodeling. H1 remains important as
a confounder-control problem because generic stress explains a substantial
fraction of the mimic landscape. H3 is a discovery lead, not yet a robust
mechanistic module. No module is promoted to a Tier 1 combination candidate at
this stage.
