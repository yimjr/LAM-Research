# Initial findings — 2026-08-22

## What the first run supports

- GSE179044 was analyzed as the complete WT/TSC2−/− × vehicle/rapamycin × plastic/hydrogel 2×2×2 design. All requested contrasts are present: TSC2-loss, rapamycin residual, genotype-dependent response, rapamycin-conditioned hydrogel-specific residual, and G×R×E environment-dependent response.
- The ratio is gated at `|d0| >= 0.5`; it is signed and descriptive. Moderated standard errors and FDR values use a pooled residual-variance prior because each cell has only two biological replicates.
- In hydrogel, the contrast classes were: 1,056 near-complete rescue, 3,099 partial rescue + residual, 668 persistent residual, 302 worsened residual, and 1,385 direction reversals among ratio-eligible genes. These are discovery strata, not independent evidence.
- The strongest rapamycin-conditioned hydrogel-vs-plastic residual effects included NNMT, FBLN5, COL8A1, MIR210HG and CAT. They are hypothesis-generating because the full environment interaction—not a significant/non-significant comparison—was used.
- Environment-dependent escape had fewer robust signals; exploratory examples included TVP23C, SPARCL1 and PRELP. This should remain below confirmed escape until supported by the external model and functional interpretation.
- GSE27982 provided the required independent mouse Tsc2×rapamycin 2×2 comparison. Top-effect sign concordance with the human factorial contrasts was 0.91 for TSC2-loss, 0.47 for residual after rapamycin, and 0.82 for genotype-dependent rapamycin response among the selected top-feature overlaps. This supports the loss phenotype more strongly than a universal residual program.
- GSE16944 reproduced its intended historical role: MMP2 was strongly elevated in TSC2-null AML-like cells and remained high after rapamycin. Because the dataset lacks TSC2-restored + rapamycin, this is not formal residual or G×R replication.

## Human-state mapping status

The copied GSE135851 snapshot yielded preliminary patient-by-state rank scores for `candidate` versus `other` labels across five donors. It does not contain the formal LAMCORE1/2/3 and LAF-seed/LAF-niche taxonomy, so it is not used as the formal GSE302356 validation. A separate processed-input contract is recorded for GSE302356 so that RNA, ATAC/gene activity and spatial evidence remain modality-aware.

## Research ideas to test

1. The hydrogel-specific residual signal suggests a possible ECM–metabolic coupling rather than a generic “rapamycin resistance” state. NNMT, FBLN5 and COL8A1 are a compact starting panel for checking whether matrix context preserves a metabolic/ECM residual after mTORC1 suppression.
2. The lower residual sign concordance between GSE179044 and low-serum GSE27982 cautions against treating every persistent program as a conserved escape mechanism. Candidate escape should require direction, baseline state, function and model support together.
3. Direction reversals are common enough to deserve a separate mechanism audit, but they are deliberately excluded from default CMap queries until over-correction, dose/time effects, low-baseline instability and true reversal are separated.

## Additional interpretation after contrast review

- The residual is not simply an unblocked mTORC1 program: the curated mTORC1 module had near-zero mean hydrogel residual effect and no FDR-supported member, while selected myogenic, ECM/invasion, metabolic and autophagy genes retained signal. This is a derived interpretation from the same discovery experiment, not independent validation.
- Hydrogel-specific residual is a selective state reshaping rather than a universal increase in resistance. NNMT, COL8A1, MIR210HG and SLC40A1 were more positive in the hydrogel-conditioned residual, while FBLN5, DCN and LUM moved toward stronger rescue relative to plastic. Direction must be interpreted together with the sign of d0; a positive G×E simple interaction does not always mean a stronger residual.
- The strongest hydrogel G×R×E signal is sparse: only eight genes passed the current effect/FDR screen, and most had q approximately 0.097. At present this argues against a broad environment-dependent escape program and favors a small, gene-specific niche response.
- NNMT and COL8A1 are the clearest current candidate compensatory-escape hypotheses in hydrogel because they have positive baseline TSC2-loss, larger positive rapamycin residual and positive genotype-dependent response. They are not yet externally replicated in GSE27982 and should not be called confirmed escape.
- MMP2 now has a useful three-layer interpretation: strong historical rapamycin-insensitive support in GSE16944, partial residual in GSE179044 hydrogel, and a descriptive residual increase in GSE27982 despite a baseline effect below the ratio gate. This makes MMP2 a good orthogonal assay control, not a Tier 1 drug target.

## New input pass

- GSE84476 is no longer only a transcript-level file: a GENCODE v24 map was added and gene-level log2 TPM contrasts were generated. STAT3 itself decreases after siSTAT3 in both 102-cell and 103-cell contexts, confirming target engagement; rapamycin responses are context-dependent rather than uniformly STAT3-like. COL8A1, DCN, MMP2 and stress genes show context-specific differences, but the dataset has limited replication and remains mechanism support.
- GSE104335 is now locally complete at the archive level. GEO metadata identify three biological replicates each for shGFP+DMSO, shGFP+100 nM rapamycin and shSRPK2+DMSO in the LAM 621-101 cell line. Gene-level CHP/CEL extraction is still pending, so no SRPK2 conclusion is claimed yet.
- GSE302356 raw matrices now support a preliminary human tissue check: LAM3/LAM4 scRNA-seq, LAM18 Visium HD and LAM20 Visium were scored against the residual and functional programs. LAM20 had the highest raw ECM, myogenic, metabolism, stress and autophagy scores; LAM3/LAM4 were more stress-weighted; LAM18 HD was lower overall after tissue-spot filtering. This is an exploratory sample/modality pattern, not formal LAMCORE state validation, because the downloaded archives lack the state labels and specimens/modalities are not paired.

## State-panel refinement

- The accessible article text makes an important limitation explicit: LAMCORE3 shares LAMCORE1 biology but lacks a unique marker/differential signature and appears to be a lower-transcription state. Therefore it is not scientifically appropriate to manufacture a conventional LAMCORE3 marker list. The project now represents it with shared LAM-core genes plus a clearly labeled translation/low-activity surrogate, and keeps it below formal state-validation claims.
- The paper-derived panel records the reported subtype markers (LAMCORE1: PCP4, PI15, LRRC7, PLIN4, GAD1; LAMCORE2: MMP11, TDO2, MLANA, SPINK13), shared LAM markers, canonical LAF genes (FAP, S100A4, VIM, IGFBP7, SPARC), and the reported LAF-seed/LAF-niche functional programs. This is a useful bridge for exploratory scoring until Supplementary Table 1 or author-processed state metadata are available.
- A practical implication is that human mapping should ask two separate questions: whether the residual program co-occurs with a LAM state panel, and whether it co-occurs with an activated fibroblast/niche program. A high ECM score alone should not be labeled LAMCORE2, because the paper describes LAMCORE2 and LAF-niche as distinct but overlapping ECM-rich compartments.

## SRPK2 mechanism pass

- GSE104335 was resolved through its processed `HTA-2_0` `sst-rma-gene-full` CHP files: 30,262 transcript clusters mapped to 23,910 gene symbols across 9 samples (three biological replicates in each of shGFP+vehicle, shGFP+rapamycin, and shSRPK2+vehicle).
- In the shGFP background, rapamycin increased COL8A1 (+1.16), NNMT (+1.23), DCN (+0.58), ACTA2 (+0.69), and MMP2 (+0.56), while strongly reducing HMGCS1, NUPR1, and HSPA5. These changes are notable because NNMT and COL8A1 are also among the clearest GSE179044 hydrogel residual/hydrogel-specific residual hypotheses.
- SRPK2 knockdown produced a different state: SRPK2 itself fell by 1.52 log2 units, COL3A1, LUM, FBLN5 and NNMT decreased, while DDIT3 and ATF4 increased. NNMT and COL8A1 moved opposite to the rapamycin response, so SRPK2 perturbation is not a simple transcriptional surrogate for rapamycin.
- New hypothesis: SRPK2 may be a regulator of an ECM–metabolic residual axis that remains after mTORC1 inhibition, but its therapeutic role is unresolved. The decisive missing experiment is not another bulk signature; it is the interaction `shSRPK2 × rapamycin` in the same LAM model, followed by viability/ECM readouts.
- This is a useful example of why mechanism support and formal escape evidence must remain separate: GSE104335 strengthens the NNMT/COL8A1 mechanism hypothesis, but cannot by itself upgrade it to confirmed escape or a Tier 1 combination.
- A direct overlap screen found 19 genes significant in both the GSE104335 shGFP rapamycin contrast and the GSE179044 rapamycin-conditioned hydrogel-specific residual; 11 were directionally concordant and 8 discordant. The concordant NNMT/COL8A1 pair is especially interesting because SRPK2 knockdown suppresses both, whereas rapamycin increases them in the shGFP background. This suggests a testable SRPK2-sensitive residual axis, while the discordant FAP/SLC40A1 examples prevent a generic “ECM resistance” interpretation.
