# Fixed hydrogel translation-residual core: functional annotation

## Analysis object

The object was fixed to the 13 genes satisfying `residual_category=hydrogel_residual_q10`, baseline-effect eligibility for both RMC-6272 and eFT-508, and `both_drugs_distance_reduced=True`:

`CDC42EP3, RND3, SERPINE2, GPR27, WWTR1, FBN2, REEP2, ZNF354C, PNMA2, CACFD1, NFATC4, SPIN4, GPC4`

Plastic residual genes were not intersected into this object.

## Functional convergence

The 13 genes were annotated with GO Biological Process, GO Cellular Component, Reactome and MSigDB Hallmark. The repeated multi-gene themes were:

- ECM: SERPINE2 and FBN2;
- cell adhesion: RND3, SERPINE2 and GPC4;
- Rho GTPase: CDC42EP3 and RND3;
- TGF-β: WWTR1 and FBN2;
- migration: CDC42EP3 and SERPINE2.

Actin cytoskeleton, focal adhesion, Hippo/YAP/TAZ and mechanotransduction were present but currently had single-gene support in the operational keyword summary. WWTR1 provided the clearest Hippo/YAP/TAZ signal, while FBN2/SERPINE2 provided the clearest ECM/EMT-related convergence.

## Interpretation

The result suggests that the shared translation-residual core is more compatible with an ECM–adhesion–Rho/TGF-β–migration axis than with a generic translation-only label. This is a functional convergence hypothesis, not proof that all 13 genes act in one cell type or pathway. The four gene-set libraries overlap substantially, and the 13-gene size makes term-level FDR secondary to repeated gene-level support.

## Next question

Test whether the repeated themes remain present in human LAMCORE/LAF states and whether the 13-gene module, rather than isolated genes, is altered by selective mTORC1/MNK1/2 or genetic perturbation.
