# GSE277844 translation-loss and residual-program analysis

## Research question

先定义 TSC2 loss 是否改变 polysome-associated mRNA 相对于 total mRNA 的翻译效率，再与 GSE179044 的 persistent/hydrogel residual 比较，最后检查 mTORC1/MNK1/2 translation-targeting treatment 是否把重叠基因的 KO-vs-WT 翻译效率差异拉回 WT 附近。

## Data and method

- input: `data/raw/GSE277844/GSE277844_raw_counts.txt.gz`; GEO supplementary source: https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277844/suppl/GSE277844_raw_counts.txt.gz
- raw genes: 19354; genes after count filter: 13985; TE pairs: 16
- TE = log2(CPM polysome + 0.5) − log2(CPM total + 0.5), paired by genotype × treatment × replicate.
- TSC2-loss translation genes are called from the conditional polysome-vs-total model at the published-style FDR cutoff; the translation effect size is retained for ranking and interpretation rather than used as the primary significance test.
- residual overlap uses the existing moderated GSE179044 contrasts. `persistent_residual_plastic` and `persistent_residual_hydrogel` are kept as separate categories; `_q10` additionally requires absolute residual effect ≥ 0.5 and moderated FDR ≤ 0.10. `hydrogel_residual_q10` is the post-rapamycin hydrogel residual with the same effect/FDR gate. No intersection between ordinary and hydrogel residual is required.
- drug rescue is assessed by signed residual ratio = post-treatment KO-vs-WT TE effect / DMSO KO-vs-WT TE effect. Ratios are descriptive, especially for RMC-6272 where each group has two replicates.

## Step 1: TSC2-loss translation program

- selected translation-up genes: 89
- selected translation-down genes: 107

## Step 2: overlap with GSE179044 residual programs

translation_direction                residual_category  n_translation_genes  n_residual_genes  n_overlap  overlap_fraction_of_translation_genes
       translation_up      persistent_residual_plastic                   89               409          2                               0.022472
       translation_up  persistent_residual_plastic_q10                   89               154          1                               0.011236
       translation_up     persistent_residual_hydrogel                   89               668          4                               0.044944
       translation_up persistent_residual_hydrogel_q10                   89               381          3                               0.033708
       translation_up            hydrogel_residual_q10                   89              1681          7                               0.078652
       translation_up   hydrogel_specific_residual_q10                   89               182          0                               0.000000
     translation_down      persistent_residual_plastic                  107               409          3                               0.028037
     translation_down  persistent_residual_plastic_q10                  107               154          0                               0.000000
     translation_down     persistent_residual_hydrogel                  107               668          6                               0.056075
     translation_down persistent_residual_hydrogel_q10                  107               381          4                               0.037383
     translation_down            hydrogel_residual_q10                  107              1681         13                               0.121495
     translation_down   hydrogel_specific_residual_q10                  107               182          0                               0.000000

The union of selected translation genes overlapping at least one residual category contains 26 genes for recordkeeping. Treatment summaries are reported separately by residual category, so this union is not used to require ordinary and hydrogel residual to overlap first.

## Step 3: translation-targeting treatment

The following summaries concern genes in each translation × residual overlap category separately. `distance_to_WT_reduced` asks whether the KO-vs-WT TE distance became smaller; it does not by itself prove a mechanistic rescue.

   drug          drug_description                residual_category  n_translation_residual_overlap_genes  n_baseline_effect_eligible  n_distance_to_WT_reduced  fraction_distance_to_WT_reduced  median_abs_baseline_effect  median_abs_after_effect  median_signed_residual_ratio  n_near_complete_rescue  n_partial_rescue_residual  n_persistent_residual  n_worsened_residual  n_direction_reversal  n_after_significant_and_effect_sized
RMC6272 mTORC1 inhibitor RMC-6272            hydrogel_residual_q10                                    20                          18                        15                         0.833333                    1.335946                 0.817711                      0.520135                       1                         12                      1                    3                     1                                     3
RMC6272 mTORC1 inhibitor RMC-6272     persistent_residual_hydrogel                                    10                           9                         8                         0.888889                    1.126790                 0.603099                      0.456994                       1                          6                      1                    1                     0                                     1
RMC6272 mTORC1 inhibitor RMC-6272 persistent_residual_hydrogel_q10                                     7                           7                         6                         0.857143                    1.256035                 0.907251                      0.480161                       0                          5                      1                    1                     0                                     1
RMC6272 mTORC1 inhibitor RMC-6272      persistent_residual_plastic                                     5                           5                         4                         0.800000                    3.553566                 4.164667                      0.679198                       1                          3                      0                    1                     0                                     1
RMC6272 mTORC1 inhibitor RMC-6272  persistent_residual_plastic_q10                                     1                           1                         1                         1.000000                    6.815757                 4.629250                      0.679198                       0                          1                      0                    0                     0                                     0
 eFT508  MNK1/2 inhibitor eFT-508            hydrogel_residual_q10                                    20                          18                        15                         0.833333                    1.335946                 0.680374                      0.503256                       1                          9                      4                    1                     3                                     0
 eFT508  MNK1/2 inhibitor eFT-508     persistent_residual_hydrogel                                    10                           9                         9                         1.000000                    1.126790                 0.452580                      0.404862                       1                          6                      0                    0                     2                                     0
 eFT508  MNK1/2 inhibitor eFT-508 persistent_residual_hydrogel_q10                                     7                           7                         7                         1.000000                    1.256035                 0.452580                      0.404862                       0                          5                      0                    0                     2                                     0
 eFT508  MNK1/2 inhibitor eFT-508      persistent_residual_plastic                                     5                           5                         2                         0.400000                    3.553566                 2.860342                      0.729964                       0                          2                      2                    0                     1                                     0
 eFT508  MNK1/2 inhibitor eFT-508  persistent_residual_plastic_q10                                     1                           1                         0                         0.000000                    6.815757                 7.260806                      1.065297                       0                          0                      1                    0                     0                                     0

## Background comparison

The direct comparison uses all selected GSE277844 translation-abnormal genes; the effect-eligible comparison applies the same |baseline translation effect| ≥ 0.5 gate used for rescue ratios. The all-gene comparison includes overlap genes, while the non-overlap comparison and Fisher exact p-value are sensitivity diagnostics, not gene-independent proof.

   drug                residual_category  n_all_selected_translation_genes  n_overlap_genes_all_selected  n_overlap_recovered_all_selected  all_selected_recovery_rate  overlap_recovery_rate_all_selected  overlap_minus_all_selected_rate  n_all_selected_translation_genes_effect_eligible  n_overlap_genes_effect_eligible  n_overlap_recovered  all_recovery_rate  overlap_recovery_rate  overlap_minus_all_rate  overlap_minus_nonoverlap_rate  fisher_p_overlap_vs_nonoverlap
RMC6272            hydrogel_residual_q10                               196                            20                                15                    0.678571                            0.750000                         0.071429                                               176                               18                   15           0.727273               0.833333                0.106061                       0.118143                        0.405215
RMC6272     persistent_residual_hydrogel                               196                            10                                 8                    0.678571                            0.800000                         0.121429                                               176                                9                    8           0.727273               0.888889                0.161616                       0.170326                        0.447396
RMC6272 persistent_residual_hydrogel_q10                               196                             7                                 6                    0.678571                            0.857143                         0.178571                                               176                                7                    6           0.727273               0.857143                0.129870                       0.135249                        0.675545
RMC6272      persistent_residual_plastic                               196                             5                                 4                    0.678571                            0.800000                         0.121429                                               176                                5                    4           0.727273               0.800000                0.072727                       0.074854                        1.000000
RMC6272  persistent_residual_plastic_q10                               196                             1                                 1                    0.678571                            1.000000                         0.321429                                               176                                1                    1           0.727273               1.000000                0.272727                       0.274286                        1.000000
 eFT508            hydrogel_residual_q10                               196                            20                                15                    0.750000                            0.750000                         0.000000                                               176                               18                   15           0.778409               0.833333                0.054924                       0.061181                        0.766544
 eFT508     persistent_residual_hydrogel                               196                            10                                 9                    0.750000                            0.900000                         0.150000                                               176                                9                    9           0.778409               1.000000                0.221591                       0.233533                        0.209965
 eFT508 persistent_residual_hydrogel_q10                               196                             7                                 7                    0.750000                            1.000000                         0.250000                                               176                                7                    7           0.778409               1.000000                0.221591                       0.230769                        0.350577
 eFT508      persistent_residual_plastic                               196                             5                                 2                    0.750000                            0.400000                        -0.350000                                               176                                5                    2           0.778409               0.400000               -0.378409                      -0.389474                        0.072980
 eFT508  persistent_residual_plastic_q10                               196                             1                                 0                    0.750000                            0.000000                        -0.750000                                               176                                1                    0           0.778409               0.000000               -0.778409                      -0.782857                        0.221591

## Gene-level agreement between translation-targeting drugs

`both_reduced` means both drugs reduced the KO-vs-WT translation distance for that gene. `RMC6272_only_reduced` and `eFT508_only_reduced` identify drug-specific support; genes failing the baseline effect gate are not assigned a recovery-support pattern.

               residual_category  n_overlap_genes  n_both_drug_effect_eligible  n_both_reduced  fraction_both_reduced  n_RMC6272_only_reduced  n_eFT508_only_reduced  n_neither_reduced  fraction_any_drug_reduced  reduced_set_jaccard                                                                   both_reduced_genes RMC6272_only_reduced_genes eFT508_only_reduced_genes neither_reduced_genes
           hydrogel_residual_q10               20                           18              13               0.722222                       2                      2                  1                   0.944444             0.764706 CACFD1;CDC42EP3;FBN2;GPC4;GPR27;NFATC4;PNMA2;REEP2;RND3;SERPINE2;SPIN4;WWTR1;ZNF354C                FIBIN;HOXC6                APOL1;TYMS                 ZWINT
    persistent_residual_hydrogel               10                            9               8               0.888889                       0                      1                  0                   1.000000             0.888889                                 FBN2;NFATC4;PNMA2;REEP2;SERPINE2;UGP2;ZNF354C;ZNF563                                                APOL1                      
persistent_residual_hydrogel_q10                7                            7               6               0.857143                       0                      1                  0                   1.000000             0.857143                                             FBN2;NFATC4;PNMA2;REEP2;SERPINE2;ZNF354C                                                APOL1                      
     persistent_residual_plastic                5                            5               2               0.400000                       2                      0                  1                   0.800000             0.500000                                                                         GPR27;RNF182                HOXC6;NETO1                                           ALCAM
 persistent_residual_plastic_q10                1                            1               0               0.000000                       1                      0                  0                   1.000000             0.000000                                                                                                           HOXC6                                                

## Interpretation limits

GSE277844 is a human NPC model and is biologically distinct from the LAM cell model in GSE179044. The comparison is therefore a cross-model program test, not a direct LAM replication. Translation efficiency is estimated from normalized bulk count fractions, and treatment-specific rescue ratios are unstable when the baseline TE effect is small. RMC-6272/eFT-508 results should be followed by gene/module-level confirmation and, where possible, independent translational or genetic perturbation data.

## Outputs

- `GSE277844_tsc2_loss_translation_effects.csv`: all filtered genes with TE effect, p/q values and up/down selection.
- `GSE277844_translation_residual_overlap_genes.csv`: gene-level overlap with persistent/hydrogel residual categories and direction relation.
- `GSE277844_translation_residual_overlap_drug_effects.csv`: treatment effects for each gene × residual category overlap, including rescue ratios and classes.
- `GSE277844_translation_residual_overlap_drug_summary.csv`: treatment × residual-category summary of recovery toward WT.
- `GSE277844_translation_residual_recovery_background_comparison.csv`: overlap recovery compared with all selected translation-abnormal genes and the non-overlap sensitivity background.
- `GSE277844_translation_residual_drug_gene_concordance.csv`: gene-level RMC-6272/eFT-508 recovery pattern for every residual-category overlap.
- `GSE277844_translation_residual_drug_concordance_summary.csv`: per-category counts and gene lists for both-drug, RMC-6272-only and eFT-508-only recovery.
