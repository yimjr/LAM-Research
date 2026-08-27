# GSE277844 hydrogel translation-residual core functional annotation

## Fixed analysis object

The analysis object is fixed by `residual_category=hydrogel_residual_q10`, baseline-effect eligibility for both RMC-6272 and eFT-508, and `both_drugs_distance_reduced=True`. It contains exactly 13 genes; ordinary/plastic residual genes are not used in this functional annotation step.

CDC42EP3; RND3; SERPINE2; GPR27; WWTR1; FBN2; REEP2; ZNF354C; PNMA2; CACFD1; NFATC4; SPIN4; GPC4

## Data and method

- concordance input: `data/processed/translation_analysis/GSE277844_translation_residual_drug_gene_concordance.csv`
- translation input/background: `data/processed/translation_analysis/GSE277844_tsc2_loss_translation_effects.csv`; background = 13985 genes with finite conditional translation effect
- libraries: GO Biological Process, GO Cellular Component, Reactome and MSigDB Hallmark from local GMT files
- term-level enrichment uses a hypergeometric over-representation test with the analyzable GSE277844 background and Benjamini–Hochberg FDR within each library.
- the primary interpretation is gene-level annotation and repeated themes across the 13 genes; non-significant FDR does not erase a coherent small-module signal.

## Per-gene functional annotation summary

    gene translation_direction  translation_effect  translation_q  GO_Biological_Process_n_terms  GO_Cellular_Component_n_terms  Reactome_n_terms  MSigDB_Hallmark_n_terms                                                focus_themes
CDC42EP3      translation_down           -2.149265       0.030456                            8.0                            2.0               8.0                      0.0                     Rho_GTPase;actin_cytoskeleton;migration
    RND3      translation_down           -1.589212       0.125281                            0.0                            2.0               5.0                      1.0 Rho_GTPase;cell_adhesion;focal_adhesion;mechanotransduction
SERPINE2      translation_down           -1.997215       0.071526                           25.0                            6.0               5.0                      1.0                cell_adhesion;extracellular_matrix;migration
   GPR27      translation_down           -3.553566       0.069683                            0.0                            0.0               9.0                      0.0                                                            
   WWTR1      translation_down           -1.293564       0.083127                           34.0                            2.0              22.0                      2.0                                      Hippo_YAP_TAZ;TGF_beta
    FBN2      translation_down           -1.758873       0.109605                           16.0                            3.0               2.0                      1.0                               TGF_beta;extracellular_matrix
   REEP2        translation_up            1.256035       0.080923                            2.0                            2.0               0.0                      0.0                                                            
 ZNF354C      translation_down          -11.882506       0.055213                            0.0                            0.0               3.0                      0.0                                                            
   PNMA2      translation_down           -0.959546       0.111693                            3.0                            0.0               0.0                      0.0                                                            
  CACFD1        translation_up            1.377495       0.109605                            1.0                            0.0               0.0                      0.0                                                            
  NFATC4        translation_up            1.117863       0.121011                           19.0                            2.0               0.0                      0.0                                                            
   SPIN4      translation_down           -1.294397       0.142389                            3.0                            0.0               0.0                      0.0                                                            
    GPC4      translation_down           -1.043283       0.143891                            8.0                            7.0              25.0                      2.0                                               cell_adhesion

## Repeated functional themes

               theme  n_core_genes              genes  n_supporting_libraries                           supporting_libraries
extracellular_matrix             2      SERPINE2;FBN2                       2                 GO_Cellular_Component;Reactome
       cell_adhesion             3 RND3;SERPINE2;GPC4                       2    GO_Biological_Process;GO_Cellular_Component
          Rho_GTPase             2      CDC42EP3;RND3                       1                                       Reactome
            TGF_beta             2         WWTR1;FBN2                       3 GO_Biological_Process;MSigDB_Hallmark;Reactome
           migration             2  CDC42EP3;SERPINE2                       1                          GO_Biological_Process

Themes are treated as repeated when at least two of the 13 genes have one or more supporting terms. Keyword-based theme labels are operational summaries; the full term-level annotations are retained in the long output.

## Term enrichment (descriptive)

              library                                                       term_name  overlap_count overlap_genes         focus_themes  p_value      fdr
GO_Biological_Process               Positive Regulation Of Osteoblast Differentiation              2    FBN2;WWTR1                      0.000494 0.054796
GO_Biological_Process                        Regulation Of Osteoblast Differentiation              2    FBN2;WWTR1                      0.001912 0.057217
GO_Biological_Process         Maintenance Of Protein Location In Extracellular Region              1          FBN2                      0.003713 0.057217
GO_Biological_Process                   Negative Regulation Of Plasminogen Activation              1      SERPINE2                      0.004640 0.057217
GO_Biological_Process              Sequestering Of Extracellular Ligand From Receptor              1          FBN2                      0.004640 0.057217
GO_Biological_Process                       Negative Regulation Of Protein Maturation              1      SERPINE2                      0.005565 0.057217
GO_Biological_Process                                     Embryonic Eye Morphogenesis              1          FBN2                      0.006490 0.057217
GO_Biological_Process                     Negative Regulation Of Platelet Aggregation              1      SERPINE2                      0.006490 0.057217
GO_Biological_Process                Positive Regulation Of Astrocyte Differentiation              1      SERPINE2                      0.006490 0.057217
GO_Biological_Process                               Negative Regulation Of Hemostasis              1      SERPINE2                      0.007414 0.057217
GO_Cellular_Component                                                     Microfibril              1          FBN2 extracellular_matrix 0.009260 0.111871
GO_Cellular_Component                        Collagen-Containing Extracellular Matrix              2 FBN2;SERPINE2 extracellular_matrix 0.022296 0.111871
GO_Cellular_Component                           Endoplasmic Reticulum Tubular Network              1         REEP2                      0.023911 0.111871
GO_Cellular_Component                                            Supramolecular Fiber              1          FBN2                      0.023911 0.111871
GO_Cellular_Component                                          Neuromuscular Junction              1      SERPINE2        cell_adhesion 0.026636 0.111871
GO_Cellular_Component                        Extracellular Membrane-Bounded Organelle              1      SERPINE2                      0.043730 0.134510
GO_Cellular_Component                                           Extracellular Vesicle              1      SERPINE2                      0.045513 0.134510
GO_Cellular_Component                                                     Golgi Lumen              1          GPC4                      0.052616 0.134510
GO_Cellular_Component                                          Platelet Alpha Granule              1      SERPINE2                      0.060547 0.134510
GO_Cellular_Component                                                 Lysosomal Lumen              1          GPC4                      0.064052 0.134510
      MSigDB_Hallmark                               Epithelial Mesenchymal Transition              2 FBN2;SERPINE2                      0.011706 0.070238
      MSigDB_Hallmark                                              TGF-beta Signaling              1         WWTR1             TGF_beta 0.049071 0.147212
      MSigDB_Hallmark                                                  UV Response Dn              1          RND3                      0.123476 0.156643
      MSigDB_Hallmark                                                      Myogenesis              1         WWTR1                      0.140612 0.156643
      MSigDB_Hallmark                                                         Hypoxia              1          GPC4                      0.153459 0.156643
      MSigDB_Hallmark                                                      Glycolysis              1          GPC4                      0.156643 0.156643
             Reactome       RUNX3 Regulates YAP1-mediated Transcription R-HSA-8951671              1         WWTR1        Hippo_YAP_TAZ 0.006490 0.087215
             Reactome                             Physiological Factors R-HSA-5578768              1         WWTR1                      0.007414 0.087215
             Reactome                          Dissolution Of Fibrin Clot R-HSA-75205              1      SERPINE2                      0.008337 0.087215
             Reactome            Common Pathway Of Fibrin Clot Formation R-HSA-140875              1      SERPINE2                      0.009260 0.087215
             Reactome  YAP1- And WWTR1 (TAZ)-stimulated Gene Expression R-HSA-2032785              1         WWTR1        Hippo_YAP_TAZ 0.011102 0.087215
             Reactome Defective EXT1 Causes Exostoses 1, TRPS2 And CHDS R-HSA-3656253              1          GPC4                      0.012022 0.087215
             Reactome         Intrinsic Pathway Of Fibrin Clot Formation R-HSA-140837              1      SERPINE2                      0.012942 0.087215
             Reactome        Defective B3GALT6 Causes EDSP2 And SEMDJL1 R-HSA-4420332              1          GPC4                      0.017526 0.087215
             Reactome                   Defective B3GAT3 Causes JDSSDHD R-HSA-3560801              1          GPC4                      0.017526 0.087215
             Reactome      Defective B4GALT7 Causes EDS, Progeroid Type R-HSA-3560783              1          GPC4                      0.017526 0.087215

## Interpretation

The repeated-theme table should be read as a functional convergence check, not as proof that every annotated term is active in the same cell. GO/Reactome/Hallmark terms are overlapping, so counts across libraries are not independent evidence. The 13 genes were selected using residual and treatment-response criteria, and GSE277844 remains a cross-model human NPC dataset rather than a direct LAM experiment.

## Outputs

- `GSE277844_hydrogel_translation_core_genes.csv`: the locked 13-gene analysis object with treatment-response fields.
- `GSE277844_hydrogel_translation_core_functional_annotations.csv.gz`: gene × library × term annotations.
- `GSE277844_hydrogel_translation_core_gene_function_summary.csv`: one row per fixed gene with library term counts and focus themes.
- `GSE277844_hydrogel_translation_core_functional_theme_by_gene.csv`: gene-level theme support matrix.
- `GSE277844_hydrogel_translation_core_functional_theme_summary.csv`: repeated-theme summary across genes.
- `GSE277844_hydrogel_translation_core_functional_enrichment.csv.gz`: descriptive term enrichment with FDR.
