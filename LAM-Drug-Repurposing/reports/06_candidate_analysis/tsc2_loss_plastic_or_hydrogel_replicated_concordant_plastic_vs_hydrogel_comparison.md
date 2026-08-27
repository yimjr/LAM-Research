# Plastic vs hydrogel LINCS signature comparison

同一批 66 个候选药物分别使用 `tsc2_loss_plastic` 与 `tsc2_loss_hydrogel` top150 up + top150 down 面板；两个 LINCS release 先独立汇总，再比较共同基因。
这里的主要问题不是重新比较同一套 LINCS perturbation effect，而是比较两个 disease panel 各自判定出的 reversal 基因集合：哪些基因在两个 panel 都被 reversal、哪些只在 plastic 或 hydrogel panel 中出现。
两边使用同一批 LINCS drug signatures，因此共同可比较基因上的 drug effect 数值本来就相同；effect correlation 只作为结构性可比性诊断，不作为生物学证据。

- plastic panel: 300 genes; hydrogel panel: 300 genes
- same-direction panel overlap: 165; panel union: 435
- median common drug-effect genes compared per release: {'GSE70138': 113.0, 'GSE92742': 113.0}

## Primary: reversal-set comparison

对每个 drug × release 定义 plastic reversal set 与 hydrogel reversal set，并报告集合大小、交集、Jaccard、overlap coefficient 及两侧特异基因。跨 release 的主表使用这些集合指标的中位数；不把同一 LINCS effect 的相关性当成复现。

 n_drug_comparison_rows  n_unique_drugs  n_gene_comparison_rows  n_datasets  median_status_agreement  median_status_agreement_both_panels  median_reversal_jaccard_release_level  median_reversal_overlap_coefficient_release_level  median_n_reversal_both_release_level  median_effect_spearman_secondary  n_direction_switch_rows
                    132              66                   57420           2                  0.57931                                  1.0                                0.37356                                           0.575758                                  33.5                               1.0                        0

### Drugs with the largest shared reversal sets

                  entity_id        pert_iname  n_release_rows  n_releases_with_jaccard_ge_0_5  median_n_reversal_plastic  median_n_reversal_hydrogel  median_n_reversal_both  median_reversal_jaccard  median_reversal_overlap_coefficient  median_n_plastic_only_reversal  median_n_hydrogel_only_reversal
compound::homoharringtonine homoharringtonine               2                               0                       92.5                       110.5                    54.5                 0.368195                             0.590204                            38.0                             56.0
     compound::lestaurtinib      lestaurtinib               2                               0                       90.5                        84.5                    47.5                 0.374475                             0.565934                            43.0                             37.0
         compound::zstk-474          ZSTK-474               2                               0                       82.5                        78.0                    47.0                 0.414749                             0.602961                            35.5                             31.0
      compound::gsk-2126458       GSK-2126458               2                               0                       91.0                        90.0                    47.0                 0.350936                             0.522410                            44.0                             43.0
       compound::wye-125132        WYE-125132               2                               0                       95.0                        89.5                    47.0                 0.341667                             0.524913                            48.0                             42.5
          compound::osi-027           OSI-027               2                               0                       81.0                        83.0                    45.5                 0.384748                             0.570383                            35.5                             37.5
           compound::pi-103            PI-103               2                               0                       84.5                        79.0                    45.0                 0.379701                             0.569551                            39.5                             34.0
         compound::gdc-0941          GDC-0941               2                               0                       76.5                        71.5                    43.5                 0.416422                             0.608470                            33.0                             28.0
        compound::fadrozole         fadrozole               2                               0                       77.5                        85.0                    43.5                 0.366750                             0.558529                            34.0                             41.5
       compound::nvp-bez235        NVP-BEZ235               2                               0                       85.5                        80.0                    43.0                 0.350597                             0.537238                            42.5                             37.0
     compound::levalbuterol      levalbuterol               2                               0                       69.0                        72.5                    42.5                 0.444141                             0.633462                            26.5                             30.0
         compound::azd-8055          AZD-8055               2                               0                       86.0                        80.0                    41.5                 0.332250                             0.518223                            44.5                             38.5
         compound::ql-x-138          QL-X-138               2                               0                       67.0                        81.5                    41.0                 0.381210                             0.612341                            26.0                             40.5
         compound::l-655240          L-655240               2                               0                       66.0                        69.5                    40.0                 0.411996                             0.606894                            26.0                             29.5
      compound::gsk-1059615       GSK-1059615               2                               0                       72.5                        75.5                    39.5                 0.364780                             0.551957                            33.0                             36.0
          compound::torin-1           torin-1               2                               0                       79.0                        78.0                    39.5                 0.335522                             0.505438                            39.5                             38.5
       compound::gossypetin        Gossypetin               2                               0                       69.5                        71.0                    39.0                 0.388086                             0.572251                            30.5                             32.0
       compound::novobiocin        novobiocin               2                               0                       73.0                        71.5                    39.0                 0.370676                             0.554558                            34.0                             32.5
     compound::ellagic-acid      ellagic-acid               2                               0                       66.0                        74.5                    38.5                 0.363124                             0.582027                            27.5                             36.0
     compound::jnj-26481585      JNJ-26481585               2                               0                       68.0                        65.5                    38.0                 0.381949                             0.577160                            30.0                             27.5

### Release-level shared reversal genes

                  entity_id        pert_iname  dataset  n_reversal_plastic  n_reversal_hydrogel  n_reversal_both  n_plastic_only_reversal  n_hydrogel_only_reversal  reversal_jaccard  reversal_overlap_coefficient                                                                                                                                                                                                                                                                                                                                                             shared_reversal_genes
        compound::fadrozole         fadrozole GSE92742                 106                  119               60                       46                        59          0.363636                      0.566038 ARID5B;ATP6V1C1;ATXN7L3B;AVPI1;CAT;CDKN1A;CLIP4;CRMP1;CSRP2;DDIT3;DIDO1;DUSP3;EDNRA;EFEMP1;FAM162A;FILIP1L;GDF15;GPC3;HBP1;HEXB;HSD17B14;HSPA9;IFNGR2;IGFBP5;ITGB5;LAMP2;LGALS3;LONP1;LUM;MBIP;MGST3;MTCH2;NEU1;NR3C1;NREP;PDE10A;PGK1;PKM;PLA2G15;PLA2G4A;PLOD3;PPT1;RSRP1;SDCBP;SRSF6;TAF9;TAF9B;TBC1D2;TNFRSF21;TOMM34;TOX;TPBG;TPI1;TSPAN31;TXN;UCHL1;VEGFB;WNT11;WSB2;ZNF711
compound::homoharringtonine homoharringtonine GSE92742                  99                  119               57                       42                        62          0.354037                      0.575758                           ARID5B;ARMC9;CABP1;CIRBP;CKB;CNTLN;CSRP2;EDNRA;EGLN3;ENO2;ERO1A;FAM162A;FILIP1L;GAPDH;GBE1;GOLPH3;HBP1;HEXB;HSD17B14;HSPA9;IFNGR2;ITGB5;LAMP2;LAMP3;LIMA1;LONP1;LOXL2;MBIP;MGST3;MTCH2;NAP1L3;NDUFA4L2;NR3C1;ODC1;P4HA2;PDE10A;PGK1;PKM;PLA2G15;PLOD3;PPT1;REXO2;SCARB1;SEMA3A;SPAG4;SRSF6;TOX;TPBG;TPI1;TRIB3;TXN;UCHL1;VEGFB;WNT11;WSB2;ZNF302;ZNF711
     compound::ellagic-acid      ellagic-acid GSE92742                  94                   96               55                       39                        41          0.407407                      0.585106                                            ATP6V1C1;CABP1;CDKN1A;CHI3L1;CKB;CNTLN;CRMP1;CSRP2;DDIT3;DUSP3;EDNRA;EDNRB;EFEMP1;ENO2;ERO1A;GBE1;GDF15;HBP1;HEXB;ITGB5;LAMB3;LAMP2;LGALS3;LONP1;LOXL2;MAGED2;MBIP;MOSPD1;NEU1;NREP;NXPH4;ODC1;P4HA2;PGK1;PKM;POPDC3;PRRX2;REXO2;SCARB1;SLC2A1;SLC38A6;STEAP3;SV2A;TAF9B;TBC1D2;TOX;TPBG;TPI1;TRIB3;TXN;VEGFA;WNT11;WSB2;ZNF302;ZNF711
     compound::levalbuterol      levalbuterol GSE92742                  94                   97               55                       39                        42          0.404412                      0.585106                                ARID5B;ARMC9;ATXN7L3B;CA9;CAT;CHI3L1;CLIP4;CNTLN;CSRP2;DDIT3;DUSP3;EDNRA;EFEMP1;EGLN3;ENO2;FAM162A;FBLN5;GBE1;GDF15;GPM6B;HBP1;HEXB;HSD17B14;HSPA9;IGFBP5;ITGB5;LAMP2;LONP1;LRRC15;LUM;MGST3;MTCH2;NR3C1;NREP;NXPH4;PFKP;PGK1;PKM;PLA2G15;PLA2G4A;PLOD3;REXO2;SCARB1;SLC2A1;SRSF6;SV2A;TAF9B;TNFRSF11B;TNFRSF21;TOMM34;TOX;TSPAN31;TXN;VEGFB;WNT11
     compound::jnj-26481585      JNJ-26481585 GSE70138                  81                   85               53                       28                        32          0.469027                      0.654321                                            ARID5B;ATXN7L3B;CAT;CCND2;CRMP1;CSRP2;DIDO1;DUSP3;EDNRA;EDNRB;EFEMP1;ERO1A;FBLN5;FILIP1L;GBE1;GOLPH3;GPC3;HBP1;HSPA9;IGFBP5;ITGB5;LAMB3;LIMA1;LONP1;LOXL2;MTCH2;NAP1L3;PFKP;PKM;PLA2G4A;PLOD3;RSRP1;SEMA3A;SLC2A1;STEAP3;SVEP1;TAF9;TBC1D2;TEAD3;TNFRSF11B;TOMM34;TOX;TPBG;TPI1;TRIB3;TSPAN31;TXN;VEGFA;VEGFB;WNT11;WSB2;ZNF302;ZNF711
compound::homoharringtonine homoharringtonine GSE70138                  86                  102               52                       34                        50          0.382353                      0.604651                                                       ARID5B;ARMC9;ATP6V1C1;CIRBP;CKB;CNTLN;CSRP2;EGLN3;ENO2;ERO1A;FAM162A;FILIP1L;GAPDH;GBE1;GOLPH3;GPM6B;HBP1;HEXB;HSD17B14;ITGB5;LAMP2;LIMA1;LOXL2;MBIP;MGST3;NAP1L3;NDUFA4L2;NR3C1;ODC1;P4HA2;PDE10A;PGK1;PKM;PLA2G15;PLOD3;PPT1;REXO2;SCARB1;SEMA3A;SPAG4;SRSF6;SV2A;TPBG;TPI1;TRIB3;TXN;UCHL1;VEGFA;VEGFB;WNT11;WSB2;ZNF711
         compound::l-655240          L-655240 GSE92742                  83                   82               51                       32                        31          0.447368                      0.621951                                               ARID5B;ATXN7L3B;AVPI1;BCL2L1;CNTLN;CRMP1;CSRP2;DIDO1;EDNRA;EDNRB;EFEMP1;ERO1A;FAM162A;FBLN5;FILIP1L;GDF15;GPM6B;HBP1;HEXB;IFNGR2;IGFBP5;ITGB5;LAMB3;LGALS3;LONP1;LUM;MBIP;MGST3;MTCH2;NAP1L3;NEU1;NREP;PDE10A;PFKP;PGK1;PKM;PLA2G15;PLA2G4A;PLOD3;REXO2;RSRP1;SCARB1;SEMA3A;TNFRSF11B;TNFRSF21;TOMM34;TPI1;TSPAN31;TXN;VEGFB;ZNF711
       compound::wye-125132        WYE-125132 GSE92742                  98                   91               49                       49                        42          0.350000                      0.538462                                                                   ARID5B;CIRBP;CNTLN;CRMP1;DIDO1;EDNRA;EFEMP1;EGLN3;ERO1A;FAM162A;FILIP1L;GAPDH;GPM6B;HBP1;HSPA9;IFNGR2;LAMP3;LONP1;LPAR1;LUM;MAGED2;MDK;MTCH2;NDUFA4L2;NR3C1;P4HA2;PFKP;PGK1;PKM;PLA2G4A;PRRX2;RSRP1;SCARB1;SEMA3A;SPAG4;STEAP3;SVEP1;TAF9;TEAD3;TNFRSF21;TOMM34;TPI1;TSPAN31;TXN;UCHL1;WNT11;WSB2;ZNF302;ZNF711
     compound::lestaurtinib      lestaurtinib GSE70138                  89                   78               48                       41                        30          0.403361                      0.615385                                                                        ARID5B;ATXN7L3B;CA9;CNTLN;CRMP1;DIDO1;EFEMP1;EGLN3;ENO2;ERO1A;FAM162A;FILIP1L;G0S2;GBE1;GDF15;HBP1;ITGB5;LAMB3;LAMP2;LGALS3;LOXL2;MAGED2;MTCH2;NAP1L3;NR3C1;P4HA2;PDE10A;PFKP;PGK1;PKM;PLA2G4A;PLOD3;POPDC3;REXO2;RSRP1;SCARB1;SLC2A1;SLC38A6;SPAG4;SVEP1;TNFRSF11B;TOMM34;TOX;TPI1;TXN;WSB2;ZNF302;ZNF711
          compound::osi-027           OSI-027 GSE92742                  89                   87               48                       41                        39          0.375000                      0.551724                                                                       ARID5B;CA9;CIRBP;CNTLN;CRMP1;DIDO1;EDNRA;EFEMP1;FAM162A;FBLN5;FILIP1L;HBP1;HSPA9;IFNGR2;LAMP3;LONP1;LPAR1;MAGED2;MDK;MT3;MTCH2;NDUFA4L2;NR3C1;P4HA2;PFKP;PGK1;PKM;PLA2G4A;PLOD3;POPDC3;PRRX2;REXO2;RSRP1;SCARB1;SEMA3A;SLC2A1;SVEP1;TAF9;TEAD3;TNFRSF11B;TNFRSF21;TOMM34;TPI1;TSPAN31;TXN;UCHL1;WSB2;ZNF302
      compound::gsk-2126458       GSK-2126458 GSE92742                  91                   89               48                       43                        41          0.363636                      0.539326                                                                           ARID5B;ATXN7L3B;CAT;CIRBP;CNTLN;CRMP1;DIDO1;EDNRA;EFEMP1;ERO1A;FAM162A;FILIP1L;GAPDH;GPC3;GPM6B;HBP1;HSPA9;IFNGR2;IGFBP5;LONP1;LUM;MAGED2;MDK;MTCH2;NR3C1;PFKP;PGK1;PKM;PLA2G4A;PLOD3;PRRX2;REXO2;RSRP1;SCARB1;SEMA3A;TAF9;TAF9B;TEAD3;TNFRSF11B;TNFRSF21;TOMM34;TOX;TPI1;TSPAN31;TXN;WNT11;WSB2;ZNF711
         compound::azd-8055          AZD-8055 GSE92742                  97                   92               48                       49                        44          0.340426                      0.521739                                                                         ARID5B;CAT;CDKN1A;CIRBP;CNTLN;CRMP1;DIDO1;EDNRA;EFEMP1;FAM162A;FILIP1L;GAPDH;HBP1;HSPA9;IFNGR2;LAMP3;LONP1;LPAR1;MAGED2;MDK;MTCH2;NAP1L3;NR3C1;P4HA2;PFKP;PGK1;PKM;PLA2G4A;PRRX2;RSRP1;SCARB1;SEMA3A;SPAG4;SVEP1;TAF9;TAF9B;TEAD3;TNFRSF11B;TNFRSF21;TOMM34;TOX;TPI1;TSPAN31;TXN;WNT11;WSB2;ZNF302;ZNF711

这些共享基因和各自特异基因的完整列表见 reversal-set comparison 输出；它们比共同基因上的 effect correlation 更直接回答两个 disease panel 是否把同一药物映射到相同的 reversal program。

### Recurrent reversal genes by panel

下面按 drug × release 统计每个 disease panel 中反复被判为 reversal 的基因；这是 panel-level reversal 频率，不是同一 LINCS effect 的重复测量。完整频率表同时保留 plastic、hydrogel、观察次数和药物数。

   panel    gene  n_drug_release_reversal  n_drug_release_observed  n_unique_drugs_reversal  n_unique_releases_reversal  fraction_reversal_among_observed
hydrogel     MIF                       79                      132                       51                           2                          0.598485
hydrogel   DDIT4                       81                      132                       49                           2                          0.613636
hydrogel    PGM1                       76                      132                       49                           2                          0.575758
hydrogel   ALDOA                       74                      132                       49                           2                          0.560606
 plastic  ARID5B                       72                      132                       49                           2                          0.545455
hydrogel  ARID5B                       72                      132                       49                           2                          0.545455
 plastic  GLT8D2                       67                      132                       49                           2                          0.507576
hydrogel SLC25A4                       67                      132                       49                           2                          0.507576
 plastic     PKM                       65                      132                       47                           2                          0.492424
hydrogel     PKM                       65                      132                       47                           2                          0.492424
 plastic     TXN                       62                      132                       47                           2                          0.469697
hydrogel     TXN                       62                      132                       47                           2                          0.469697
hydrogel   PGAM1                       70                      132                       46                           2                          0.530303
hydrogel    VAT1                       69                      132                       46                           2                          0.522727
 plastic   ERO1A                       66                      132                       46                           2                          0.500000
hydrogel   ERO1A                       66                      132                       46                           2                          0.500000
hydrogel  GABRA4                       58                      132                       46                           2                          0.439394
 plastic    TPI1                       65                      132                       45                           2                          0.492424
hydrogel    TPI1                       65                      132                       45                           2                          0.492424
 plastic    TBX1                       64                      132                       45                           2                          0.484848
 plastic   NFIL3                       59                      132                       45                           2                          0.446970
 plastic    GRK5                       54                      132                       45                           2                          0.409091
 plastic   LONP1                       65                      132                       44                           2                          0.492424
hydrogel   LONP1                       65                      132                       44                           2                          0.492424
 plastic  HSPA4L                       63                      132                       44                           2                          0.477273
 plastic     GSS                       61                      132                       44                           2                          0.462121
hydrogel  ZNF254                       61                      132                       44                           2                          0.462121
hydrogel   PTBP2                       60                      132                       44                           2                          0.454545
 plastic FILIP1L                       58                      132                       44                           2                          0.439394
hydrogel FILIP1L                       58                      132                       44                           2                          0.439394

#### Genes recurrent in both panels

   gene  n_unique_drugs_reversal_plastic  n_unique_drugs_reversal_hydrogel  fraction_reversal_among_observed_plastic  fraction_reversal_among_observed_hydrogel  min_fraction_reversal
 ARID5B                               49                                49                                  0.545455                                   0.545455               0.545455
  ERO1A                               46                                46                                  0.500000                                   0.500000               0.500000
    PKM                               47                                47                                  0.492424                                   0.492424               0.492424
   TPI1                               45                                45                                  0.492424                                   0.492424               0.492424
  LONP1                               44                                44                                  0.492424                                   0.492424               0.492424
   HBP1                               43                                43                                  0.492424                                   0.492424               0.492424
  P4HA2                               43                                43                                  0.484848                                   0.484848               0.484848
    TXN                               47                                47                                  0.469697                                   0.469697               0.469697
FAM162A                               42                                42                                  0.469697                                   0.469697               0.469697
   PGK1                               40                                40                                  0.469697                                   0.469697               0.469697
 SCARB1                               41                                41                                  0.454545                                   0.454545               0.454545
 TOMM34                               40                                40                                  0.454545                                   0.454545               0.454545
  NR3C1                               40                                40                                  0.446970                                   0.446970               0.446970
  MTCH2                               37                                37                                  0.446970                                   0.446970               0.446970
FILIP1L                               44                                44                                  0.439394                                   0.439394               0.439394
 EFEMP1                               40                                40                                  0.431818                                   0.431818               0.431818
 IFNGR2                               39                                39                                  0.424242                                   0.424242               0.424242
PLA2G4A                               44                                44                                  0.416667                                   0.416667               0.416667
  MGST3                               42                                42                                  0.416667                                   0.416667               0.416667
  HSPA9                               37                                37                                  0.416667                                   0.416667               0.416667
   WSB2                               36                                36                                  0.416667                                   0.416667               0.416667
  WNT11                               42                                42                                  0.409091                                   0.409091               0.409091
  GDF15                               35                                35                                  0.409091                                   0.409091               0.409091
  EDNRA                               43                                43                                  0.401515                                   0.401515               0.401515
PLA2G15                               39                                39                                  0.401515                                   0.401515               0.401515
  TEAD3                               39                                39                                  0.401515                                   0.401515               0.401515
   PFKP                               41                                41                                  0.393939                                   0.393939               0.393939
 CDKN1A                               37                                37                                  0.393939                                   0.393939               0.393939
  CNTLN                               37                                37                                  0.393939                                   0.393939               0.393939
  PLOD3                               41                                41                                  0.386364                                   0.386364               0.386364

## Secondary diagnostics

status agreement、共同基因上的 effect Spearman 和 effect sign concordance 仅用于确认数据拼接/可比性；由于 drug effect 来自同一套 LINCS perturbation，不能把它们解释为独立生物学复现。

解释边界：reversal 集合重叠仍不等于 hydrogel-specific causal mechanism。候选是否真正依赖 3D 环境，仍需结合 G×E / G×R×E、人体 niche 和选择性扰动验证。