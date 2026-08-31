# 20 个 consensus state 的生物学解释

以下 analogue 是基于现有 DE、program、LAMCORE、LOO、normal/boundary 和 latent 邻域的综合解释，不是 upstream 已验证 `cell_type`。由于 `state_by_cell_type.csv` 的 5,378 个细胞均为 `unknown`，所有命名都保留 evidence boundary。

## State 1 — ciliated/airway-like epithelial

- 解释类别：relatively_clear_normal-lineage analogue。
- 支持证据：Top DE markers include DNAI1, DNAH9, DCDC1 and other cilia genes; strong ciliated geometry/identity signal in Stage16 diagnostic summary.
- 冲突证据：Upstream cell_type is unknown; Stage13 programs also contain hormone/HOX and interstitial signals.
- 与 LAM 的关系：Not supported as a LAM-core state; likely a normal airway-lineage contaminant within the high-recall pool.
- 不确定性：Moderate: cell-type label is inferred from expression, not an upstream verified annotation.
- 细胞/覆盖/结构证据：703 cells；10 patients；4 datasets；structural=0.9628848873718564；biological=0.3893193992572874。
- 表达/程序摘要：top DE markers=ZBBX(8.84), TTC29(8.51), DNAI1(8.23), ADGB(7.95), DNAH9(7.94), VWA3B(7.92), DCDC1(7.78), ARMC4(7.69)；top program deltas=hormone_related=0.301; HOX_PBX=0.268; SLS_stem_like=0.182; IL6_AT2_repair=0.044; TGFbeta_fibroblast=-0.017。

## State 2 — undetermined rare substate

- 解释类别：insufficient evidence。
- 支持证据：One cell only; mTOR/MDK program signal.
- 冲突证据：No patient or dataset-level reproducibility; no formal DE support.
- 与 LAM 的关系：Cannot establish LAM relationship.
- 不确定性：Very high; do not assign a biological analogue.
- 细胞/覆盖/结构证据：1 cells；1 patients；1 datasets；structural=0.5006857804766124；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=mTOR_translation=1.114; MDK_dormancy_persistence=0.619; SLS_stem_like=0.256; hypoxia_stress=0.225; IL6_AT2_repair=0.116。

## State 3 — myeloid/inflammatory-like, provisional

- 解释类别：insufficient evidence。
- 支持证据：Four cells; macrophage_TREM2_TYROBP and inflammatory programs.
- 冲突证据：Only two patients and no formal DE/LOO biological support.
- 与 LAM 的关系：Not sufficient for LAM interpretation.
- 不确定性：Very high.
- 细胞/覆盖/结构证据：4 cells；2 patients；2 datasets；structural=0.3414992897788756；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=macrophage_TREM2_TYROBP=0.373; IS_inflammatory=0.281; protease_ECM_niche=0.140; hypoxia_stress=0.127; IL6_AT2_repair=0.055。

## State 4 — mixed immune–mesenchymal/interstitial

- 解释类别：mixed or uncertain。
- 支持证据：LAM-myogenic/uterine-smooth and normal-interstitial programs; top markers include IGLV3-21, TNFRSF17, F13A1 and CD163L1.
- 冲突证据：Immune/plasma-cell-like markers conflict with a simple LAM interpretation; upstream cell_type unknown.
- 与 LAM 的关系：May contain LAM-shared myogenic signal but is not LAM-specific.
- 不确定性：High; mixed profile and no single decisive lineage assignment.
- 细胞/覆盖/结构证据：550 cells；10 patients；4 datasets；structural=0.7964477901985803；biological=0.4653757956746926。
- 表达/程序摘要：top DE markers=IGLV3-21(5.71), TNFRSF17(4.89), CLEC1B(4.57), SDS(4.24), F13A1(4.20), ITGAD(4.10), CD163L1(4.03), MIR137HG(3.98)；top program deltas=CORE2=0.113; CORE3_identity=0.039; protease_ECM_niche=-0.014; cell_cycle=-0.045; TGFbeta_fibroblast=-0.049。

## State 5 — macrophage/myeloid-like

- 解释类别：relatively clear normal-lineage analogue。
- 支持证据：FABP4, APOC1, MARCO, RETN and CCL23; macrophage_TREM2_TYROBP program delta 2.238, the strongest program signal in the atlas.
- 冲突证据：LAM-myogenic, mTOR and protease programs are also elevated, reflecting shared tissue programs or admixture.
- 与 LAM 的关系：Not LAM-core; LAM-associated programs should not override the strong myeloid analogue.
- 不确定性：Low-to-moderate for the broad analogue; exact macrophage subtype is not resolved.
- 细胞/覆盖/结构证据：864 cells；10 patients；4 datasets；structural=0.950408815502366；biological=0.559758193072894。
- 表达/程序摘要：top DE markers=FABP4(5.90), APOC1(5.40), AGRP(5.39), RETN(5.27), MARCO(5.16), CD52(5.11), CCL23(5.10), IGSF6(4.94)；top program deltas=macrophage_TREM2_TYROBP=2.238; mTOR_translation=0.531; protease_ECM_niche=0.493; MDK_dormancy_persistence=0.279; cell_cycle=0.150。

## State 6 — rare myeloid/AT2-mixed substate

- 解释类别：insufficient evidence。
- 支持证据：Three cells; macrophage and IL6/AT2-repair programs.
- 冲突证据：Only two patients, no formal biological reproducibility.
- 与 LAM 的关系：Cannot distinguish LAM biology from a mixed rare state.
- 不确定性：Very high.
- 细胞/覆盖/结构证据：3 cells；2 patients；2 datasets；structural=0.4086910368110567；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=macrophage_TREM2_TYROBP=0.643; mTOR_translation=0.221; MDK_dormancy_persistence=0.203; protease_ECM_niche=0.122; SLS_stem_like=0.060。

## State 7 — AT2-like alveolar epithelial/repair

- 解释类别：relatively clear normal-lineage analogue。
- 支持证据：SFTPC, SFTPA1, SFTPA2, KRT78 and SFTA2; IL6_AT2_repair is the leading program; Stage22 labels its adjacency ordinary lineage adjacency.
- 冲突证据：Some myogenic/macrophage program overlap is present in the candidate pool.
- 与 LAM 的关系：Not a LAM state on current evidence.
- 不确定性：Low-to-moderate; exact repair-state biology remains unresolved.
- 细胞/覆盖/结构证据：576 cells；11 patients；4 datasets；structural=0.6826819842967636；biological=0.3978258755809657。
- 表达/程序摘要：top DE markers=PLA2G1B(4.94), AGTR2(4.92), SFTPC(4.75), SFTPA2(4.75), KRT78(4.70), SFRP5(4.70), SFTA2(4.58), SFTPA1(4.55)；top program deltas=IL6_AT2_repair=0.602; hormone_related=0.064; TGFbeta_fibroblast=-0.020; HOX_PBX=-0.032; SLS_stem_like=-0.033。

## State 8 — undetermined TGFβ/interstitial rare substate

- 解释类别：insufficient evidence。
- 支持证据：One cell; TGFbeta_fibroblast, hormone and HOX/PBX signals.
- 冲突证据：No replication or DE support.
- 与 LAM 的关系：Cannot establish LAM relationship.
- 不确定性：Very high.
- 细胞/覆盖/结构证据：1 cells；1 patients；1 datasets；structural=0.501159222532781；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=TGFbeta_fibroblast=0.501; hormone_related=0.348; HOX_PBX=0.148; MDK_dormancy_persistence=0.115; cell_cycle=-0.071。

## State 9 — mixed repair/interstitial-like

- 解释类别：mixed or uncertain。
- 支持证据：IL6_AT2_repair, hormone and TGFβ programs; 406 cells across four datasets.
- 冲突证据：Only five patients, no formal patient-level DE support and weak independent biological evidence.
- 与 LAM 的关系：May reflect tissue-repair programs rather than LAM identity.
- 不确定性：High.
- 细胞/覆盖/结构证据：406 cells；5 patients；4 datasets；structural=0.7286135596598886；biological=0.0。
- 表达/程序摘要：top DE markers=；top program deltas=IL6_AT2_repair=0.374; hormone_related=0.230; TGFbeta_fibroblast=0.164; IS_inflammatory=0.025; hypoxia_stress=-0.025。

## State 10 — rare myeloid/repair mixed state

- 解释类别：insufficient evidence。
- 支持证据：Four cells; macrophage and IL6/AT2-repair programs.
- 冲突证据：Single dataset and no patient support.
- 与 LAM 的关系：No reliable LAM conclusion.
- 不确定性：Very high.
- 细胞/覆盖/结构证据：4 cells；1 patients；1 datasets；structural=0.4467594784634767；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=macrophage_TREM2_TYROBP=0.647; IL6_AT2_repair=0.245; TGFbeta_fibroblast=0.126; protease_ECM_niche=0.101; CORE2=0.080。

## State 11 — rare fibroblast/HOX-hormone-like state

- 解释类别：insufficient evidence。
- 支持证据：Three cells; HOX/PBX, hormone and TGFβ programs.
- 冲突证据：Single dataset and no biological replication.
- 与 LAM 的关系：No reliable LAM conclusion.
- 不确定性：Very high.
- 细胞/覆盖/结构证据：3 cells；1 patients；1 datasets；structural=0.4719788426043272；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=HOX_PBX=0.495; hormone_related=0.359; TGFbeta_fibroblast=0.311; LAF_niche=0.186; normal_lung_interstitial=0.132。

## State 12 — endothelial/lymphatic endothelial-like

- 解释类别：relatively clear normal-lineage analogue。
- 支持证据：MMRN1, TM4SF18, SELE, STAB2, CCL21 and FLT4; Stage22 dominant adjacent direction is endothelial.
- 冲突证据：Some LAM-myogenic/normal-interstitial programs overlap and high-recall candidate selection contributes to inclusion.
- 与 LAM 的关系：Not LAM-core; current branch null supports ordinary lineage adjacency.
- 不确定性：Moderate; endothelial and lymphatic components are not separated into new states.
- 细胞/覆盖/结构证据：605 cells；9 patients；4 datasets；structural=0.681160600023974；biological=0.3264022951522951。
- 表达/程序摘要：top DE markers=MMRN1(4.53), TM4SF18(4.35), SELE(4.35), UNC5A(4.26), STAB2(4.15), CCL21(4.05), FLT4(3.95), SCN3B(3.76)；top program deltas=IL6_AT2_repair=0.170; hormone_related=0.154; TGFbeta_fibroblast=0.130; HOX_PBX=0.055; biomarker_VEGFD_PMEL_CCL14_MMP8=-0.024。

## State 13 — rare HOX/CORE3-mixed substate

- 解释类别：insufficient evidence。
- 支持证据：Six cells; HOX/PBX and CORE3_identity signals.
- 冲突证据：Only three datasets, no patient-level DE support and mixed programs.
- 与 LAM 的关系：A possible LAM-like signal cannot be separated from sampling noise.
- 不确定性：Very high.
- 细胞/覆盖/结构证据：6 cells；3 patients；3 datasets；structural=0.3388228531216182；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=HOX_PBX=0.387; hormone_related=0.266; CORE3_identity=0.251; cell_cycle=0.100; TGFbeta_fibroblast=0.079。

## State 14 — rare LAM-myogenic/contractile-like substate

- 解释类别：insufficient evidence。
- 支持证据：Three cells with high LAM-myogenic, uterine-smooth, CORE1 and CORE3 program deltas.
- 冲突证据：One dataset, no patient support and no formal DE replication.
- 与 LAM 的关系：Interesting LAM-like signal, but not a reproducible state claim.
- 不确定性：Very high.
- 细胞/覆盖/结构证据：3 cells；1 patients；1 datasets；structural=0.4193117513949105；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=LAM_myogenic_contractile=1.894; lineage_uterine_smooth_muscle=1.894; CORE1=1.374; MDK_dormancy_persistence=0.875; mTOR_translation=0.765。

## State 15 — LAM-rich contractile/mesenchymal candidate; provisional LAM-core anchor

- 解释类别：LAM-associated candidate。
- 支持证据：200 cells; 7 patients/4 datasets; LAMCORE median 0.5125 in Stage18; CORE1/CORE3/LAM-myogenic/ECM programs; 49 author-style cells where author labels are available; removal of LAM1163 leaves 73 cells and LAMCORE median 0.4683.
- 冲突证据：LAM1163 contributes 127/200 with 6.8301-fold composition enrichment; author labels are available only in GSE135851; Stage11 biological reproducibility is 0.374758; Stage18 did not formally upgrade the anchor.
- 与 LAM 的关系：Strongest LAM-rich state in this project, but still a provisional reference candidate rather than a diagnostic classifier or formal cross-patient anchor.
- 不确定性：Moderate-to-high due to patient composition and incomplete independent author annotation.
- 细胞/覆盖/结构证据：200 cells；7 patients；4 datasets；structural=0.8541110574273367；biological=0.37475757007007。
- 表达/程序摘要：top DE markers=HTN3(7.12), MMP11(6.47), PAGE4(6.25), PGM5-AS1(6.25), HOXA11(6.10), EMX2(5.64), PLAT(5.63), SFRP1(5.54)；top program deltas=LAM_myogenic_contractile=2.051; lineage_uterine_smooth_muscle=2.051; CORE1=1.456; CORE3_identity=1.145; ECM_remodeling=1.035。

## State 16 — immune/T-NK-adjacent mixed state; no confirmed transition

- 解释类别：mixed or uncertain。
- 支持证据：396 cells; CORE1/CORE3/LAM-myogenic programs; per-patient slopes were negative in 7/7 represented patients and LOPO slopes in 10/10 omissions, but the direction-matched left-tail test did not separate this branch from its local matched null.
- 冲突证据：Only 8 patients contribute State16 cells in the direct 1-hop branch criterion; immune signal is heterogeneous, upstream cell_type is unknown, and Stage11 biological reproducibility is 0.323088.
- 与 LAM 的关系：Geometrically adjacent to State15 with a directional LAMCORE decrease, but current corrected evidence does not support a LAM-preserving transition label; ordinary adjacency and mixed biology remain plausible.
- 不确定性：High; no temporal or lineage-transition evidence.
- 细胞/覆盖/结构证据：396 cells；10 patients；4 datasets；structural=0.7196327799664204；biological=0.3230880484706342。
- 表达/程序摘要：top DE markers=PCP4(5.73), HOXC10(5.42), TMEM196(5.32), CD3G(5.16), LINC00906(4.99), OR2L5(4.98), LINC00402(4.95), FABP7(4.94)；top program deltas=LAM_myogenic_contractile=1.286; lineage_uterine_smooth_muscle=1.286; CORE1=0.948; CORE3_identity=0.696; CORE2=0.220。

## State 17 — mesothelial/secretory epithelial-like, uncertain

- 解释类别：mixed or uncertain。
- 支持证据：121 cells; ITLN1, CALB2, CPB1, CPA4 and ANXA8 among top markers; mTOR/ECM/inflammatory programs.
- 冲突证据：Only 3 supported patients, 3 datasets, and the marker profile is not a single clean mesothelial signature.
- 与 LAM 的关系：No evidence for LAM-core; shared ECM/inflammatory signals are non-specific.
- 不确定性：High.
- 细胞/覆盖/结构证据：121 cells；5 patients；3 datasets；structural=0.9131796417279768；biological=0.4048718447231132。
- 表达/程序摘要：top DE markers=TMEM151A(10.01), CALB2(9.95), ITLN1(9.90), CPB1(9.45), IL20(8.81), BNC1(8.39), CPA4(8.19), ANXA8(8.02)；top program deltas=mTOR_translation=1.163; ECM_remodeling=0.857; IS_inflammatory=0.509; normal_lung_interstitial=0.499; hypoxia_stress=0.376。

## State 18 — pericyte/VSMC/smooth-muscle-like

- 解释类别：relatively clear normal-lineage analogue。
- 支持证据：COX4I2, FOXC2, CASQ2, KCNA5, FHL5 and HIGD1B; high uterine-smooth/LAM-myogenic and CORE1 signal; strong structural stability 0.940033.
- 冲突证据：LAM and smooth-muscle programs overlap biologically, so ACTA2/myogenic signal alone cannot identify LAM; State18 did not meet the direct-connection branch rule.
- 与 LAM 的关系：Important LAM mimic/comparator; current evidence favors ordinary VSMC/pericyte adjacency rather than a LAM branch.
- 不确定性：Moderate.
- 细胞/覆盖/结构证据：174 cells；9 patients；4 datasets；structural=0.9400330318823972；biological=0.4089603116630143。
- 表达/程序摘要：top DE markers=AC093390.1(8.55), CASQ2(8.41), KCNA5(8.39), COX4I2(7.38), FHL5(7.36), HIGD1B(7.36), FOXC2(7.30), ATP1A2(6.96)；top program deltas=lineage_uterine_smooth_muscle=2.437; LAM_myogenic_contractile=2.437; CORE1=1.934; mTOR_translation=0.590; ECM_remodeling=0.287。

## State 19 — undetermined rare interstitial/hormone-like substate

- 解释类别：insufficient evidence。
- 支持证据：Two cells; hormone, normal-interstitial and LAF programs.
- 冲突证据：Two patients, two datasets and no biological replication.
- 与 LAM 的关系：Cannot establish LAM relationship.
- 不确定性：Very high.
- 细胞/覆盖/结构证据：2 cells；2 patients；2 datasets；structural=0.4063892162606694；biological=0.0。
- 表达/程序摘要：top DE markers=not available；top program deltas=hormone_related=0.642; normal_lung_interstitial=0.596; LAF_niche=0.578; TGFbeta_fibroblast=0.452; ECM_remodeling=0.408。

## State 20 — fibroblast/lung interstitial-like

- 解释类别：relatively clear normal-lineage analogue。
- 支持证据：PI16, SFRP2, SCARA5, MYOC, DPT, C7, COMP and CXCL14; normal-lung-interstitial/LAF/ECM programs; Stage22 ordinary lineage adjacency.
- 冲突证据：Macrophage/MDK/LAM-myogenic shared programs and candidate enrichment can blur boundaries.
- 与 LAM 的关系：Not LAM-core; a principal normal interstitial/fibroblast comparator.
- 不确定性：Low-to-moderate for broad analogue; fibroblast subtype not resolved.
- 细胞/覆盖/结构证据：756 cells；10 patients；4 datasets；structural=0.9125450360451768；biological=0.4207180402657939。
- 表达/程序摘要：top DE markers=PI16(7.90), SFRP2(7.75), SCARA5(6.78), MYOC(6.77), DPT(6.58), C7(6.37), COMP(6.09), CXCL14(5.97)；top program deltas=normal_lung_interstitial=1.701; LAF_niche=1.607; ECM_remodeling=0.996; MDK_dormancy_persistence=0.396; mTOR_translation=0.375。
