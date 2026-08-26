# LINCS 逐基因与稳定模块分析

## 分析范围

- 主面板：GSE179044 `tsc2_loss_plastic` top 150 up + top 150 down，共 300 genes。
- 药物：29 个 `replicated_concordant` compound，拆分为 20 个 reversal-only 和 9 个 mimic-only。
- LINCS release：GSE92742、GSE70138 分别分析，再做 cross-phase/cross-release comparison。
- 方向主指标：drug effect、reversal/mimic/neutral 状态、方向一致比例和 cell-line consistency；weighted contribution 仅作为辅助指标。
- 靶点遗传扰动按 `pert_type` 分开，未合并 shRNA、sh.cgs、过表达或 ligand。

## 当前发现

1. 300-gene 面板在两个 LINCS release 中共同可分析的基因为 204 个，而不是将缺失基因当作 neutral。其余基因保留为 `not_available`，没有参加跨 release 的复现判断。

2. 跨 release 的 gene clustering 得到 3 对达到第一版模块匹配标准的模块：
   - reversal-only：86 个共同基因，Jaccard 0.717；
   - mimic-only：70 个共同基因，Jaccard 0.614；
   - 另一对 reversal-only：27 个共同基因，Jaccard 0.519。
   这些是稳定模块候选，不等同于已确定的生物学通路。

3. Lestaurtinib 在两个 release 都有较多稳定 reversal genes（GSE92742 92 个，GSE70138 89 个），且与 PI3K/mTOR 药物存在部分重叠，但并非完全等同于 PI3K/mTOR reversal pattern。其共同跨 release 方向为 93 genes，14 genes 方向相反，其余主要为 neutral/weak。

4. Lestaurtinib 靶点轴的本地证据具有明显 modality 和 release 限制：
   - GSE70138 中没有可用的对应遗传扰动，因此没有 cross-release target validation；
   - GSE92742 中 RET 的 shRNA 和 `sh.cgs` 方向一致性较好，暂列为 supportive；
   - NTRK3 接近但尚属 weak；
   - JAK2、FLT3 和其他 NTRK 结果多数为 weak，FLT3 shRNA 有 discordant 信号。
   这支持“RET/NTRK 轴值得优先验证”的研究假说，但不能证明 Lestaurtinib 由单一靶点介导。

5. QL-X-138 的 drug pattern 与 BTK/MNK 遗传扰动有部分重叠：
   - BTK `sh.cgs` 暂列 supportive；
   - BTK 其他 modality 主要为 weak；
   - MKNK1/MKNK2 目前主要为 weak，尚不能支持 MNK 为主要介导轴。

6. 9 个 mimic-only 药物在两个 release 中形成了一个较稳定的共同模块（70 个共同基因，Jaccard 0.614）。其中包含 DDIT3、SQSTM1、HSPA9、LAMP3、SLC2A1、TRIB3 等 stress/proteostasis 相关基因，提示 mimic 信号可能包含共同的应激与蛋白稳态反应；但这还不能区分普通细胞毒性和 TSC2-loss 特异脆弱性。

7. 既有八个功能模块的事后 overlap 分析尚未得到可靠 FDR 显著结果；因此当前不把 ECM、mTOR、autophagy 或 stress 之一称为已确认主导程序。

## 下一步

- 对 3 对跨 release 稳定模块进行 GO/Reactome/MSigDB 富集；background 使用实际进入聚类的 204 个共同可分析基因。
- 对 Lestaurtinib 的 RET、NTRK3 和 QL-X-138 的 BTK 进行选择性药物或 CRISPR/shRNA 外部验证。
- 加入 generic cytotoxicity、proteasome/microtubule/NEDD8 stress controls，判断 mimic 模块是否只是泛应激。
- 将稳定 drug × gene modules 映射到 GSE135851/GSE302356 的 LAMCORE/LAF 状态；目前不生成 Tier 1 联合治疗候选。
