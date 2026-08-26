# 已知程序阳性对照 Benchmark（首轮，已被 v2 取代）

> 本文件保留用于历史追溯。请以 [`benchmark_report_v2_zh.md`](benchmark_report_v2_zh.md) 和 `known_program_extended_benchmark_summary.csv` 为当前结论。

## 目的

在解释“不同 PatientID 之间未知程序匹配很弱”之前，先测试当前 top-gene/Jaccard matching 能否找回本来就已知存在的程序。若已知程序也经常匹配不上，未知程序的弱匹配不能被解释为患者特异性。

## 首轮结果

在 top-50、当前近似 gene set 和输出程序基因并集作为 null universe 的条件下，阳性恢复并不均匀：

- ECM_remodeling：pooled program detectable fraction 0.667；
- SLS_stem_like：0.444；
- CORE1：0.333；
- IS_inflammatory：0.333；
- LAM_myogenic_contractile / contractile：0.111；
- MDK_dormancy_persistence：0.111；
- CORE2：0.000。

这不是“CORE2 不存在”的证据，而是说明当前 top-gene matching 对程序定义、候选池、模块数和患者/assay 变化敏感。尤其 CORE2 不能稳定找回，因此现阶段不能把跨 PatientID 未匹配的未知程序称为患者特异性新生物学。

## 历史结论

方向一在这份首轮报告中暂时停留在方法校准阶段。后续 v2 已加入 donor-level expression score 和 NMF loading similarity；v2 显示 CORE2 可通过表达和 loading 层恢复，因此不能继续沿用本文件的“CORE2=0”作为当前结论。regulon/TF 与 pathway 层仍待补充；同时保留 pooled 与 donor-wise 两条发现路径。

## 可复核文件

- `known_program_matching_benchmark.csv`
- `known_program_matching_summary.csv`
- `benchmark_manifest.json`

注：当前 null universe 是“已生成程序基因 + 已知程序基因”的并集，不是完整 AnnData 基因宇宙，因此该结果是保守的校准结果，不应被解释为正式统计检验。
