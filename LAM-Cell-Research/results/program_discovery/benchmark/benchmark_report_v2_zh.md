# 已知程序阳性对照 Benchmark（扩展版）

top-gene Jaccard 与 expression benchmark 给出了不同信息。原先 CORE2 top-50 Jaccard recovery 为 0，但扩展分析显示 CORE2 在 LAM1–LAM4 的 donor-level expression enrichment 均为正，并且存在 NMF loading association。因此，top-gene list 不应单独决定 CORE2 是否被找回。

当前可稳定支持的 benchmark 类别包括：CORE1、CORE2、CORE3 identity、LAM myogenic/contractile、ECM remodeling。SLS、IS、MDK 等 treatment-associated program 的 expression 或 loading recovery 不均匀，不能直接用于未知 treatment-state 的新颖性判断。

TF/regulon/pathway 层仍未在当前 runtime 中完成，因此未知程序的最终新颖性判定仍受限。当前方向一可以继续处理 CORE/ECM 类候选，但 treatment-associated 候选仍保持 gate。

文件：

- `known_program_donor_expression_scores.csv`
- `known_program_nmf_loading_similarity.csv`
- `known_program_extended_benchmark_summary.csv`
- `extended_benchmark_manifest.json`
