# LAM Research Hypothesis Card：high_confidence / program_3

## 当前分级

探索性假说（仅来自 GSE135851 同一队列；尚未达到独立验证标准）。

## 观察到了什么

在 `high_confidence` 候选池的 pooled NMF 中出现了 `program_3`。当前排名靠前的基因为：MALAT1, IL32, PROM1, CCL14, CTSK, ARHGDIB, STK17A, DES, HSP90AB1, CDC37, CCL4L2, HSP90AA1。

## donor 重复性

逐 donor 独立程序发现后，达到当前基因重叠标准的 donor：。这不是“在 pooled 模型中打分”，而是独立 donor 发现的初步匹配；若为空，说明目前没有足够证据称其跨 donor 重复。

## 与已知框架的关系

最强的已知程序匹配：protease_ECM_niche (microenvironment; overlap=4)。部分重叠不会自动淘汰该候选；还需要判断它是否是已知状态在 LAMCORE 中的新实现、是否有 LAM-specific 基因或新的 TF/regulon。

## 外部证据

当前可用的外部 AnnData：GSE190260, GSE217108, GSE302356。ATAC、空间或蛋白证据尚未用于本卡片的结论。

## 替代解释

- donor、assay 或批次特异信号；
- cell cycle、doublet、测序深度或低质量造成的程序；
- 现有 CORE/SLS/IS/ECM 等状态的部分投影；
- 宽松候选集中仍未确认的细胞身份。

## 下一步验证

在 GSE190260、GSE217108 和 GSE302356 中按 PatientID 重复发现；按 donor 独立提取程序；比较已知程序解释比例；再检查 ATAC、空间或蛋白支持。未经这些步骤，不将其命名为新 LAMCORE 亚型或机制。

## 新颖性 / 可信度 / 优先级

- 新颖性：未评估；
- 当前可信度：低到中等，仅为同队列候选；
- 推荐优先级：中，取决于外部 donor 是否独立重现。
