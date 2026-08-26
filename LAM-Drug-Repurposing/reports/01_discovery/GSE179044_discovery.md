# GSE179044：TSC2 loss、rapamycin residual 与环境交互

## 研究问题

TSC2 loss 产生了什么异常，rapamycin 修复了什么，还留下了什么，hydrogel 是否改变 residual 或 genotype-dependent response？

## 数据与方法

- 数据：GSE179044，WT/TSC2−/− × vehicle/rapamycin × plastic/hydrogel 完整 2×2×2 设计；
- 脚本：`scripts/analyze_factorial.py`；
- 主要输出：`results/tables/GSE179044_factorial_contrasts.csv`；
- ratio 使用 shrinkage/moderated effect 后的 signed residual ratio，并设置 `|d0|` 门槛；
- hydrogel-specific residual 使用 rapamycin 条件下的 G×E simple interaction；
- environment-dependent escape 使用 G×R×E。

## 当前结果

hydrogel 中 ratio-eligible 基因包括 near-complete rescue、partial rescue + residual、persistent residual、worsened residual 和 direction reversal 五类。curated mTORC1 module 没有显示广泛的 hydrogel residual，保留信号更多涉及 ECM/invasion、myogenic、代谢和 autophagy。

NNMT、COL8A1、MIR210HG、SLC40A1、FBLN5、DCN 和 LUM 是探索性 hydrogel-conditioned residual 例子，但方向并不统一。G×R×E 信号稀少，当前不支持广泛的 environment-dependent escape program。

## 限制与下一步

GSE179044 每个条件只有两个生物学重复；内部 contrasts 是正交比较，不是独立复现。后续需要结合 GSE27982、机制数据和人体 LAM 状态判断 residual 是否具有跨模型和人体相关性。
