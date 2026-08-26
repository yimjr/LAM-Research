# 研究线索卡 06（修订版）：Rapamycin 后 ECM/protease 程序保留

## 定义

不使用“rapamycin 后不显著下降”定义 retained。对每个基因计算：

```text
E_pre  = TSC2-loss_vehicle - WT_vehicle
E_post = TSC2-loss_rapamycin - WT_rapamycin
suppression_fraction = 1 - E_post/E_pre
```

## 当前观察

GSE179044 中 ECM/protease 相关程序存在 effect-size 定义的 retained 或 enhanced genes。ELANE 在 hydrogel 和 plastic 中均满足重复方向一致的 partial-retention 条件，是目前首个跨环境 protease 候选；MMP2 仅在 hydrogel 中满足条件，plastic 中不满足资格条件。

## 当前结论

可以保留为高价值候选机制：mTOR 抑制可能控制生长，但不完全清除 matrix-related pathology。ELANE 只是扰动模型中的候选 retained gene，不能称为患者级 sirolimus persistence，也不能称为已证实的 protease-resistant mechanism。

## QC 稳健性

GSE179044 是处理后的扰动表达矩阵，不具备与单细胞候选池同义的 baseline/strict-QC 细胞过滤层；该卡片的模型内稳健性由两个 biological replicate 的方向一致性和 suppression fraction 敏感性分析承担。人类 LAM 交叉验证仍须通过 baseline 140 与 strict 85 候选池检查。

## 下一步

将 retained genes 与人类 LAMCORE、空间 source attribution、proteolytic balance 和 lesion/niche 数据交叉验证。
