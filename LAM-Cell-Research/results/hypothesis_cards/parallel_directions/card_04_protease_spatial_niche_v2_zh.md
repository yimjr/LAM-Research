# 研究线索卡 04（修订版）：多细胞 protease–antiprotease 空间生态位

## 观察

Visium、Visium HD 和 Xenium 中，LAMCORE-like spatial signal 与 protease signal 具有方向一致的空间关联。现在进一步使用实际单细胞表达估计 source state，再独立测量 protease/antiprotease。

## 重要修正

不预先规定 CTSK、MMP、ELANE 或 CTSS 属于 LAMCORE、LAF 或 immune。首轮单细胞 source attribution 显示不同基因的贡献模式不同，并存在较大的 unclassified 部分，因此当前不能声称某一细胞群是确定来源。

## Balance 定义

使用：

```text
proteolytic_balance_z = standardized protease activity - standardized antiprotease activity
```

不使用 protease/antiprotease 比值。三种空间技术分别分析，不合并 raw units、score 或 p 值。

## 当前等级

高价值探索性假说。尚未完成可复核 cyst wall/lesion-edge mask，也尚未证明多细胞共同贡献比单一来源更能预测病灶位置。

## QC 稳健性

全局候选池的固定敏感性检查为 baseline 140 个候选、strict-QC 85 个候选；LAM1/2/3/4 分别为 31/4/84/21 与 28/1/52/4。该空间 source attribution 尚未在两套候选池上分别重跑，因此本卡片暂不升级为高可信机制。

## 下一步验证

1. 提高 source reference 的细胞状态注释质量；
2. 对每个 protease gene 输出 source contribution 和 donor-specific variation；
3. 分别比较 protease ↑、antiprotease ↓ 和二者空间分离；
4. 在有可靠 lesion mask 后再进行 lesion-edge enrichment。
