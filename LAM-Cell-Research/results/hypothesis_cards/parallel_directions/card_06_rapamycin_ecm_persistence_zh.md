# 研究线索卡 06：mTOR 抑制后 ECM/蛋白酶程序仍然保留

## 当前观察

在 GSE179044 的 TSC2-knockout 模型中，rapamycin 对 ECM remodeling 和 protease/ECM niche 的平均程序分数抑制很弱，且 hydrogel 与 plastic 两种环境方向相近。GSE84476 的 TSC2-null LAM 来源细胞中，rapamycin 相对 siCtrl 的 ECM remodeling 和 protease/ECM niche 分数反而升高。

## 解释边界

这是扰动模型之间方向相容的候选机制，不是患者级 sirolimus 耐药证据。GSE179044 有两个生物学重复，适合继续做 replicate-aware gene-level 分析；GSE84476 每个条件的样本量有限，当前只能作为描述性支持。程序分数来自不同矩阵和尺度，不能直接比较绝对数值。

## 候选机制

mTOR 抑制可能压低部分生长/炎症程序，但 ECM 环境或蛋白酶相关程序未被同步清除，甚至在某些 TSC2-null 模型中相对增强。这可能代表一种“抑制生长而保留基质适应”的候选持续机制。

## 证据等级

高价值探索性假说。尚不能称为患者治疗耐受，也不能称为新机制；需要基因级差异、重复一致性、患者状态对应和 ECM 条件实验验证。

## 可检验预测

1. rapamycin 后保留的核心 ECM/protease 基因应在两个 GSE179044 重复中方向一致；
2. ECM hydrogel 条件下该程序的保留应强于 plastic，或出现明确的环境交互；
3. 该程序应能映射到人类 LAMCORE/空间 protease niche，而不是只出现在体外模型。

## 文件

`results/perturbation/GSE179044_program_contrasts.csv`、`GSE84476_program_contrasts.csv` 及对应 analysis manifest。
