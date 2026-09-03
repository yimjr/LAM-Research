# LAM Research

一个基于公开生物医学数据的 **LAM（Lymphangioleiomyomatosis，淋巴管平滑肌瘤病）计算研究项目**。

本项目尝试重新连接公开的转录组、单细胞、空间组学、药物扰动等数据，探索与 **LAM 发病机制、肺组织破坏、药物治疗、免疫识别和细胞异质性** 有关的研究问题。

目前主要包含四个相对独立、同时可以相互提供证据的研究方向：

| 方向 | 核心问题 |
| --- | --- |
| [LAM Cell Research](LAM-Cell-Research/) | LAM 细胞具有哪些生物学程序？它们如何适应肺部环境、参与组织破坏？ |
| [LAM Drug Repurposing](LAM-Drug-Repurposing/) | 已有药物中，有没有可能逆转 LAM/TSC2 缺失造成的异常状态？西罗莫司没有完全改变的部分是什么？ |
| [LAM Immune Visibility](LAM-Immune-Visibility/) | LAM 细胞表达了哪些可能被免疫系统识别的特征？这些特征与抗原呈递之间有什么关系？ |
| [LAM State Modeling](LAM-State-Modeling/) | 能否在多个单细胞数据集中更稳定地识别 LAM 细胞，并进一步研究这些细胞的状态和生物学特征？ |

LAM State Modeling 在方法上引入了前三个方向中没有使用的 **神经网络单细胞模型 scVI**，用于整合不同患者、不同数据集中的高维基因表达信息。

这一方向的核心目标是提高 LAM-rich 细胞群识别的可靠性，并进一步研究这些细胞的共同特征和内部异质性。

[English](README.md)

---

## 为什么做这个项目？

LAM 是一种罕见病。一个重要分子基础是 **TSC1/TSC2 功能异常及其导致的 mTOR 通路持续活跃**。

目前 sirolimus（西罗莫司/雷帕霉素）已经显著改变了 LAM 的治疗，但仍有很多问题值得继续研究：

* LAM 细胞并不是完全相同的一群细胞，它们可能具有不同状态；
* mTOR 抑制能够控制很多异常，但部分细胞程序仍可能保留下来；
* LAM 为什么尤其容易在肺内形成破坏性病变，仍有很多机制可以继续探索；
* 肺囊性破坏可能涉及 LAM 细胞、成纤维细胞、免疫细胞和蛋白酶系统共同作用；
* 已有药物中可能存在能够作用于 mTOR 之外异常程序的候选；
* LAM 细胞具有一些比较特殊的黑色素细胞/间叶细胞相关特征，但这些特征与免疫识别之间的关系仍有探索空间。

与此同时，过去几年已经积累了不少公开的 LAM 单细胞、空间组学、扰动实验和药物数据库。

因此，这个项目的基本想法是：

> **把已经公开、但原本服务于不同研究问题的数据重新连接起来，从中寻找新的研究问题、机制线索和可进一步实验验证的假说。**

这个仓库更关注“还能从现有数据中发现什么值得继续研究的问题”，希望最终得到能够交给实验研究者继续验证的候选机制、候选药物和研究方向。

---

# 1. LAM Cell Research

📁 [LAM-Cell-Research](LAM-Cell-Research/)

## 为什么研究 LAM 细胞状态？

单细胞测序已经表明，LAM 病灶中存在被称为 **LAMCORE** 的核心 LAM 细胞群，并进一步观察到了不同的细胞状态。

但“LAM 细胞是什么”仍然可以继续向下追问：

* LAMCORE 内部是否还有没有被充分描述的状态？
* 不同患者的 LAM 细胞是否共享某些隐藏的表达程序？
* 这些状态是几个清晰的亚型，还是一个连续变化的状态空间？
* LAM 细胞为什么能够在肺中生存和扩张？
* LAM 病灶周围的其他细胞是否和 LAM 细胞一起参与肺组织破坏？
* sirolimus 控制 mTOR 后，哪些与 ECM、侵袭或蛋白酶有关的程序还存在？

因此，这一方向从公开单细胞数据出发，同时结合空间转录组、药物扰动和其他组织参考数据研究这些问题。

---

## 基本研究思路

```text
公开 LAM 单细胞数据
        ↓
重新识别 LAMCORE-like 细胞
        ↓
复现已有 LAMCORE 状态
        ↓
研究患者之间和患者内部的状态差异
        ↓
寻找新的表达程序
        ↓
结合空间组学、药物扰动和其他数据解释这些程序
```

分析中既会研究单个基因，也会更多关注 **gene program（基因程序）**：

> 一组经常一起变化、可能共同反映某种细胞功能或状态的基因。

---

## 当前主要探索方向

### 1. LAMCORE 状态是否更像连续变化？

目前重新实现得到约 **140 个 LAMCORE-like 候选细胞**，数量与原研究报告的约 125 个处于相近范围；采用更加严格的 QC 后保留约 **85 个**。

进一步使用不同聚类方法分析时，没有看到一个能够稳定覆盖多个患者的清晰新亚型。

例如 HDBSCAN 分析中：

* 140 个候选细胞中形成 2 个局部密度簇；
* 只有 28 个细胞进入这些簇；
* 主要来自 LAM3/LAM4；
* 没有形成覆盖所有 donor 的稳定离散亚型。

因此目前更值得考虑的是：

> **LAMCORE 内部的异质性可能更多体现为连续或局部变化，而不是几个边界非常清楚的新细胞亚型。**

---

### 2. 寻找尚未被充分描述的 LAM 表达程序

项目使用 pooled NMF、donor-wise NMF、meta-program matching 等方式寻找在多个细胞中共同出现的表达程序。

首轮候选程序涉及：

* ECM / extracellular matrix；
* protease；
* migration；
* proliferation；
* hormone-related signals；
* LAM lineage；
* microenvironment interaction。

目前没有出现一个能够在多个独立患者中稳定重新发现、并可直接定义为“全新 LAMCORE 程序”的结果。

这提示重要差异可能来自已知程序的重新组合、强弱变化和微环境依赖，而不一定表现为全新的细胞类型。

---

### 3. Protease–antiprotease 空间生态位

LAM 最重要的病理特征之一是肺组织逐渐形成大量囊性破坏。

因此项目进一步研究：

> **破坏肺组织的蛋白酶信号究竟来自 LAM 细胞本身，还是 LAMCORE、成纤维细胞、免疫细胞等多种细胞共同形成一个 proteolytic niche（蛋白水解生态位）？**

目前已经在 Visium、Visium HD 和 Xenium 等空间数据中观察到 LAMCORE-like spatial signal 与 protease signal 方向一致的空间关联。

关注基因包括：

* `CTSK`
* `MMP` family
* `ELANE`
* `PRTN3`
* `CTSS`

同时分析 antiprotease，从而形成：

```text
protease activity
        -
antiprotease activity
        ↓
proteolytic balance
```

长期问题是：LAM 肺囊性破坏是否来自一个由多种细胞共同维持的局部蛋白水解环境。

---

### 4. LAM 细胞进入肺后是否获得新的适应程序？

项目比较：

```text
Normal uterus
LAM uterus
Normal lung
Pulmonary LAM
```

并计算类似：

```text
Pulmonary LAM - Normal lung
        vs.
LAM uterus - Normal uterus
```

的差异。

目前结果提示：

* `lung_adaptation` 程序更符合进入肺后获得的状态；
* ECM 程序在 LAM uterus 和 pulmonary LAM 中都增强。

因此可以将 pulmonary LAM 暂时理解为：

```text
LAM transformation / lineage program
        +
lung-acquired adaptation program
        =
pulmonary LAM state
```

---

### 5. Rapamycin 后仍保留的 ECM / protease 程序

在 TSC2-loss 扰动数据中，部分 ECM/protease 相关基因在 rapamycin 后仍存在一定程度的异常。

当前例子包括：

* `ELANE` 在 plastic 和 hydrogel 两种环境下均表现出方向一致的 partial retention；
* `MMP2` 在 hydrogel 环境中表现出 retention 信号。

这支持继续研究：

> **mTOR inhibition 可能很好地控制细胞生长，同时仍有部分 matrix-related pathology 值得单独研究。**

这条线也与药物再利用项目直接产生联系。

---

# 2. LAM Drug Repurposing

📁 [LAM-Drug-Repurposing](LAM-Drug-Repurposing/)

## 为什么研究药物再利用？

LAM 已经存在有效治疗药物 sirolimus，但 sirolimus 主要针对 mTOR 轴。

因此这里提出的问题是：

> **TSC2 缺失造成的整个异常细胞状态中，还有哪些部分没有被 sirolimus 完全改变？**

以及：

> **现有药物中，有没有药物能够把这些异常状态向正常方向拉回？**

---

## 研究思路的演变

```text
TSC2-loss transcriptional state
        ↓
分别建立 plastic / hydrogel disease signature
        ↓
与 LINCS/CMap 中大量药物扰动进行比较
        ↓
寻找能够产生相反表达变化的药物
        ↓
得到候选药物
        ↓
再分析药物靶点、机制、人体 LAM 状态和其他实验
```

这里研究的不是单独恢复 `TSC2` 基因表达，而是 **TSC2 缺失后整个细胞的转录状态**。

---

## 核心数据：TSC2 × Rapamycin × Environment

GSE179044 包含：

```text
WT / TSC2-null
        ×
Vehicle / Rapamycin
        ×
Plastic / Hydrogel
```

共分析 16 个样本、59,055 个基因。

在 hydrogel 条件中，目前分类包括：

| 类型 | 基因数 |
| --- | ---: |
| 接近完全恢复 | 1,056 |
| 部分恢复但仍残留 | 3,099 |
| 持续残留 | 668 |
| 进一步恶化 | 302 |
| 方向反转 | 1,385 |

经典 mTORC1 program 本身在 residual 中已经相对较弱，更持久的信号主要涉及 ECM / invasion、myogenic programs、metabolism、autophagy 和 stress-related programs。

这使研究逐渐转向：

> **sirolimus 后剩下的问题可能不只是“mTOR 没压干净”，还可能来自其他相对独立的细胞程序。**

---

## SRPK2 / ECM 线索

目前一个值得继续追踪的方向是：

```text
TSC2 loss
   ↓
mTOR inhibition
   ↓
部分 ECM / metabolic state 仍保留
   ↓
SRPK2-sensitive program ?
```

其中 `NNMT`、`COL8A1` 等成为重要研究线索。

---

## 本地 LINCS / CMap 药物筛选

项目已经本地分析 GSE92742 和 GSE70138 两套 LINCS Level 5 数据。

基本逻辑是寻找与 LAM/TSC2-loss disease signature 方向相反的药物扰动。

目前形成：

* **258 条候选记录**；
* 去除不同数据集和 query size 的重复后得到 **66 个唯一候选药物**。

后续继续结合药物靶点、generic stress/cytotoxicity、遗传 perturbation、人体 LAM 数据和与 sirolimus 的机制互补性进行筛选。

---

## Translation program：另一个异常层次

GSE277844 用于研究 TSC2 loss 是否不仅改变 RNA abundance，也改变哪些 RNA 更容易被翻译成蛋白质。

目前识别出约：

* 89 个 translation-up genes；
* 107 个 translation-down genes；
* 合计 **196 个 translation-abnormal genes**。

在 18 个可直接比较的 hydrogel residual genes 中：

* RMC-6272：15/18 向 WT 靠近；
* eFT-508：15/18 向 WT 靠近。

两种 perturbation 共同支持 13 个基因，主要指向 ECM、cell adhesion、Rho GTPase、TGF-β 和 migration 等主题。

---

# 3. LAM Immune Visibility

📁 [LAM-Immune-Visibility](LAM-Immune-Visibility/)

## 为什么研究“免疫可见性”？

LAM 细胞具有一些特色表达，包括：

* `PMEL`
* `MLANA`
* `MITF`
* `GPNMB`
* `TYRP1`
* `DCT`

其中不少基因同时与 melanocytic lineage 和肿瘤免疫研究有关。

因此提出：

> **LAM 细胞是否已经表达了一些可以成为免疫识别线索的分子，但抗原加工、呈递或周围免疫环境并没有形成相对应的有效反应？**

项目把几个层次分开研究：

```text
抗原相关表达
       ↓
抗原加工
       ↓
HLA / antigen presentation machinery
       ↓
免疫环境
       ↓
潜在 immune visibility
```

---

## 首轮结果

目前已经生成：

* **6 个候选抗原/lineage marker 排序**；
* **735 条患者级模块结果**；
* **40 条状态关联结果**；
* **20 条免疫上下文关联结果**。

首轮候选表达情况：

| Gene | 患者中至少一次检出的比例 |
| --- | ---: |
| MITF | 100.0% |
| GPNMB | 100.0% |
| PMEL | 100.0% |
| MLANA | 84.6% |
| TYRP1 | 69.2% |
| DCT | 53.8% |

其中 `MITF`、`GPNMB` 和 `PMEL` 在当前纳入患者中具有尤其稳定的检出。

后续可继续追问：

```text
这些蛋白是否真正产生抗原肽？
            ↓
这些抗原肽是否进入 HLA presentation？
            ↓
不同患者的 HLA genotype 是否能够呈递？
            ↓
是否存在对应的 T-cell recognition？
```

---

# 4. LAM State Modeling

📁 [LAM-State-Modeling](LAM-State-Modeling/)

## 为什么单独进行 LAM 状态建模？

LAM 单细胞研究首先面临一个基础问题：

> **如何在复杂的肺组织单细胞数据中更可靠地识别 LAM 细胞？**

LAM 细胞数量较少，不同患者之间存在表达差异，同时单个标志物还会受到测序深度、掉零和细胞状态变化等因素影响。

为了尽量减少真正 LAM 细胞的遗漏，较宽的候选筛选通常具有较高召回率，但其中也可能包含具有普通肺细胞特征的细胞。

因此这一方向尝试：

> **整合多个 LAM 单细胞数据集，寻找能够跨患者重复出现的 LAM-rich 细胞群，并进一步分析这些细胞的生物学特征和内部结构。**

---

## 使用神经网络整合单细胞数据

项目使用 **scVI** 进行单细胞数据整合。

scVI 是一种基于神经网络的概率模型。

每个单细胞通常包含数千个基因的表达信息。scVI 将这些高维信息压缩为更紧凑的表示，同时处理不同数据集之间较大的技术差异。

基本过程可以概括为：

```text
每个细胞数千个基因的表达
        ↓
scVI 神经网络整合
        ↓
形成每个细胞的低维表示
        ↓
反复进行邻域和聚类分析
        ↓
建立稳定的共识状态
        ↓
验证 LAM 身份和生物学特征
```

这种方法能够利用整个转录表达模式比较细胞，而不局限于少数标志物。

---

## 总体研究思路

当前分析从上游单细胞研究得到的高置信度候选细胞池出发。

整体流程为：

```text
多个数据集中的 LAM 候选细胞
        ↓
质量控制与数据整理
        ↓
scVI 神经网络整合
        ↓
在不同参数下反复聚类
        ↓
建立共识状态
        ↓
跨患者、跨数据集稳健性分析
        ↓
重新验证各状态的 LAM 身份
        ↓
识别 LAM-rich 状态
        ↓
分析这些细胞的生物学特征
        ↓
研究不同状态之间的关系
```

分析同时考虑：

* 不同聚类参数下的稳定性；
* 不同患者和数据集之间的复现程度；
* 正式 LAMCORE 基因程序；
* 已知 LAM 相关标志物和表达程序；
* 原始研究提供的细胞注释；
* 普通肺细胞谱系信号；
* 患者级敏感性分析。

这些证据共同用于判断一个细胞群是否值得进行进一步的 LAM 生物学解释。

---

## 从五千多个候选细胞建立共识状态

当前整合分析包含来自多个患者和数据集的五千多个高置信度候选细胞。

scVI 首先将这些细胞映射到统一的低维空间。

随后项目在多组邻居数量、聚类分辨率和随机种子下重复进行聚类。

如果一批细胞在不同设置下仍然经常被分到一起，就说明这种结构具有较好的稳定性。

这些重复出现的结构最终被整理为 **共识状态**。

之后又分别进行患者和数据集的留一验证，检查去掉某一个患者或数据来源以后，相近的状态结构是否仍然能够保留。

这一过程减少了由单一患者、单一实验或单次聚类参数造成的偶然结构。

---

## 重新验证哪些状态真正富集 LAM 特征

分析进一步发现，最初的候选细胞池具有明显的生物学异质性。

不同共识状态之间的 LAM 特征强度差异较大，其中部分状态同时具有内皮、成纤维、免疫、上皮等普通肺细胞谱系特征。

因此项目增加了专门的 LAM 身份分析。

主要比较：

* 正式 LAMCORE 基因程序；
* 黑色素细胞样和其他 LAM 相关标志物；
* LAM 相关表达程序；
* 普通肺细胞谱系信号；
* 原始研究中的细胞注释；
* 不同患者中的复现情况。

这一步使研究从单纯描述不同共识状态，进一步转向识别其中最具有 LAM 生物学特征的细胞群。

---

## 识别出一个 LAM-rich 共识细胞群

在 **LAM State Modeling** 内部，一个共识状态被编号为 **State15**。

在当前所有共识状态中，State15 具有最集中的 LAM 相关证据。

主要表现为：

* 正式 LAMCORE 程序较强；
* 在具有原始 LAM 注释的数据中明显富集；
* 多种 LAM 相关标志物和表达程序较高；
* 患者匹配比较中仍然保持较强的 LAM 特征。

其中一个患者贡献了相对较多的 State15 细胞，因此又进行了去除该患者后的敏感性分析。

State15 的 LAM-rich 特征仍然能够保留。

这些结果支持将 State15 作为当前模型中最具有代表性的 **LAM-rich 共识细胞群**。

其患者分布仍不完全均衡，因此目前主要将其作为后续分析中的高置信度 LAM-rich 参考群体。

---

## 分析 LAM-rich 细胞周围的状态结构

识别出 State15 后，项目进一步研究其周围细胞是否呈现系统性的状态变化。

首先分析了从 State15 向外延伸的整体结构，并观察 LAM 相关表达是否随着距离变化。

合并分析中可以观察到一定梯度。

随后进一步检查：

* 独立 LAMCORE 指标；
* 不同患者中的变化方向；
* 不同数据集中的一致性；
* 匹配对照；
* 局部细胞邻接关系。

结果显示，LAM 相关表达在局部空间中确实存在连续变化，但不同患者和不同方向之间具有明显异质性。

因此后续分析进一步聚焦于 State15 周围的局部结构。

---

## 局部邻近状态

在 **LAM State Modeling** 中，State15 周围主要形成四个局部邻近状态：

* State16；
* State12；
* State20；
* State7。

其中 State16 表现出最清楚的患者级方向一致性。

在具有足够细胞的患者中，随着与 State15 距离增加，LAMCORE 信号均呈下降趋势；留一患者分析中也保持相同方向。

随后项目建立了更严格的局部匹配对照，同时控制：

* 患者组成；
* 数据集来源；
* 局部图距离；
* 与 State15 的连续距离。

在完成这些匹配以后，State16 的变化幅度处于相似局部结构可以出现的范围内。

同样的方法也被应用于其他邻近状态。

目前这些结果更支持将它们视为 **LAM-rich 细胞周围稳定存在的局部状态结构**，其具体生物学含义仍需要空间组学和其他独立数据进一步解释。

---

## 当前认识

这一方向目前形成了几项主要认识：

* 多数据集整合可以识别出具有较强 LAM 特征的稳定共识细胞群；
* 原始宽候选池中的 LAM 身份分布并不均一；
* 共识状态的统计稳定性和 LAM 生物学身份需要分别验证；
* LAM-rich 状态周围存在可重复的局部状态结构；
* LAM 相关表达在局部空间中呈现连续变化；
* 这些邻近状态的具体生物学含义仍有进一步研究空间。

这些结果为后续连接 LAM 细胞识别、状态异质性、空间位置和功能研究提供了一个统一框架。

---

# 四个方向之间的关系

四个项目分别研究 LAM 的不同层次：

```text
                         LAM biology
                             │
        ┌────────────────────┼────────────────────┐
        │                    │                    │
        ▼                    ▼                    ▼
   细胞程序与机制        LAM细胞识别与状态       药物反应
        │                    │                    │
 LAM细胞表达什么       哪些细胞形成稳定的       什么能够改变
 以及具有哪些功能        LAM-rich 群体          异常程序
        │                    │                    │
        └─────────────┬──────┴─────────────┬─────┘
                      │                    │
                      ▼                    ▼
                  免疫可见性            可验证机制
```

**LAM Cell Research** 主要研究 LAM 细胞中的生物学程序，包括肺适应、细胞外基质重塑、蛋白酶活动和微环境相互作用。

**LAM State Modeling** 主要研究如何在多个单细胞数据集中更稳定地识别 LAM-rich 细胞群，并描述这些细胞内部和周围的状态结构。

建立起来的状态框架可以进一步支持其他方向：

* 将 **LAM Cell Research** 中发现的基因程序定位到更稳定的 LAM-rich 细胞群；
* 在 **LAM Drug Repurposing** 中分析不同 LAM-rich 状态可能对应的药物反应；
* 在 **LAM Immune Visibility** 中比较这些细胞的抗原相关表达和抗原呈递特征。

四个方向最终可以连接为：

> **更可靠地识别 LAM 细胞，理解这些细胞的生物学程序和异质性，寻找能够改变异常程序的干预方式，并研究它们与肺组织和免疫系统之间的关系。**

---

# 目前比较值得继续追踪的线索

截至目前，比较值得继续研究的假说包括：

1. **LAM 细胞的异质性可能来自多种共享生物学程序的不同组合和连续变化，而不只表现为少数边界清晰的亚型；**
2. **更严格的跨患者 LAM 身份模型可能从当前高召回候选池中进一步分离出新的疾病特异状态；**
3. **LAM-rich 细胞周围可重复出现的局部状态可能反映 LAM 细胞与特定基质、内皮、上皮或免疫微环境之间的关系；**
4. **空间组学可能揭示这些局部状态是否在 LAM 病灶周围具有稳定的空间组织方式；**
5. **肺部可能存在由多种细胞共同形成的蛋白酶—抗蛋白酶空间生态位；**
6. **LAM 可能同时包含疾病本身的谱系/转化程序与进入肺以后获得的适应程序；**
7. **雷帕霉素治疗后可能仍保留部分细胞外基质、蛋白酶和代谢程序；**
8. **ELANE、MMP2 可能连接雷帕霉素后残留程序与肺组织破坏机制；**
9. **NNMT、COL8A1 以及 SRPK2 敏感的细胞外基质/代谢程序可能代表值得进一步验证的非 mTOR 机制；**
10. **TSC2 缺失造成的转录异常可能用于发现与 mTOR 抑制具有机制互补性的药物再利用候选；**
11. **翻译调控可能构成 TSC2 缺失异常状态的另一个重要层次；**
12. **LAM 相关谱系抗原及抗原呈递程序可能定义免疫可见性不同的 LAM 细胞状态。**

这些方向的共同目标是把大规模公开数据逐渐压缩成少量、明确、可以继续实验验证的问题。

---

# 数据来源与项目形式

项目主要使用已经公开的处理后数据，包括：

* GEO transcriptomics；
* single-cell / single-nucleus RNA-seq；
* spatial transcriptomics；
* Visium / Visium HD / Xenium；
* TSC2-loss perturbation datasets；
* rapamycin perturbation datasets；
* LINCS / CMap Level 5 perturbation data；
* 部分公开的 translation / polysome 数据。

典型分析流程为：

```text
公开数据
    ↓
标准化处理
    ↓
gene / program level analysis
    ↓
跨数据集比较
    ↓
candidate mechanism
    ↓
Hypothesis Card
```

每个子项目保存自己的分析脚本、数据 manifest、中间结果、最终表格、图表、research log 和 hypothesis cards。

---

# Repository Structure

```text
LAM-Research/
│
├── LAM-Cell-Research/
│   ├── 单细胞 / 空间组学
│   ├── LAMCORE 状态
│   ├── 基因程序发现
│   ├── 蛋白酶空间生态位
│   ├── 肺适应
│   └── 雷帕霉素后残留的 ECM / protease
│
├── LAM-Drug-Repurposing/
│   ├── TSC2 缺失 / 雷帕霉素分析
│   ├── 残留异常程序
│   ├── LINCS / CMap
│   ├── 候选药物分析
│   ├── 机制整合
│   └── 翻译调控分析
│
├── LAM-Immune-Visibility/
│   ├── 抗原相关表达
│   ├── 抗原呈递
│   ├── 免疫环境
│   ├── 候选抗原排序
│   └── 研究假说
│
└── LAM-State-Modeling/
    ├── 神经网络单细胞整合
    ├── 共识状态建模
    ├── 跨患者稳健性分析
    ├── LAM 身份分析
    └── 状态特征分析
```

详细的方法、运行方式和结果文件请进入各子目录查看。

---

# Project Status

这个仓库仍在持续演进。

当前四个方向大致处于：

| 项目 | 当前阶段 |
| --- | --- |
| LAM Cell Research | 已建立复现基线，正在推进多条生物学问题 |
| LAM Drug Repurposing | 已完成主要候选生成，正在进行机制和跨数据集验证 |
| LAM Immune Visibility | 已完成首轮计算分析，正在设计进一步验证问题 |
| LAM State Modeling | 已完成主要状态建模分析，正在进行生物学解释和与其他方向的整合 |

项目会随着新的公开 LAM 数据、已有数据的重新分析和新的研究问题继续更新。

---

## 关于结果

本仓库中的结果主要来自公开数据的计算分析，用于产生 **research hypotheses（研究假说）** 和 **candidate mechanisms（候选机制）**。

它们更适合作为后续实验研究、机制研究和药物研究的起点，不构成临床治疗建议。

---

## 许可证

本仓库中的源代码采用 [Apache License 2.0](LICENSE) 许可证发布。

项目的版权与署名信息见 [NOTICE](NOTICE)。

本项目中引用、使用或涉及的第三方数据集、软件、论文及其他材料，仍分别遵循其原始许可证和使用条款。
