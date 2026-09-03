# LAM Research

一个基于公开生物医学数据的 **LAM（Lymphangioleiomyomatosis，淋巴管平滑肌瘤病）计算研究项目**。

本项目尝试重新连接公开的转录组、单细胞、空间组学、药物扰动等数据，探索与 **LAM 发病机制、肺组织破坏、药物治疗、免疫识别和细胞状态组织方式** 有关的研究问题。

目前主要包含四个相对独立、同时可以相互提供证据的研究方向：

| 方向 | 核心问题 |
| --- | --- |
| [LAM Cell Research](LAM-Cell-Research/) | LAM 细胞有哪些状态和表达程序？它们如何适应肺部环境、参与组织破坏？ |
| [LAM Drug Repurposing](LAM-Drug-Repurposing/) | 已有药物中，有没有可能逆转 LAM/TSC2-loss 的异常状态？sirolimus 没有完全改变的部分是什么？ |
| [LAM Immune Visibility](LAM-Immune-Visibility/) | LAM 细胞表达了哪些可能被免疫系统识别的特征？抗原表达、呈递和免疫环境之间有什么关系？ |
| [LAM State Modeling](LAM-State-Modeling/) | 整合单细胞数据后，能否得到稳健的 LAM-rich latent state？这些状态之间是否真的形成可复现的连续轨迹或局部转变？ |

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
* LAM 细胞具有一些比较特殊的黑色素细胞/间叶细胞相关特征，但这些特征与免疫识别之间的关系仍有探索空间；
* 整合单细胞后的 latent-state 模型还需要区分真正可复现的 LAM 结构、患者效应、普通谱系邻接和潜在空间本身的几何结构。

与此同时，过去几年已经积累了不少公开的 LAM 单细胞、空间组学、扰动实验和药物数据库。

因此，这个项目的基本想法是：

> **把已经公开、但原本服务于不同研究问题的数据重新连接起来，从中寻找新的研究问题、机制线索和可进一步实验验证的假说。**

这个仓库更关注“还能从现有数据中发现什么值得继续研究的问题”，希望最终得到能够交给实验研究者继续验证的候选机制、候选药物、研究方向，以及有价值的限制性/阴性结果。

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

## 为什么单独做状态建模？

`LAM-Cell-Research` 更关注 LAM 细胞的 gene programs、空间生态位和生物学机制；这个方向更专门研究另一个问题：

> **把多个单细胞数据整合进同一个 latent space 后，能否定义真正可复现的 LAM-rich state？这些状态之间的几何关系是否足以支持“连续轨迹”或“状态转变”这样的更强结论？**

该项目继承 `LAM-Cell-Research` 已经处理好的 AnnData、患者映射、candidate pool、LAMCORE、program 和 upstream state 信息，再建立专门的 latent-state 分析流程。

目前已经完成 **Stage 1–24**，主要包括：

* harmonization 和 QC；
* PCA/NMF baseline；
* scVI latent-space modeling；
* consensus state construction；
* patient/dataset leave-one-out robustness；
* state hierarchy 和 biology annotation；
* candidate identity audit；
* State15 anchor validation；
* global manifold 和 local branch 检验；
* matched-null、patient-level 和 LOPO robustness analysis。

---

## 当前结论

### State15 是目前最可信的 LAM-rich consensus state

在当前高召回 candidate pool 中，**State15 是最 LAM-rich 的 frozen consensus state**。

它得到 formal LAMCORE、可用数据中的 author-label enrichment、patient-matched comparison，以及去除占比较高的 LAM1163 后敏感性分析的支持。

但 State15 还不能被视为一个完全独立、患者分布均衡的正式 reference anchor。

### 目前不支持单一的 State15-centered global manifold

早期 pooled analysis 中可以看到从 State15 向外的梯度，但更严格验证后，这种结构并不能稳定表现为一个统一的 LAM manifold。

因此当前**不把 latent space 解释为一条已经建立的发育或时间轨迹**。

### State16/12/20/7 的局部分支不支持 LAM-to-lineage transition

Stage22 在 State15 周围识别出四条主要局部邻接 branch：**State16、State12、State20 和 State7**。

在修正 branch eligibility、统一 real/null 的 1–3 hop local scope、增加距离结构匹配、使用经验尾概率、进行 FDR 校正，并加入 patient-level/LOPO 检查以后，没有任何一条 branch 显示出超越 matched local null 的特异 LAM-to-lineage transition 证据。

State16 仍然是一个值得记录的稳定局部邻接状态：患者内 slope 方向一致，但这种变化并没有超过 distance-matched null 的预期，因此目前仅保留为 **ordinary lineage adjacency / mixed neighboring state**，而不是 transition candidate。

当前最稳妥的总体结论是：

> **可以识别一个可复现的 LAM-rich 核心状态，但现有数据既没有建立一条统一的 global LAM trajectory，也没有建立一个特异的局部 LAM-to-lineage transition。**

项目也不声称已经证明时间转化、建立诊断 classifier，或证明所有 candidate state 都是真正的 LAM。

---

# 四个方向之间的关系

四个方向可以独立研究，但也从不同层次观察同一种疾病：

```text
                         LAM biology
                             │
       ┌─────────────────────┼─────────────────────┐
       │                     │                     │
       ▼                     ▼                     ▼
  Cell Programs         State Geometry        Drug Response
       │                     │                     │
细胞表达什么、如何互作   哪些 latent state      什么能改变
                        真正可复现              异常程序
       │                     │                     │
       └──────────────┬──────┴──────────────┬──────┘
                      │                     │
                      ▼                     ▼
             Immune Visibility      Testable Hypotheses
                      │
                 免疫系统能看到什么
```

`LAM-Cell-Research` 与 `LAM-State-Modeling` 尤其互补：

* **LAM Cell Research** 更强调 gene programs、空间生态位、lung adaptation 和生物学机制；
* **LAM State Modeling** 更强调 integrated latent states 是否足够稳健，以及 state-to-state geometry 是否真的支持更强的结构或 transition 假说。

例如，在 **LAM Cell Research** 中发现的 rapamycin-persistent ECM program，可以进一步：

* 映射到 **LAM State Modeling** 中较稳健的 state；
* 在 **Drug Repurposing** 中寻找能够逆转它的药物；
* 在空间数据中观察它是否位于肺组织破坏区域；
* 在 **Immune Visibility** 中研究这种状态是否伴随不同的免疫环境。

四个方向最终逐渐汇总成一个共同框架：

> **识别可复现的 LAM 状态和程序，寻找能够改变这些异常的因素，并进一步理解这些状态如何与肺组织和免疫系统相互作用。**

---

# 目前比较值得继续追踪的线索

1. **LAMCORE 异质性可能包含连续和局部变化，但当前状态建模不支持单一、稳健的 global LAM trajectory；**
2. **State15 是目前最强的 LAM-rich frozen consensus state，但患者分布仍限制其作为正式 reference anchor 的强度；**
3. **State16 是 State15 的可复现局部邻接状态，但 distance-matched null 不支持把它解释为 LAM-to-lineage transition；**
4. **肺部可能存在由多种细胞共同形成的 protease–antiprotease spatial niche；**
5. **LAM 可能同时包含疾病本身的 lineage/transformation program 与进入肺后获得的 lung-adaptation program；**
6. **rapamycin 后仍可能保留部分 ECM / protease / metabolic programs；**
7. **ELANE、MMP2 等可能连接 rapamycin persistence 与肺组织破坏研究；**
8. **NNMT、COL8A1 及 SRPK2-sensitive ECM/metabolic program 值得继续研究；**
9. **TSC2-loss transcriptional state 已经产生一批可继续筛选的药物再利用候选；**
10. **translation regulation 可能是 TSC2-loss abnormal state 的另一个重要层次；**
11. **PMEL、MITF、GPNMB 等 LAM-associated lineage signals 在多个患者中具有较稳定表达；**
12. **LAM antigen expression、antigen presentation 与 immune context 之间的关系值得进一步连接研究。**

这些方向的共同目标不是提前确定一个答案，而是把大规模公开数据逐渐压缩成少量、明确、可以继续实验验证的问题。

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
gene / program / latent-state analysis
    ↓
跨数据集比较和稳健性检验
    ↓
candidate mechanism 或限制性结论
    ↓
可实验验证的问题
```

每个子项目保存自己的分析脚本、数据 manifest、中间和最终结果、结果表格、图表以及研究记录/报告。

---

# Repository Structure

```text
LAM-Research/
│
├── LAM-Cell-Research/
│   ├── 单细胞 / 空间组学
│   ├── LAMCORE 状态
│   ├── gene program discovery
│   ├── protease spatial niche
│   ├── lung adaptation
│   └── rapamycin-persistent ECM/protease
│
├── LAM-Drug-Repurposing/
│   ├── TSC2-loss / rapamycin factorial analysis
│   ├── residual programs
│   ├── LINCS / CMap
│   ├── drug candidate analysis
│   ├── mechanism integration
│   └── translation analysis
│
├── LAM-Immune-Visibility/
│   ├── antigen-related expression
│   ├── antigen presentation
│   ├── immune context
│   ├── candidate antigen ranking
│   └── hypothesis cards
│
└── LAM-State-Modeling/
    ├── scVI latent-state modeling
    ├── consensus 与 robustness analysis
    ├── candidate-identity audits
    ├── State15 anchor validation
    ├── manifold / local-branch testing
    └── final state-modeling reports
```

详细的方法、运行方式和结果文件请进入各子目录查看。

---

# Project Status

这个仓库仍在持续演进。

当前四个方向大致处于：

| 项目 | 当前阶段 |
| --- | --- |
| LAM Cell Research | 已建立复现基线，进入多条新生物学问题探索 |
| LAM Drug Repurposing | 已完成主要 TSC2-loss/LINCS 候选生成，进入机制与跨数据验证 |
| LAM Immune Visibility | 已完成首轮计算和候选抗原排序，进入进一步验证问题设计 |
| LAM State Modeling | Stage 1–24 已完成；State15 保留为最强 LAM-rich consensus state，而 global manifold 与 local transition 均未得到最终 robustness analysis 支持 |

项目会随着新的公开 LAM 数据、已有数据的重新分析和新的研究问题继续更新。

---

## 关于结果

本仓库中的结果主要来自公开数据的计算分析，用于产生 **research hypotheses（研究假说）**、**candidate mechanisms（候选机制）**，以及定义清楚的阴性或限制性结果。

它们更适合作为后续实验研究、机制研究和药物研究的起点，不构成临床治疗建议。

---

## 许可证

本仓库中的源代码采用 [Apache License 2.0](LICENSE) 许可证发布。

项目的版权与署名信息见 [NOTICE](NOTICE)。

本项目中引用、使用或涉及的第三方数据集、软件、论文及其他材料，仍分别遵循其原始许可证和使用条款。
