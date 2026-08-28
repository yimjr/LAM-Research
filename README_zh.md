# LAM Research

一个基于公开生物医学数据的 **LAM（Lymphangioleiomyomatosis，淋巴管平滑肌瘤病）计算研究项目**。

本项目尝试从公开的转录组、单细胞、空间组学、药物扰动等数据中重新提出问题，寻找与 **LAM 发病机制、肺组织破坏、药物治疗和免疫识别** 有关的研究线索。

目前主要包含三个相对独立、同时可以相互提供证据的研究方向：

| 方向                                              | 核心问题                                                      |
| ----------------------------------------------- | --------------------------------------------------------- |
| [LAM Cell Research](LAM-Cell-Research/)         | LAM 细胞究竟有哪些状态？它们如何适应肺部环境、参与组织破坏？                          |
| [LAM Drug Repurposing](LAM-Drug-Repurposing/)   | 已有药物中，有没有可能逆转 LAM/TSC2-loss 的异常状态？sirolimus 没有完全改变的部分是什么？ |
| [LAM Immune Visibility](LAM-Immune-Visibility/) | LAM 细胞表达了哪些可能被免疫系统识别的特征？免疫识别与抗原呈递之间是否存在值得研究的差异？           |

[English](README.md)

---

## 为什么做这个项目？

LAM 是一种罕见病。一个重要的分子基础是 **TSC1/TSC2 功能异常及其导致的 mTOR 通路持续活跃**。

目前 sirolimus（西罗莫司/雷帕霉素）已经显著改变了 LAM 的治疗，但很多问题仍然值得继续研究：

* LAM 细胞并不是完全相同的一群细胞，它们可能具有不同状态；
* mTOR 抑制能够控制很多异常，但部分细胞程序仍可能保留下来；
* LAM 为什么尤其容易在肺内形成破坏性病变，仍有很多机制可以继续探索；
* 肺囊性破坏可能涉及 LAM 细胞、成纤维细胞、免疫细胞和蛋白酶系统共同作用；
* 已有药物中可能存在能够作用于 mTOR 之外异常程序的候选；
* LAM 细胞具有一些比较特殊的黑色素细胞/间叶细胞相关特征，但这些特征与免疫识别之间的关系仍有探索空间。

与此同时，过去几年已经积累了不少公开的 LAM 单细胞、空间组学、扰动实验和药物数据库。

因此，这个项目的基本想法是：

> **把已经公开、但原本服务于不同研究问题的数据重新连接起来，从中寻找新的研究问题、机制线索和可进一步实验验证的假说。**

这个仓库更关注“还能从现有数据中发现什么值得继续研究的问题”，希望最终得到的是能够交给实验研究者继续验证的候选机制、候选药物和研究方向。

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

大致流程为：

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

> 一组经常一起变化、可能共同反映某种细胞功能的基因。

这样能够避免只围绕某一个 marker 展开研究，而是从整个细胞状态出发寻找规律。

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

这使一个研究方向逐渐变得有趣：

> **LAMCORE 内部的异质性可能更多体现为连续变化的细胞状态，而不是几个边界非常清楚的新细胞亚型。**

项目因此也在分析连续状态、局部状态以及 donor-specific state，而不只寻找新的 cluster。

---

### 2. 寻找尚未被充分描述的 LAM 表达程序

项目使用 pooled NMF、donor-wise NMF、meta-program matching 等方式寻找在多个细胞中共同出现的表达程序。

首轮已经得到一些候选程序，其中涉及：

* ECM / extracellular matrix；
* protease；
* migration；
* proliferation；
* hormone-related signals；
* LAM lineage；
* microenvironment interaction。

目前没有出现一个非常明显、能够在多个独立患者中重新发现并可以直接定义为“全新 LAMCORE 程序”的结果。

这本身也提供了一条线索：

> LAM 的重要差异可能来自已知程序的重新组合、强弱变化和微环境依赖，而不一定表现为一个完全独立的新细胞类型。

部分候选程序仍在通过其他 LAM 数据集继续比较。

---

### 3. Protease–antiprotease 空间生态位

LAM 最重要的病理特征之一是肺组织逐渐形成大量囊性破坏。

因此项目进一步研究：

> **破坏肺组织的蛋白酶信号究竟来自 LAM 细胞本身，还是 LAMCORE、成纤维细胞、免疫细胞等多种细胞共同形成一个 proteolytic niche（蛋白水解生态位）？**

目前已经在：

* Visium；
* Visium HD；
* Xenium

等空间数据中观察到：

**LAMCORE-like spatial signal 与 protease signal 存在方向一致的空间关联。**

目前关注的蛋白酶/相关基因包括：

* `CTSK`
* `MMP` family
* `ELANE`
* `PRTN3`
* `CTSS`

同时也分析 antiprotease，从而形成：

```text
protease activity
        -
antiprotease activity
        ↓
proteolytic balance
```

这一方向希望最终解释：

> LAM 肺囊性破坏是否来自一个由多种细胞共同维持的局部蛋白水解环境。

---

### 4. LAM 细胞进入肺后是否获得新的适应程序？

LAM 具有一个很特别的问题：

> LAM 细胞可能来自肺外组织，但最终能够在肺内长期生存并形成病灶。

因此项目比较：

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

目前得到一个有趣的初步结果：

* `lung_adaptation` 程序更符合 **进入肺后获得的状态**；
* ECM 程序在 LAM uterus 和 pulmonary LAM 中都增强。

这意味着两类信号可能需要区分：

```text
LAM transformation / lineage program
        +
lung-acquired adaptation program
        =
pulmonary LAM state
```

因此，“LAM 为什么特别能够适应肺”成为这个方向继续研究的问题之一。

---

### 5. Rapamycin 后仍保留的 ECM / protease 程序

sirolimus/rapamycin 能够很好地抑制 mTOR 相关生长信号。

项目进一步追问：

> **如果细胞生长受到控制，与组织破坏、ECM 或蛋白酶有关的程序是不是全部同时消失？**

在 TSC2-loss 扰动数据中，目前已经观察到部分 ECM/protease 相关基因在 rapamycin 后仍存在一定程度的异常。

其中：

* `ELANE` 在 plastic 和 hydrogel 两种环境下均表现出方向一致的 partial retention；
* `MMP2` 在 hydrogel 环境中表现出 retention 信号。

这支持继续研究一个潜在机制：

> **mTOR inhibition 可能很好地控制细胞生长，同时仍有部分 matrix-related pathology 值得单独研究。**

这条线也与药物再利用项目直接产生了联系。

---

# 2. LAM Drug Repurposing

📁 [LAM-Drug-Repurposing](LAM-Drug-Repurposing/)

## 为什么研究药物再利用？

LAM 已经存在有效治疗药物 sirolimus。

但 sirolimus 主要针对 mTOR 轴。

因此这里提出的问题是：

> **TSC2 缺失造成的整个异常细胞状态中，还有哪些部分没有被 sirolimus 完全改变？**

以及：

> **现有药物中，有没有药物能够把这些异常状态向正常方向拉回？**

如果能够找到这样的药物，它们可能成为：

* 新的机制研究工具；
* sirolimus 联合治疗研究的候选；
* 已有药物再利用的研究方向。

---

## 研究思路的演变

这个项目最初从：

```text
TSC2 loss
    ↓
rapamycin
    ↓
还有哪些异常没有恢复？
```

开始。

后来分析发现，plastic 与 hydrogel 两种环境之间的 residual 并不完全一致。

因此目前使用更宽的策略：

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

这里研究的不是“TSC2 这个基因表达量能否恢复”。

使用的是：

> **TSC2 缺失以后整个细胞发生的转录变化。**

因此一个药物即使不直接作用 TSC2，也可能从下游把部分异常状态向 WT 靠近。

---

## 核心数据：TSC2 × Rapamycin × Environment

GSE179044 提供了一个非常适合这个问题的数据结构：

```text
WT / TSC2-null
        ×
Vehicle / Rapamycin
        ×
Plastic / Hydrogel
```

共分析：

* 16 个样本；
* 59,055 个基因。

这样可以分别观察：

1. TSC2 缺失造成了什么；
2. rapamycin 修复了什么；
3. rapamycin 后还剩下什么；
4. extracellular environment 是否影响这些变化。

---

## Rapamycin 后的表达状态

在 hydrogel 条件中，目前得到的分类包括：

| 类型       |   基因数 |
| -------- | ----: |
| 接近完全恢复   | 1,056 |
| 部分恢复但仍残留 | 3,099 |
| 持续残留     |   668 |
| 进一步恶化    |   302 |
| 方向反转     | 1,385 |

一个比较有意思的观察是：

经典 **mTORC1 program 本身在 residual 中已经很弱**。

相对更容易留下来的信号出现在：

* ECM / invasion；
* myogenic programs；
* metabolism；
* autophagy；
* stress-related programs。

这使后续研究逐渐转向：

> **sirolimus 后剩下的问题可能不只是“mTOR 没压干净”，还可能来自其他相对独立的细胞程序。**

在 hydrogel 相关分析中还出现了：

* `NNMT`
* `COL8A1`
* `MIR210HG`
* `SLC40A1`
* `FBLN5`
* `DCN`
* `LUM`

等值得进一步研究的信号。

---

## SRPK2 / ECM 线索

另外的数据集进一步提供了一个很有意思的机制连接。

在相关细胞模型中，rapamycin 后：

* `COL8A1`
* `NNMT`
* `DCN`
* `ACTA2`
* `MMP2`

等基因增加。

而 SRPK2 knockdown 可以降低包括：

* `NNMT`
* `COL3A1`
* `LUM`
* `FBLN5`

在内的一部分 ECM / metabolic signals。

因此形成了一个目前比较值得继续追踪的方向：

```text
TSC2 loss
   ↓
mTOR inhibition
   ↓
部分 ECM / metabolic state 仍保留
   ↓
SRPK2-sensitive program ?
```

其中 `NNMT / COL8A1` 等成为比较重要的研究线索。

---

## 本地 LINCS / CMap 药物筛选

项目已经下载并本地分析：

* GSE92742
* GSE70138

两套 LINCS Level 5 数据。

基本逻辑非常直观：

```text
LAM / TSC2-loss：A↑ B↑ C↓ D↓

某种药物：       A↓ B↓ C↑ D↑

                ↓

这个药物可能具有 reversal potential
```

目前候选生成不要求一个药物必须同时在 plastic 和 hydrogel 两种环境下成功。

只要在其中一个具有可信的 TSC2-loss reversal signal，就可以进入后续考虑。

目前：

* 形成 **258 条候选记录**；
* 去除不同数据集和 query size 的重复后；
* 得到 **66 个唯一候选药物**。

这些药物随后继续按照：

* 药物作用靶点；
* 是否主要产生 generic stress / cytotoxicity；
* 遗传 perturbation 是否支持相同机制；
* LAM 人体数据中靶点和相关程序是否存在；
* 与 sirolimus 是否可能具有机制互补性

进行分析。

候选中既可以看到 mTOR/PI3K/AKT 相关方向，也出现了值得进一步研究的非 mTOR 机制，因此后续重点并不局限于再次寻找另一种 mTOR inhibitor。

---

## Translation program：另一个异常层次

项目还加入了 GSE277844，用来研究：

> TSC2 loss 是否除了改变“产生多少 RNA”，还改变“哪些 RNA 更容易被翻译成蛋白质”。

目前识别出约：

* 89 个 translation-up genes；
* 107 个 translation-down genes；

合计 **196 个 translation-abnormal genes**。

其中一个值得继续研究的现象是：

与 hydrogel residual 重叠的一组基因，对两个 translation-targeting perturbations：

* RMC-6272；
* eFT-508

都出现较明显的向 WT 靠近趋势。

在 18 个可以直接比较的 hydrogel residual genes 中：

* RMC-6272：15/18 向 WT 靠近；
* eFT-508：15/18 向 WT 靠近。

其中 **13 个基因被两种药物共同支持**：

```text
CDC42EP3
RND3
SERPINE2
GPR27
WWTR1
FBN2
REEP2
ZNF354C
PNMA2
CACFD1
NFATC4
SPIN4
GPC4
```

这些基因进一步指向几个可能相关的主题：

* ECM；
* cell adhesion；
* Rho GTPase；
* TGF-β；
* migration。

因此 translation regulation 也成为解释 TSC2-loss residual state 的一个潜在方向。

---

# 3. LAM Immune Visibility

📁 [LAM-Immune-Visibility](LAM-Immune-Visibility/)

## 为什么研究“免疫可见性”？

LAM 细胞具有一些非常有特色的表达特征。

例如：

* `PMEL`
* `MLANA`
* `MITF`
* `GPNMB`
* `TYRP1`
* `DCT`

其中不少基因同时与 melanocytic lineage 和肿瘤免疫研究有关。

这产生了一个很自然的问题：

> **LAM 细胞是否已经表达了一些可以成为免疫识别线索的分子，但抗原加工、呈递或周围免疫环境并没有形成相对应的有效反应？**

因此这个项目没有简单寻找“哪个基因表达最高”，而是把几个步骤分开研究：

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

## 基本研究方法

这一项目复用了前两个项目已经整理的公开单细胞和空间数据。

分别分析：

### Antigen-related expression

LAM 细胞是否表达潜在抗原相关基因。

### Antigen presentation

细胞是否同时表达：

* HLA；
* antigen processing；
* presentation machinery

相关基因。

### Immune context

观察候选状态周围是否同时存在：

* T cell；
* NK；
* macrophage；
* immune suppression

等相关信号。

### Treatment persistence

进一步观察部分特征在 rapamycin perturbation 后是否仍然保留。

---

## 首轮结果

目前已经生成：

* **6 个候选抗原/lineage marker 排序**；
* **735 条患者级模块结果**；
* **40 条状态关联结果**；
* **20 条免疫上下文关联结果**。

首轮候选表达情况如下：

| Gene  | 患者中至少一次检出的比例 |
| ----- | -----------: |
| MITF  |       100.0% |
| GPNMB |       100.0% |
| PMEL  |       100.0% |
| MLANA |        84.6% |
| TYRP1 |        69.2% |
| DCT   |        53.8% |

其中 `MITF`、`GPNMB` 和 `PMEL` 在当前纳入患者中具有尤其稳定的检出。

这说明 melanocytic / lineage-related signals 在 LAM 中并不是个别细胞偶然出现的现象。

因此后续可以继续追问：

```text
这些蛋白是否真正产生抗原肽？
            ↓
这些抗原肽是否进入 HLA presentation？
            ↓
不同患者的 HLA genotype 是否能够呈递？
            ↓
是否存在对应的 T-cell recognition？
```

这也使 `PMEL/gp100` 等已有免疫研究基础的抗原成为很有价值的参考点，同时可以继续寻找新的 LAM-associated antigen candidates。

---

# 三个方向之间的关系

这三个项目可以独立研究。

同时，它们也刚好从三个不同角度观察同一个问题。

```text
                LAM biology
                    │
       ┌────────────┼────────────┐
       │            │            │
       ▼            ▼            ▼
   Cell State    Drug Response   Immune Visibility
       │            │            │
 LAM细胞在做什么   什么能改变它   免疫系统能看到什么
       │            │            │
       └────────────┼────────────┘
                    ▼
             Testable Hypotheses
```

例如：

一个在 **Cell Research** 中发现的 rapamycin-persistent ECM program，可以进一步：

* 在 **Drug Repurposing** 中寻找能够逆转它的药物；
* 在空间数据中观察它是否位于肺组织破坏区域；
* 在 **Immune Visibility** 中研究这种状态是否同时伴随不同的免疫环境。

反过来，药物筛选发现的候选机制，也可以回到单细胞数据中定位：

> 到底是哪一种 LAM 状态最可能受到这个药物影响？

因此三个方向最终可以逐渐形成一个共同的研究框架：

> **从“细胞处于什么状态”，走向“什么能够改变这个状态”，再进一步研究这些状态如何与肺组织和免疫系统相互作用。**

---

# 目前比较值得继续追踪的线索

截至目前，几个尤其值得继续推进的问题包括：

1. **LAMCORE 可能具有明显的连续状态结构，而不是只有少数固定亚型；**
2. **肺部可能存在由多种细胞共同形成的 protease–antiprotease spatial niche；**
3. **LAM 可能同时包含疾病本身的 lineage/transformation program 与进入肺后获得的 lung-adaptation program；**
4. **rapamycin 后仍可能保留部分 ECM / protease / metabolic programs；**
5. **ELANE、MMP2 等可能连接 rapamycin persistence 与肺组织破坏研究；**
6. **NNMT、COL8A1 及 SRPK2-sensitive ECM/metabolic program 值得继续研究；**
7. **TSC2-loss transcriptional state 已经产生一批可继续筛选的药物再利用候选；**
8. **translation regulation 可能是 TSC2-loss abnormal state 的另一个重要层次；**
9. **PMEL、MITF、GPNMB 等 LAM-associated lineage signals 在多个患者中具有较稳定的表达；**
10. **LAM antigen expression、antigen presentation 与 immune context 之间的关系值得进一步连接研究。**

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

分析结果通常整理为：

```text
原始公开数据
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

每个子项目保存自己的：

* 分析脚本；
* 数据 manifest；
* 中间结果；
* 最终表格；
* 图表；
* research log；
* hypothesis cards。

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
└── LAM-Immune-Visibility/
    ├── antigen-related expression
    ├── antigen presentation
    ├── immune context
    ├── candidate antigen ranking
    └── hypothesis cards
```

详细的方法、运行方式和结果文件请进入各子目录查看。

---

# Project Status

这个仓库仍在持续演进。

当前三个方向大致处于：

| 项目                    | 当前阶段                                  |
| --------------------- | ------------------------------------- |
| LAM Cell Research     | 已建立复现基线，进入多条新生物学问题探索                  |
| LAM Drug Repurposing  | 已完成主要 TSC2-loss/LINCS 候选生成，进入机制与跨数据验证 |
| LAM Immune Visibility | 已完成首轮计算和候选抗原排序，进入进一步验证问题设计            |

项目会随着新的公开 LAM 数据、已有数据的重新分析和新的研究问题继续更新。

---

## 关于结果

本仓库中的结果主要来自公开数据的计算分析，用于产生 **research hypotheses（研究假说）和 candidate mechanisms（候选机制）**。

它们更适合作为后续实验研究、机制研究和药物研究的起点，不构成临床治疗建议。

---

## 许可证

本仓库中的源代码采用 [Apache License 2.0](LICENSE) 许可证发布。

项目的版权与署名信息见 [NOTICE](NOTICE)。

本项目中引用、使用或涉及的第三方数据集、软件、论文及其他材料，仍分别遵循其原始许可证和使用条款。

