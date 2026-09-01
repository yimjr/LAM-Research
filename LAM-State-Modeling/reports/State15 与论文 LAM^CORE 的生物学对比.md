## State15 与原作者 LAM^CORE 的生物学对比

### 1. 基本关系

State15 与原作者定义的 LAM^CORE 在多项生物学特征上高度重合，包括：

- HOXA11、EMX2 等子宫 / 发育相关基因
- 平滑肌与收缩程序
- LAM 相关转录程序
- ECM 重塑
- 间充质样特征

State15 的主要差异基因为 **HTN3、MMP11、PAGE4、PGM5-AS1、HOXA11、EMX2、PLAT、SFRP1**，同时具有很强的 LAM-myogenic、uterine smooth muscle、CORE1、CORE3 identity 和 ECM-remodeling 程序。

原作者的 LAM^CORE 同样表现出子宫发育、平滑肌、ECM、LAM marker、激素和其他多种转录特征。

两者的高度相似说明：

**State15 和 LAM^CORE 很可能捕捉到了同一个或高度重叠的 LAM 生物学结构。**

---

### 2. 两者并不完全等价

LAM^CORE 和 State15 来自不同的数据组织方式，因此可以呈现同一底层生物学的不同切面。

#### 原作者 LAM^CORE 中较突出的特征
- PMEL、MLANA 等黑色素细胞样 LAM marker
- VEGFD 及淋巴管相关信号
- ESR1 / PGR 等激素相关特征
- HOXA11、HOXD11、EMX2 等子宫发育程序
- PCP4、UNC5D、PITX2 等 neural-like 特征
- 平滑肌 / 收缩
- ECM、迁移和分泌程序

#### State15 中较突出的特征
- HOXA11
- EMX2
- MMP11
- PLAT
- SFRP1
- 强 uterine smooth muscle
- 强 myogenic / contractile
- 明显 ECM remodeling

因此，两种表示的主要区别可以概括为：

**LAM^CORE 更明显地同时包含多种 LAM 相关生物学程序；State15 中 uterine–myogenic–ECM 这一组特征更加集中。**

这种差异本身值得研究。

---

### 3. 差异可能对应几种生物学解释

目前不能预设其中任何一种解释成立。

#### 解释 A：同一细胞群的不同表示方式

LAM^CORE 和 State15 可能主要描述同一类细胞。

不同的数据集、特征选择、降维方法和聚类方式，使某些表达程序在一个分析中更加突出。

此时：

**State15 ≈ LAM^CORE 的另一种高维表示。**

---

#### 解释 B：LAM^CORE 包含多个内部状态

LAM^CORE 可能是一个具有内部异质性的细胞群。

其中部分细胞更偏：

- uterine / HOX
- myogenic / contractile
- ECM remodeling

另一些细胞更偏：

- PCP4 / neural-like
- VEGFD / secretory
- hormone-related
- 其他 LAM 程序

此时 State15 可能对应其中一个富集区域，但其他区域同样属于 LAM biology。

---

#### 解释 C：State15 捕捉的是一条连续生物学轴

LAM 相关细胞可能沿若干连续转录轴变化，例如：

**uterine / HOX**
↕  
**myogenic / contractile**
↕  
**ECM remodeling**
↕  
**PCP4 / neural-like**
↕  
**VEGFD / secretory**

State15 代表其中 uterine–myogenic–ECM 权重较高的一端或区域。

此时离散的“LAM^CORE”和“State15”都只是对连续结构的一种划分方式。

---

#### 解释 D：两种划分分别捕捉了不同层级的生物学结构

LAM^CORE 可能更接近一种较宽的 **cell identity**。

State15 可能更接近其中一种 **cell state**。

也可能存在相反情况：State15 捕捉到一个更稳定的转录身份，而 LAM^CORE 将多个相关状态合并在一起。

目前仅靠这些转录组结果无法决定哪一种层级划分更合适。

---

### 4. 一个尤其值得关注的差异：State15 与 PCP4 相关状态

原始 LAM^CORE 中同时存在：

- HOXA11 / EMX2 等 uterine-developmental 特征
- PCP4 等 neural-like 特征

当前无监督结果中：

- **State15** 更突出 HOXA11、EMX2、MMP11
- **State16** 更突出 PCP4、HOXC10、FABP7

这一现象可以产生一个中性的研究问题：

> **HOXA11/EMX2 与 PCP4/HOXC10 是否代表 LAM 生物学中两个可部分独立变化的转录维度？**

若成立，原作者的 LAM^CORE 和当前高维 state 都是在以不同方式描述这一内部结构。

---

### 5. 当前最合适的科学表述

现阶段可以将结果表述为：

**无监督高维聚类得到的 State15 与已发表 LAM^CORE 在 uterine-developmental、myogenic/contractile 和 ECM-remodeling 等方面高度重合，但两者的转录特征并不完全一致。State15 更突出 HOXA11/EMX2–myogenic–ECM 组合，而已发表 LAM^CORE 同时覆盖更广泛的 melanocytic、hormonal、PCP4/neural-like 和 VEGFD/secretory 特征。**

**这种差异可能来源于同一细胞群的不同表示、LAM 内部异质性、连续状态轴，或不同层级的细胞身份划分。当前结果无法预设哪一种分类更接近真实生物学结构。**

---

## 6. 由此产生的核心研究问题

1. **State15 与 LAM^CORE 是否主要描述同一组底层细胞状态？**
2. **HOXA11/EMX2、PCP4、VEGFD 等程序是否可以在 LAM 细胞中独立变化？**
3. **LAM 的转录结构更接近离散亚状态，还是连续表达谱？**
4. **LAM^CORE 应理解为单一细胞身份，还是多个相关状态的集合？**
5. **不同聚类框架得到的边界差异，是否对应真实的功能差异、空间位置或疾病阶段？**

这些问题本身构成了 State15 与 LAM^CORE 对比的主要研究价值。