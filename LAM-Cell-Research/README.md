# LAM Cell Research

## 本项目范围

- 公开 LAM 单细胞/单核组学数据的获取与整理
- Python/Scanpy 环境配置
- LAMCORE 细胞状态、跨供体稳定性和表达程序分析
- 可复现的分析脚本、数据清单、图表和研究报告

完整研究计划见 [RESEARCH_PLAN.md](RESEARCH_PLAN.md)。

## 项目阶段与当前重点

项目按三个层次组织：先进行原论文核心结果复现，再进行必要的稳健性检验，最后进入以新生物学问题为中心的探索阶段。前两步用于建立可信基线；当前工作的重点已经转向探索，不再以重复方法比较作为主要目标。对于探索中出现的具体候选，再按需补充验证。

目前探索重点包括：LAMCORE 状态异质性、肺部 protease–antiprotease 空间生态位、rapamycin 后仍保留的 ECM/protease 程序、LAM 细胞的肺适应程序，以及肺内状态与血浆/EV 蛋白之间的联系。

## 工作区与环境

项目使用 Python 3.12.13 和本地虚拟环境 `.venv`。不使用 conda，也不向系统 Python 安装依赖。

```bash
./.venv/bin/python -m pip install --index-url https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple -r requirements.txt
./.venv/bin/python environment/verify_environment.py
./.venv/bin/python scripts/run_pipeline.py
```

Python 包通过清华 PyPI 镜像安装；`requirements.txt` 是依赖范围，`requirements.lock` 保存实际安装版本。

## 数据与分析入口

核心复现使用 GEO `GSE135851` 的 LAM1–LAM4、Donor1，以及作者肺部对照模块中的 GSE122960 六个正常肺供体。辅助数据包括 GSE118180 wild-type 小鼠子宫。所有数据均为公开处理后矩阵，不下载 FASTQ。

```bash
./.venv/bin/python scripts/download_geo.py
./.venv/bin/python scripts/prepare_matrix.py
./.venv/bin/python scripts/qc_and_preprocess.py
./.venv/bin/python scripts/analyze_lam_states.py
```

按“复现与稳健性基础 → 新生物学探索”的顺序运行：

```bash
./.venv/bin/python scripts/build_reproduction_baseline.py
./.venv/bin/python scripts/reproduce_core.py
./.venv/bin/python scripts/run_targeted_robustness.py
./.venv/bin/python scripts/explore_lam_hypotheses.py
# HDBSCAN 是状态异质性探索的独立方法分支，不覆盖 KMeans 结果
./.venv/bin/python scripts/explore_lam_hypotheses_hdbscan.py
```

数据来源、下载哈希和运行记录分别保存在 `manifests/` 中。

## 当前状态

论文核心复现和必要稳健性检验已经形成当前分析基线。新的作者风格重实现得到 140 个候选（LAM1/2/3/4 为 31/4/84/21），接近论文报告的约 125 个，但仍标记为独立重实现而非严格复现。

当前项目重点是基于这条基线开展新生物学探索：空间 protease 来源与 proteolytic balance、rapamycin 后 ECM/protease 保留、四组肺适应 interaction，以及 pooled/donor-wise 程序发现。LAMCORE 状态异质性还增加了 HDBSCAN 独立分支；当前主设定在 140 个候选中识别出 2 个局部密度簇、28 个非噪声细胞，未支持覆盖四个 donor 的稳定离散亚型。详见 `results/discovery/`、`results/discovery_hdbscan/`、`results/spatial/`、`results/perturbation/`、`results/adaptation/` 和相应的 `results/hypothesis_cards*/`。
