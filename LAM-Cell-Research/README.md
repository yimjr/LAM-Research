# LAM Cell Research

## 本项目范围

- 公开 LAM 单细胞/单核组学数据的获取与整理
- Python/Scanpy 环境配置
- LAMCORE 细胞状态、跨供体稳定性和表达程序分析
- 可复现的分析脚本、数据清单、图表和研究报告

首阶段研究计划见 [RESEARCH_PLAN.md](RESEARCH_PLAN.md)。

## 工作区与环境
```

项目使用 Python 3.12.13 和本地虚拟环境 `.venv`。不使用 conda，也不向系统 Python 安装依赖。

```bash
./.venv/bin/python -m pip install --index-url https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple -r requirements.txt
./.venv/bin/python environment/verify_environment.py
./.venv/bin/python scripts/run_pipeline.py
```

Python 包通过清华 PyPI 镜像安装；`requirements.txt` 是依赖范围，`requirements.lock` 保存实际安装版本。

## 首阶段数据

核心复现使用 GEO `GSE135851` 的 LAM1–LAM4、Donor1，以及作者肺部对照模块中的 GSE122960 六个正常肺供体。辅助数据包括 GSE118180 wild-type 小鼠子宫。所有数据均为公开处理后矩阵，不下载 FASTQ。

```bash
./.venv/bin/python scripts/download_geo.py
./.venv/bin/python scripts/prepare_matrix.py
./.venv/bin/python scripts/qc_and_preprocess.py
./.venv/bin/python scripts/analyze_lam_states.py
```

当前计划的分阶段入口：

```bash
./.venv/bin/python scripts/build_reproduction_baseline.py
./.venv/bin/python scripts/reproduce_core.py
./.venv/bin/python scripts/run_targeted_robustness.py
./.venv/bin/python scripts/explore_lam_hypotheses.py
```

数据来源、下载哈希和运行记录分别保存在 `manifests/` 中。

## 当前状态

原有 `results/report/` 结果已归档为同一数据集上的独立再分析。新的作者风格核心重实现得到 140 个候选（LAM1/2/3/4 为 31/4/84/21），接近论文报告的约 125 个，但当前仍标记为独立重实现而非严格复现。第二阶段稳健性验证与第三阶段新发现已并行生成，详见 `results/reproduction_core/`、`results/robustness/`、`results/discovery/` 和 `results/hypothesis_cards/`。
