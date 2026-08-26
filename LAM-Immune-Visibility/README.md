# LAMCORE Immune Visibility

独立的 LAMCORE 免疫可见性计算研究目录。

本项目只读复用 `../LAM-Cell-Research` 和
`../LAM-Drug-Repurposing` 中的公开数据与已有结果，
所有新结果写入本目录。首轮不下载新数据；数据缺口在报告中标记，后续再决定是否补充。

## 研究问题

LAMCORE 是否存在一种细胞状态：抗原相关表达已经检出，但抗原加工/呈递不足，或免疫逃逸程序增强？

本项目严格区分：

- `not_assayed`：当前矩阵或空间 panel 没有该基因；
- `not_detected`：基因在矩阵中，但当前细胞/spot 原始 count 为 0；
- `detected_low`：有非零 count，已经检出，但处于预先定义的低表达区间；
- `detected_high`：有非零 count，已经检出，且处于预先定义的高表达区间。

`not_detected` 不解释为真实生物学不表达，`detected_low` 不解释为未检出。

## 运行

使用已有 LAM-Cell-Research 环境：

```bash
../LAM-Cell-Research/.venv/bin/python \
  scripts/run_visibility_pipeline.py --stage all
```

仅运行主细胞评分：

```bash
../LAM-Cell-Research/.venv/bin/python \
  scripts/run_visibility_pipeline.py --stage score
```

## 解释边界

- antigen module 是抗原相关表达，不等同于免疫原性或真实 HLA 肽呈递；
- presentation module 是呈递机器表达，不等同于 immunopeptidomics 证据；
- spatial/ligand–receptor 结果是候选关联，不是已证明通信；
- rapamycin 结果是扰动支持的候选机制，不是患者级治疗因果；
- 只有至少两个独立 PatientID、至少两个研究体系方向一致的结果，才升级为高优先级候选。
