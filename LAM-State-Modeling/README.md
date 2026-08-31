# LAM-State-Modeling

本项目实现 LAM latent-state 研究的阶段 1–6 和 Step 7–13。它与 `LAM-Cell-Research` 并行，优先继承原项目的已转换 AnnData、patient/donor mapping、candidate pool、LAMCORE、program 和既有状态结果，不重新执行 GEO 转换或 program discovery。

## 运行顺序

```bash
PY=/mnt/py-env/venvs/LAM-State-Modeling/bin/python
$PY 01_inventory_inputs.py --config config/state_modeling.yaml
$PY 02_inherit_and_prepare.py --config config/state_modeling.yaml
$PY 03_qc_and_harmonize.py --config config/state_modeling.yaml
$PY 04_baseline_pca_nmf.py --config config/state_modeling.yaml
$PY 05_train_scvi.py --config config/state_modeling.yaml
$PY 06_stage6_checkpoint.py --config config/state_modeling.yaml
$PY 07_consensus_stability.py --config config/state_modeling.yaml
$PY 08_loo_robustness.py --config config/state_modeling.yaml
$PY 09_state_hierarchy.py --config config/state_modeling.yaml
$PY 10_biology_annotation.py --config config/state_modeling.yaml
$PY 11_patient_reproducibility.py --config config/state_modeling.yaml
$PY 12_boundary_normal_validation.py --config config/state_modeling.yaml
$PY 13_state_atlas.py --config config/state_modeling.yaml
$PY 14_merge_consensus_upstream.py --config config/state_modeling.yaml
$PY 16_rebuild_lam_identity_gate.py --config config/state_modeling.yaml
$PY 17_identity_calibration_audit.py --config config/state_modeling.yaml
$PY 18_validate_state15_anchor.py --config config/state_modeling.yaml --block-size 8192
$PY 19_state15_cross_patient_audit.py --config config/state_modeling.yaml --block-size 8192
$PY 20_state15_centered_manifold.py --config config/state_modeling.yaml --block-size 4096
$PY 21_validate_state15_manifold.py --config config/state_modeling.yaml --block-size 4096
$PY 22_state15_local_branch_analysis.py --config config/state_modeling.yaml --block-size 4096
$PY 23_visualize_state15_latent_space.py
```

虚拟环境位置固定为仓库外的 `/mnt/py-env/venvs/LAM-State-Modeling`，以减少 WSL 跨文件系统磁盘开销；仓库内不再使用 `.venv`。依赖安装建议使用清华镜像，PyTorch CUDA wheel 单独从官方 CUDA 12.8 源安装：

```bash
$PY -m pip install --index-url https://download.pytorch.org/whl/cu128 'torch==2.8.0'
$PY -m pip install -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple -r environment/requirements.txt
$PY -m pip install -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple -r environment/requirements-dev.txt
```

项目根目录的 `pyrightconfig.json` 已按实际虚拟环境配置：Python 3.12、`/mnt/py-env/venvs/LAM-State-Modeling`、`basic` 类型检查，并排除数据、结果、报告和缓存目录。安装 Pyright 后在项目根目录运行：

```bash
$PY -m pyright
```

阶段 5 的 `accelerator=auto` 会在 `torch.cuda.is_available()` 为真时使用 CUDA，否则记录并回退 CPU。大 AnnData 文件请逐个脚本运行，不要并行启动多个阶段。

若只需验证 GPU/流程，可运行 `05_train_scvi.py --max-epochs 1`；正式 Stage 6 结论应运行配置中的完整 200 epochs。本次正式结果已完成 200 epochs，并由 Stage 6 checkpoint 重新评估。Step 7–13 只复用既有 `data/processed/state_model_scvi.h5ad` 中的 `X_scVI`，不会重训 scVI。

默认输入根目录依次为 sibling `LAM-Cell-Research`、`data-temp` 和本项目的 `data/upstream`。大文件和结果目录均被 gitignore；输入清单和校验结果写入 `results/stage1_6/`。

核心 candidate 固定为 `pool_high_confidence`。`pool_broad_lam_like - pool_high_confidence` 只作为 boundary，`pool_unrestricted_lam` 只作为审计标签。

Stage 6 会保留阶段 5 在 candidate+boundary+normal cohort 上生成的全体 scVI Leiden 标签，但不使用它们定义 LAM 状态；正式 Go/No-Go 只基于已有 `X_scVI` 中 high-confidence candidate 的 LAM-only neighbors/Leiden。Stage 6 使用自己专属的 `lam_n_neighbors`、`lam_leiden_resolution`，并运行 15/30/50 × 0.2/0.4/0.6 小网格；不以 cluster 数量单独选择参数，而是比较碎片化指标、患者覆盖和跨参数 partition ARI。boundary 和 normal 仅产生辅助连接/邻域结果。

NMF 使用 `counts → normalize_total → log1p → HVG → NMF`；scVI 只使用 `layers["counts"]`，并且第一版仅使用 `batch_key="dataset"`，不使用 `assay` categorical covariate。

Step 7 使用 9 个 grid configuration 等权；其中 `n_neighbors=30` 的三个 configuration 各自先在 5 个 seed 内平均，因此 21 个 raw partitions 不会等权。最终允许生成一份 5,378×5,378 float32 co-assignment，并对完整距离执行 average-linkage。Step 8 的 LOO overlap 只在 retained cells 上计算，并保存 full 30/0.4 reference→consensus baseline。

Step 10 对每个 consensus state 单独执行 patient × group pseudobulk：`State_k vs Rest_of_LAM`，设计为 `~ patient_id + group`；低患者支持的 state 只输出描述性结果。Step 12 的 boundary/normal 只作辅助验证。Step 13 不设置硬性的 high/medium/low confidence 门槛。

Stage 16 是独立的 candidate identity gate 重构，不要求先运行或重跑 Step 7–13。它读取 `state_model_prepared.h5ad` 的全部 `condition=LAM` 细胞，使用 PMEL/MLANA/MITF、LAMCORE/CORE2/CORE3 连续证据作为 identity anchors，ACTA2/ESR1/VEGFD/CTSK 作为 supportive evidence，并加入 ciliated、AT2、myeloid、endothelial、fibroblast、mesothelial 和 conditional pericyte/VSMC competing-lineage scores。`FIGF` 在证据计算中统一为 `VEGFD`；不使用“任意两个 marker > 0”作为 gate，也不使用旧 consensus state 调参。Stage 16 仅生成 `results/stage16/`，不重训 scVI、不重新聚类、不覆盖 Stage 1–13 或 scVI artifact。正式运行命令为：

```bash
$PY 16_rebuild_lam_identity_gate.py --config config/state_modeling.yaml
```

Stage 18 (`18_validate_state15_anchor.py`) 将现有 consensus 的 State 15 精确冻结为 200 个细胞，只做 LAM-core reference anchor 的独立验证，不修改 candidate gate、不重新聚类、不重训 scVI。它使用可用的 777-gene formal LAMCORE，汇总 State 15、全部其他 consensus states、boundary 和 normal/control 的 profile；计算 author-style enrichment、State 15 与 State 18/20/12/7/5 的 marker/program 对照、patient × group pseudobulk、跨患者一致性、已有 `X_scVI` 的邻域以及距离—identity 梯度。结果只写入 `results/stage18/`，其中 `state15_anchor_summary.json` 保存冻结细胞 ID 的 SHA-256 和输入清单，`state15_anchor_report.md` 保存审计摘要。该阶段的结论是诊断性的，不会自动把 State 15 写回 candidate gate 或 atlas。

## Step 7–13 输出

- `results/stage7/`：co-assignment、完整距离层次、consensus assignments 和 stability
- `results/stage8/`：full reference baseline、patient/dataset LOO matching 和 additional loss
- `results/stage9/`：state distance、PAGA/connectivity、split/merge 和 boundary transitions
- `results/stage10/`：逐 state pseudobulk、独立 DE、markers 和 upstream program scores
- `results/stage11/`：patient×state 与连续复现指标
- `results/stage12/`：boundary/normal 辅助邻域结果
- `results/stage13/`：state atlas、hypothesis candidates 和 atlas AnnData
- `results/stage7/state_consensus_with_upstream_annotations.csv`：逐细胞 consensus + upstream annotation
- `results/stage7/state_consensus_state_summary.csv`：按 consensus state 的汇总
- `results/stage16/`：逐细胞 identity evidence、四数据集校准、leave-one-dataset-out 验证、新 candidate assignment 和旧 state 诊断汇总
- `results/stage17/`：GSE190260 upstream positive-reference 漏检的 component、raw-count dropout、competing lineage、反事实标准化和 LODO 审计；不生成新的 candidate assignment
- `results/stage18/`：冻结 State 15 的 777-gene LAMCORE、author-label、marker/program、patient-level pseudobulk、一致性和 scVI latent-neighborhood anchor 验证；不改变任何既有 artifact
- `results/stage19/`：State 15 的 candidate-pool 患者组成基线、author annotation availability、patient-matched profile、7 患者 LOPO、去除 LAM1163 敏感性和跨患者结论；不重新聚类、不重训 scVI、不修改 gate
- `results/stage20/`：以冻结 State 15 为 reference 的 candidate/boundary latent distance、距离分箱、identity/lineage gradient、State 16 LAM/immune 共表达、patient/dataset consistency、boundary 投射和 normal 远端对照；不重训 scVI、不重新聚类、不修改 candidate gate
- `results/stage21/`：将 State 15 仅作为 anchor，在其余 22,061 个细胞上审计 777-gene LAMCORE 与 scVI HVG/旧 gate 的重叠，计算 independent LAMCORE gradient、dataset/patient 一致性、500 次 patient×dataset matched-anchor null、State16/ boundary 梯度和 connectivity；不重训 scVI、不重新聚类、不修改 candidate gate
- `results/stage22/`：在固定 State15-centered k=30 局部图上分解 1–3 hop 方向，自动筛选主要 branches，分析 State16/其他 branch 的 near/mid/far gradient、patient consistency、boundary 投射和每条 branch 的 500 次 matched null；不重训 scVI、不重新聚类、不修改 candidate gate

## Stage 23：State 15 latent-space visualization

Stage 23 只读取现有 `X_scVI`、Stage 20 距离/score、Stage 21 independent scores 和 Stage 22 局部图/branch/null 结果，固定使用 5,378 个 candidate 加 16,883 个 boundary。State 15 仍是冻结的 200-cell anchor；本阶段不重训 scVI、不重新 Leiden/consensus、不修改 candidate gate 或 State 1–20 标签。为避免 Stage 21 排除 anchor 造成可视化缺列，200 个 anchor 的分数在内存中用 Stage 21 的完全相同 score modules 补齐，并优先用 Stage 20 已保存的 program scores；不回写任何输入 artifact。

`23_visualize_state15_latent_space.py` 生成全局 2D UMAP、State 15 局部 kNN 图、distance × independent LAMCORE、State 16 near/mid/far heatmap 和 branch matched-null 图，并生成可按 State、LAMCORE、patient、dataset、candidate/boundary、branch 切换颜色的 3D UMAP/PCA HTML。2D 若已有 `state_model_scvi.h5ad:obsm[X_umap]` 则直接复用；3D UMAP 从既有 20-dimensional `X_scVI` 单进程计算，PCA 作为正交对照。局部图边为同一 `X_scVI`、k=30 邻域的确定性重建，默认最多绘制 100,000 条边以控制交互文件大小。

输出目录：

```text
results/stage23_visualization/
├── 01_global_latent_umap.png/pdf
├── 02_state15_local_knn_states.png/pdf
├── 02_state15_local_knn_lamcore.png/pdf
├── 03_distance_lamcore_gradient.png/pdf
├── 04_state16_program_heatmap.png/pdf
├── 05_branch_matched_null.png/pdf
├── 3d_global_umap.html
├── 3d_state15_local_graph.html
├── 3d_global_pca.html
└── visualization_manifest.json
```

## Stage 22：State 15 局部分支分解

Stage 22 不再判断全局 manifold，而是用既有 `X_scVI` 在 22,261 个 candidate+boundary cells 上建立固定 k=30 局部图。图采用无向化 kNN 的 1–3 hop，State 15 的 200 个细胞始终只作中心 anchor；保留 State 1–20 标签，不重新聚类。branch 自动选择规则为：外部已有 state 与 State15 直接相连的细胞数至少 10，且至少来自 2 个患者；boundary 不被提升为 branch。

本次选出 4 个局部方向：State16（1/2/3-hop 为 96/234/65，10 patients）、State12（22/165/327）、State20（22/320/384）和 State7（11/87/308）。State18 直接连接不足，因此未强行定义为 branch。State16 的 independent LAMCORE 在 near/mid/far 为 0.1820/0.1128/0.1146，7 个患者中 5 个呈 near→far 下降；State16 branch matched-null 经验双侧 p=0.02595。State12、State20、State7 的 matched-null 分别为 0.67465、0.31537、1.0，当前更适合解释为 ordinary lineage adjacency。

Boundary 只在 1–3 hop 范围内按直接 branch-neighbor 优势进行方向投射，10,645 个 boundary 保持 unresolved，不生成新的 LAM 标签。Stage 22 checkpoint 为 `local_branched_lam_manifold_candidate`：State16 是当前最有支持的 LAM-preserving/lineage-transition candidate，其他方向暂不升级为 LAM branch。

## Stage 21：State 15-centered manifold 独立验证

Stage 21 将冻结的 State 15 200 个细胞与验证对象严格分开：State 15 只作 reference anchor，所有主要 gradient 检验在 Stage 20 的 22,061 个非-State15细胞上完成，并直接复用 Stage 20 的 `distance_to_state15`。对 777-gene formal LAMCORE 逐基因标记其是否属于 scVI 4000 HVG、旧 candidate gate marker 及当前表达矩阵，构造 `LAMCORE_full`、`LAMCORE_no_gate`、`LAMCORE_outside_scVI` 和最严格的 `LAMCORE_independent`。

Stage 21 同时输出 anchor-excluded distance bins、dataset/patient adjusted gradient、continuous smooth table、State 16 near/mid/far profile、boundary evidence ranking 和同一 `X_scVI` candidate+boundary scope 的 k=30 connectivity。500 次 null 从非-State15 candidate 中按 State15 的 patient×dataset 组成抽取假 anchor，真实 candidate-only independent slope 与 null 分布比较。结果为：554 个 LAMCORE 基因同时位于 scVI 4000 HVG 外且不属于旧 gate marker；真实 candidate-only independent patient-adjusted slope 为 -0.01547，matched-anchor null empirical two-sided p=0.001996，但 pooled non-State15 rank gradient 的 Spearman rho 为 +0.0757、最近/最远分箱不呈单调下降，患者方向存在异质性。因此 checkpoint 为 `state15_lam_rich_gradient_but_not_robust_manifold`：存在独立 LAM-rich gradient，但暂不足以升级为稳健统一 manifold。

## Stage 20：State 15-centered manifold

Stage 20 固定使用当前 consensus State 15 的 200 个细胞和既有 `state_model_scvi.h5ad` 的 `X_scVI`。主几何 cohort 为 5,378 个 high-confidence candidates 加 16,883 个 inherited boundary cells；normal 只产生远端摘要，不进入 State 15 邻域图。距离、centroid、30-NN State15 neighbor fraction 完全由 latent geometry 计算，State 1–20 只用于描述邻域组成；LAMCORE、program 和 competing-lineage score 在几何轴建立后映射。

该阶段不重训 scVI、不重新 Leiden/consensus clustering、不修改 candidate gate。`results/stage20/stage20_manifest.json` 保存冻结 State 15 cell ID 哈希、latent artifact 和 cohort scope；`state15_centered_manifold.csv` 是后续 State 15-centered 分析的逐细胞主表，`stage20_manifold_report.md` 保存 checkpoint 与解释。

## Stage 24：最终项目材料包

Stage 24 是只读整理阶段。它冻结 Stage 1–23 已有的 state、State15、candidate gate、scVI 和 branch artifact，不重训、不重新聚类、不修改任何既有结果。脚本 `24_finalize_project.py` 只读取已有结果、报告、脚本和声明的 upstream 路径，并将最终写作所需材料写入 `results/stage24_final/`。

它生成：

- `stage_index.csv` / `stage_summary.md`：Stage 1–23 的问题、输入、方法、checkpoint、下一步影响及后续修正；
- `artifact_index.csv`：项目脚本、结果、报告、AnnData/model、配置和声明 upstream 输入的存在性、大小、CSV 行列数、路径与小文件 SHA-256；
- `state_human_cell_analogue.csv/md`：20 个 state 的表达/program/结构证据、支持/冲突证据、LAM 关系和不确定性；
- `other_findings_registry.csv/md`、`narrative_audit.csv`、`terminology_glossary.md`：主线外发现、历史结论变化、定义冲突和文字逻辑审计；
- `final_project_source_materials.md`：完整最终报告原材料包；另有报告摘要、局限、未来方向和 methods/statistics/provenance 附录。

运行命令：

```bash
$PY 24_finalize_project.py
```

Stage 24 最终采用的叙事边界是：State15 为当前最 LAM-rich 的 provisional reference-anchor candidate；Stage20 的 global manifold 结论经 Stage21 限定为不稳健统一 manifold；Stage22 只保留 State15→State16 的局部 transition candidate，State12/20/7 暂为 ordinary lineage adjacency。历史版本仍在对应 Stage 报告和 `narrative_audit.csv` 中保留。
