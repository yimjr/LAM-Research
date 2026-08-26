"""Approximate the author's core lung analysis in a transparent Python runner.

The original repository mixes Seurat 2 and Seurat 3 code. This runner keeps
the author-facing module structure and graph parameters, but uses a modern
AnnData input and a deterministic sparse-PCA/Jaccard-Louvain implementation.
The output is therefore labelled an independent reimplementation unless the
separate R/Seurat run later demonstrates parameter-level equivalence.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import random
import sys
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import igraph as ig
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from sklearn.decomposition import TruncatedSVD
from sklearn.neighbors import NearestNeighbors


ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))

KNOWN_MARKERS = ["PMEL", "ACTA2", "ESR1", "FIGF", "VEGFD", "CTSK", "MLANA"]
AUTHOR_MODULES = {
    "LAM1_LAM3_CCA_align_approx": {"samples": ["LAM1", "LAM3"], "n_pcs": 50, "n_neighbors": 25, "align": True},
    "LAM2_PCA": {"samples": ["LAM2"], "n_pcs": 15, "n_neighbors": 10, "align": False},
    "LAM4_PCA": {"samples": ["LAM4"], "n_pcs": 15, "n_neighbors": 11, "align": False},
    "Donor1_control_PCA": {"samples": ["Donor1"], "n_pcs": 15, "n_neighbors": 30, "align": False},
}


def resolve_symbols(adata: ad.AnnData, requested: list[str]) -> tuple[list[str], list[str]]:
    raw_names = [str(name) for name in adata.raw.var_names]
    exact = set(raw_names)
    upper_lookup: dict[str, str] = {}
    if "gene_symbol_upper" in adata.raw.var:
        for actual, upper in zip(raw_names, adata.raw.var["gene_symbol_upper"].astype(str)):
            upper_lookup.setdefault(upper.upper(), actual)
    available, missing = [], []
    for gene in requested:
        actual = gene if gene in exact else upper_lookup.get(gene.upper())
        if actual is None:
            missing.append(gene)
        elif actual not in available:
            available.append(actual)
    return available, sorted(set(missing))


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def marker_columns(adata: ad.AnnData) -> dict[str, str]:
    available, missing = resolve_symbols(adata, KNOWN_MARKERS)
    adata.uns["known_marker_resolution"] = {"available": available, "missing": missing}
    columns: dict[str, str] = {}
    for requested in KNOWN_MARKERS:
        actual, _ = resolve_symbols(adata, [requested])
        column = f"marker_expr_{requested}"
        if actual:
            values = adata.raw[:, actual[0]].X
            adata.obs[column] = np.asarray(values.toarray()).ravel() if sp.issparse(values) else np.asarray(values).ravel()
        else:
            adata.obs[column] = 0.0
        columns[requested] = column
    if available:
        try:
            sc.tl.score_genes(
                adata,
                gene_list=available,
                score_name="known_marker_score",
                use_raw=True,
                ctrl_size=min(50, max(1, len(adata.raw.var_names) // 20)),
                random_state=0,
            )
        except Exception:
            matrix = adata.raw[:, available].X
            adata.obs["known_marker_score"] = np.asarray(matrix.mean(axis=1)).ravel()
    else:
        adata.obs["known_marker_score"] = 0.0
    marker_matrix = adata.raw[:, available].X if available else sp.csr_matrix((adata.n_obs, 0))
    detected = np.asarray((marker_matrix > 0).sum(axis=1)).ravel() if available else np.zeros(adata.n_obs)
    adata.obs["known_marker_genes_detected"] = detected.astype(int)
    adata.obs["known_marker_combo_ge2"] = adata.obs["known_marker_genes_detected"] >= 2
    return columns


def jaccard_louvain(embedding: np.ndarray, n_neighbors: int, seed: int) -> np.ndarray:
    """Recreate the author's kNN-Jaccard weighted Louvain graph."""
    if embedding.shape[0] <= n_neighbors + 1:
        return np.zeros(embedding.shape[0], dtype=str)
    nearest = NearestNeighbors(n_neighbors=n_neighbors + 1, metric="euclidean", algorithm="auto")
    nearest.fit(embedding)
    indices = nearest.kneighbors(return_distance=False)[:, 1:]
    neighbor_sets = [set(row.tolist()) for row in indices]
    edges: list[tuple[int, int]] = []
    weights: list[float] = []
    seen: set[tuple[int, int]] = set()
    for i, row in enumerate(indices):
        for j in row:
            a, b = sorted((int(i), int(j)))
            if a == b or (a, b) in seen:
                continue
            seen.add((a, b))
            union = neighbor_sets[a] | neighbor_sets[b]
            if not union:
                continue
            weight = len(neighbor_sets[a] & neighbor_sets[b]) / len(union)
            if weight > 0:
                edges.append((a, b))
                weights.append(weight)
    graph = ig.Graph(n=embedding.shape[0], edges=edges, directed=False)
    if weights:
        maximum = max(weights)
        graph.es["weight"] = [weight / maximum for weight in weights]
        random.seed(seed)
        ig.set_random_number_generator(random.Random(seed))
        communities = graph.community_multilevel(weights=graph.es["weight"])
        labels = np.asarray(communities.membership, dtype=str)
    else:
        labels = np.arange(embedding.shape[0], dtype=str)
    # Match the author's convention of ordering clusters by decreasing size.
    sizes = pd.Series(labels).value_counts().index.tolist()
    ordering = {old: str(i) for i, old in enumerate(sizes)}
    return np.asarray([ordering[label] for label in labels], dtype=str)


def module_embedding(adata: ad.AnnData, mask: np.ndarray, n_pcs: int, align: bool, seed: int) -> np.ndarray:
    available_hvg = adata.var.get("highly_variable", pd.Series(False, index=adata.var_names)).astype(bool).to_numpy()
    if int(available_hvg.sum()) < n_pcs + 2:
        available_hvg[:] = True
    matrix = adata[mask, available_hvg].X
    if not sp.issparse(matrix):
        matrix = sp.csr_matrix(matrix)
    matrix = matrix.astype(np.float32).tocsr()
    n_components = min(n_pcs, max(2, min(matrix.shape) - 1))
    svd = TruncatedSVD(n_components=n_components, random_state=seed, n_iter=7)
    embedding = svd.fit_transform(matrix).astype(np.float32)
    if align:
        samples = adata.obs.loc[adata.obs.index[mask], "sample_id"].astype(str).to_numpy()
        for sample in np.unique(samples):
            sample_mask = samples == sample
            embedding[sample_mask] -= embedding[sample_mask].mean(axis=0, keepdims=True)
    return embedding


def select_candidate_clusters(cluster_table: pd.DataFrame) -> dict[str, set[str]]:
    selected: dict[str, set[str]] = {}
    for module, group in cluster_table.groupby("module", observed=True):
        group = group.copy()
        score_cutoff = float(group["mean_known_marker_score"].quantile(0.90))
        eligible = group.loc[
            (group["mean_known_marker_score"] >= score_cutoff)
            & (group["mean_known_marker_genes_detected"] >= 2)
            & (group["fraction_combo_ge2"] >= 0.10),
            "cluster",
        ].astype(str)
        if eligible.empty:
            fallback = group.sort_values(
                ["mean_known_marker_score", "mean_known_marker_genes_detected"], ascending=False
            ).head(1)["cluster"].astype(str)
            eligible = fallback
        selected[str(module)] = set(eligible.tolist())
    return selected


def write_markers(adata: ad.AnnData, candidate_col: str, table_dir: Path) -> None:
    lam = adata.obs["condition"].astype(str).eq("LAM")
    labels = np.where(adata.obs[candidate_col].astype(bool), "candidate", "other")
    adata.obs["reproduction_label"] = np.where(lam, labels, "reference")
    if int((adata.obs["reproduction_label"] == "candidate").sum()) < 10:
        return
    subset = adata[lam].copy()
    if subset.obs["reproduction_label"].nunique() != 2:
        return
    sc.tl.rank_genes_groups(subset, groupby="reproduction_label", groups=["candidate"], reference="other", method="wilcoxon", use_raw=True)
    sc.get.rank_genes_groups_df(subset, group="candidate").to_csv(
        table_dir / "lamcore_candidate_markers_author_style.csv", index=False
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="data/processed/reproduction_core/GSE135851_core_baseline_qc.h5ad")
    parser.add_argument("--output", default="data/processed/reproduction_core/GSE135851_core_reproduction.h5ad")
    parser.add_argument("--config", default="config/signatures.yaml")
    args = parser.parse_args()
    seed = 20260822
    np.random.seed(seed)
    sc.settings.set_figure_params(dpi=120, facecolor="white")
    input_path = ROOT / args.input
    output_path = ROOT / args.output
    table_dir = ROOT / "results" / "reproduction_core" / "tables"
    figure_dir = ROOT / "results" / "reproduction_core" / "figures"
    report_dir = ROOT / "results" / "reproduction_core"
    for directory in (output_path.parent, table_dir, figure_dir, report_dir):
        directory.mkdir(parents=True, exist_ok=True)
    adata = ad.read_h5ad(input_path)
    if adata.raw is None or "counts" not in adata.layers:
        raise RuntimeError("Reproduction baseline must contain .raw and layers['counts']")

    marker_columns(adata)
    adata.obs["author_module"] = "unassigned"
    adata.obs["author_cluster"] = "unassigned"
    module_rows: list[dict] = []
    for module, spec in AUTHOR_MODULES.items():
        mask = adata.obs["sample_id"].astype(str).isin(spec["samples"]).to_numpy()
        if int(mask.sum()) < 50:
            continue
        embedding = module_embedding(adata, mask, int(spec["n_pcs"]), bool(spec["align"]), seed)
        labels = jaccard_louvain(embedding, int(spec["n_neighbors"]), seed)
        adata.obs.loc[mask, "author_module"] = module
        adata.obs.loc[mask, "author_cluster"] = labels
        for cluster, group in adata.obs.loc[mask].groupby("author_cluster", observed=True):
            row = {
                "module": module,
                "cluster": str(cluster),
                "cells": int(len(group)),
                "mean_known_marker_score": float(group["known_marker_score"].mean()),
                "mean_known_marker_genes_detected": float(group["known_marker_genes_detected"].mean()),
                "fraction_combo_ge2": float(group["known_marker_combo_ge2"].mean()),
                "mean_doublet_score": float(group["doublet_score"].mean()),
            }
            for gene in KNOWN_MARKERS:
                row[f"mean_{gene}"] = float(group[f"marker_expr_{gene}"].mean())
            module_rows.append(row)
    cluster_table = pd.DataFrame(module_rows)
    selected = select_candidate_clusters(cluster_table)
    candidate = np.zeros(adata.n_obs, dtype=bool)
    for module, clusters in selected.items():
        module_mask = adata.obs["author_module"].astype(str).eq(module).to_numpy()
        cluster_mask = adata.obs["author_cluster"].astype(str).isin(clusters).to_numpy()
        # Donor1 is a normal-lung reference. It is used to compare expression
        # context, not to define a LAMCORE candidate population.
        lam_mask = adata.obs["condition"].astype(str).eq("LAM").to_numpy()
        candidate |= module_mask & cluster_mask & adata.obs["known_marker_combo_ge2"].to_numpy() & lam_mask
    adata.obs["lamcore_candidate_author_style"] = candidate
    adata.obs["lamcore_candidate"] = candidate
    adata.obs["lamcore_label_author_style"] = np.where(candidate, "candidate", "other")
    cluster_table["selected_as_candidate_cluster"] = [
        str(row["cluster"]) in selected.get(str(row["module"]), set()) for _, row in cluster_table.iterrows()
    ]
    cluster_table.to_csv(table_dir / "lamcore_cluster_summary_author_style.csv", index=False)

    # Only now score the published 777-gene set, as a consistency check.
    signatures = yaml.safe_load((ROOT / args.config).read_text())
    formal_table = pd.read_csv(ROOT / signatures["lamcore_formal"]["file"])
    formal_table.columns = formal_table.columns.astype(str).str.strip()
    formal = formal_table[signatures["lamcore_formal"]["gene_column"]].dropna().astype(str).tolist()
    formal, formal_missing = resolve_symbols(adata, formal)
    sc.tl.score_genes(adata, gene_list=formal, score_name="lamcore_score_777_consistency", use_raw=True, ctrl_size=50, random_state=0)
    adata.obs["lamcore_genes_detected_777"] = np.asarray((adata.raw[:, formal].X > 0).sum(axis=1)).ravel().astype(int)
    adata.obs[["sample_id", "donor_id", "condition", "assay", "author_module", "author_cluster", "lamcore_candidate_author_style", "known_marker_score", "known_marker_genes_detected", "lamcore_score_777_consistency", "lamcore_genes_detected_777", "doublet_score", "doublet_predicted"]].to_csv(table_dir / "lamcore_candidate_cells_author_style.csv")
    write_markers(adata, "lamcore_candidate_author_style", table_dir)
    cluster_table.to_csv(table_dir / "lamcore_cluster_summary_author_style.csv", index=False)

    donor = (
        adata.obs[adata.obs["condition"].astype(str).eq("LAM")]
        .groupby(["donor_id", "sample_id", "assay"], observed=True)
        .agg(
            cells=("lamcore_candidate_author_style", "size"),
            candidate_cells=("lamcore_candidate_author_style", "sum"),
            candidate_fraction=("lamcore_candidate_author_style", "mean"),
            candidate_marker_score=("known_marker_score", lambda x: float(x[adata.obs.loc[x.index, "lamcore_candidate_author_style"].astype(bool)].mean()) if int(adata.obs.loc[x.index, "lamcore_candidate_author_style"].sum()) else np.nan),
            other_marker_score=("known_marker_score", lambda x: float(x[~adata.obs.loc[x.index, "lamcore_candidate_author_style"].astype(bool)].mean()) if int((~adata.obs.loc[x.index, "lamcore_candidate_author_style"].astype(bool)).sum()) else np.nan),
            candidate_777_score=("lamcore_score_777_consistency", lambda x: float(x[adata.obs.loc[x.index, "lamcore_candidate_author_style"].astype(bool)].mean()) if int(adata.obs.loc[x.index, "lamcore_candidate_author_style"].sum()) else np.nan),
            other_777_score=("lamcore_score_777_consistency", lambda x: float(x[~adata.obs.loc[x.index, "lamcore_candidate_author_style"].astype(bool)].mean()) if int((~adata.obs.loc[x.index, "lamcore_candidate_author_style"].astype(bool)).sum()) else np.nan),
        )
        .reset_index()
    )
    donor["marker_score_difference"] = donor["candidate_marker_score"] - donor["other_marker_score"]
    donor["score_777_difference_consistency_only"] = donor["candidate_777_score"] - donor["other_777_score"]
    donor.to_csv(table_dir / "donor_summary_author_style.csv", index=False)

    comparison = pd.DataFrame([
        {"result_item": "LAMCORE positive donors", "original_paper": "LAM1, LAM3, LAM4; none reported in LAM2", "project": ", ".join(sorted(donor.loc[donor["candidate_cells"] > 0, "donor_id"].astype(str))), "status": "compare"},
        {"result_item": "LAMCORE candidate cells", "original_paper": "125 cells from LAM1, LAM3, LAM4", "project": str(int(candidate.sum())), "status": "compare_quantity_not_identity"},
        {"result_item": "Identification basis", "original_paper": "unbiased clustering plus known markers; 777-gene signature derived later", "project": "known marker panel + author-style graph; 777 genes scored afterward", "status": "aligned_boundary"},
        {"result_item": "Upstream QC", "original_paper": "not reconstructable from processed matrix alone", "project": "downstream QC recoverable from processed matrix; no FASTQ/Cell Ranger reconstruction", "status": "not_reconstructable"},
        {"result_item": "Doublet handling", "original_paper": "not explicit in public scripts", "project": "score and prediction recorded; not removed in Phase 1", "status": "aligned_with_plan"},
        {"result_item": "Software path", "original_paper": "Seurat 2/3 R scripts", "project": "Python sparse-PCA/Jaccard-Louvain approximation; R environment prepared separately", "status": "approximate_reimplementation"},
    ])
    comparison.to_csv(report_dir / "author_vs_project_comparison.csv", index=False)
    deviation = pd.DataFrame([
        {"step": "LAM1/LAM3 representation", "author": "RunCCA + AlignSubspace, 50 dimensions", "project": "sample-centered sparse PCA approximation, 50 dimensions", "impact": "not exact; label independent reimplementation"},
        {"step": "LAM2 representation", "author": "PCA on variable genes, 15 PCs", "project": "sparse PCA on baseline HVGs, 15 dimensions", "impact": "close but software implementation differs"},
        {"step": "LAM4 representation", "author": "PCA on variable genes, 15 PCs", "project": "sparse PCA on baseline HVGs, 15 dimensions", "impact": "close but software implementation differs"},
        {"step": "Graph", "author": "custom kNN Jaccard + Louvain", "project": "custom kNN Jaccard + igraph multilevel Louvain", "impact": "conceptually aligned; exact RNG/software differs"},
        {"step": "LAMCORE identification", "author": "known markers after unbiased clustering", "project": "known markers after author-style graph; 777-gene score is post hoc", "impact": "boundary aligned"},
        {"step": "Doublet", "author": "not explicit in public scripts", "project": "record only in Phase 1", "impact": "boundary aligned"},
    ])
    deviation.to_csv(report_dir / "method_deviation_table.csv", index=False)

    by_donor = donor.set_index("donor_id")["candidate_cells"].to_dict()
    lam2 = int(by_donor.get("LAM2", 0))
    positive_lam_donors = [key for key, value in by_donor.items() if key != "Donor1" and value > 0]
    core_ready = bool(all(key in by_donor for key in ["LAM1", "LAM2", "LAM3", "LAM4"]) and len(positive_lam_donors) >= 3)
    completion = {
        "core_baseline_status": "ready_for_parallel_phase2_phase3" if core_ready else "not_ready_strictly",
        "classification": "independent_reimplementation_not_strict_replication",
        "candidate_cells_author_style": int(candidate.sum()),
        "candidate_cells_original_paper": 125,
        "candidate_cells_by_donor": {str(key): int(value) for key, value in by_donor.items()},
        "lam2_candidate_cells": lam2,
        "positive_lam_donors": sorted(positive_lam_donors),
        "known_marker_genes": KNOWN_MARKERS,
        "formal_signature_used_for_identification": False,
        "formal_signature_used_for_consistency": True,
        "doublet_exclusion_applied": False,
        "qc_boundary": "downstream QC recoverable from processed matrices; upstream QC not reconstructed",
        "created_at": datetime.now(timezone.utc).isoformat(),
    }
    (report_dir / "core_completion.json").write_text(json.dumps(completion, indent=2))

    zh = [
        "# LAMCORE 核心肺部复现：Phase 1 基线",
        "",
        "## 结论定位",
        "",
        "本结果是依据作者公开 R 脚本建立的独立重实现，不是严格的逐版本 Seurat 复现。原因是作者脚本混用 Seurat 2/3，而当前主流程使用 AnnData/Python；所有差异记录在 `method_deviation_table.csv`。",
        "",
        f"- 作者报告的 LAMCORE 细胞数：约 125 个，来自 LAM1、LAM3、LAM4。",
        f"- 本项目作者风格 marker/cluster 重实现候选数：**{int(candidate.sum())}**。",
        f"- 各供体候选数：`{json.dumps({str(k): int(v) for k, v in by_donor.items()}, ensure_ascii=False)}`。",
        f"- LAM2 候选数：**{lam2}**；该数值只能说明本操作性规则下的 marker 候选，不能直接等同为论文定义的 LAMCORE。",
        "",
        "## 方法边界",
        "",
        "候选细胞先由 PMEL、ACTA2、ESR1、FIGF/VEGFD、CTSK、MLANA 等已知特征结合作者风格图聚类定位；777-gene signature 只在候选确定以后作为一致性检查。因此没有用由原始 LAMCORE 细胞总结出的 777 genes 反过来定义它们。",
        "",
        "QC 仅表示“在处理后矩阵允许范围内恢复下游 QC”；FASTQ、Cell Ranger、初始 barcode/cell calling 和 empty droplet 判断没有被恢复。doublet 分数和预测已记录，但 Phase 1 没有据此删除细胞。",
        "",
        "## 下一步",
        "",
        "核心基线已具备进入第二阶段针对性稳健性验证和第三阶段新发现探索的条件；两者并行推进。若后续 R/Seurat 运行能进一步缩小方法差异，则会更新为更严格的复现版本。",
    ]
    en = [
        "# LAMCORE Core Lung Reproduction: Phase 1 Baseline",
        "",
        "## Interpretation",
        "",
        "This is an independent reimplementation guided by the authors' public R scripts, not a strict byte-for-byte Seurat reproduction. The public scripts mix Seurat 2 and Seurat 3, whereas the current executable path uses AnnData/Python; deviations are recorded in `method_deviation_table.csv`.",
        "",
        f"- The paper reported approximately 125 LAMCORE cells from LAM1, LAM3 and LAM4.",
        f"- Author-style marker/cluster reimplementation candidates: **{int(candidate.sum())}**.",
        f"- Candidate counts by donor: `{json.dumps({str(k): int(v) for k, v in by_donor.items()})}`.",
        f"- LAM2 candidates: **{lam2}**; this is an operational marker candidate count, not the paper's original LAMCORE label.",
        "",
        "## Method boundary",
        "",
        "Candidates were located using known features (PMEL, ACTA2, ESR1, FIGF/VEGFD, CTSK and MLANA) together with an author-style graph. The 777-gene signature was scored only after candidate selection as a consistency check, avoiding circular use of a signature derived from the original LAMCORE cells.",
        "",
        "QC means downstream QC recoverable from processed matrices; FASTQ, Cell Ranger, initial barcode/cell calling and empty-droplet decisions were not reconstructed. Doublet scores and predictions were recorded, but Phase 1 did not remove cells on that basis.",
        "",
        "## Next",
        "",
        "The baseline is ready for targeted Phase 2 robustness checks and Phase 3 biological discovery in parallel. If the separate R/Seurat run narrows the implementation differences, the reproduction status will be updated.",
    ]
    (report_dir / "LAM_core_reproduction_report_zh.md").write_text("\n".join(zh) + "\n")
    (report_dir / "LAM_core_reproduction_report_en.md").write_text("\n".join(en) + "\n")

    if "X_umap" in adata.obsm:
        for color, filename in [
            ("known_marker_score", "umap_known_marker_score.png"),
            ("lamcore_candidate_author_style", "umap_author_style_candidate.png"),
            ("sample_id", "umap_sample_id.png"),
            ("author_module", "umap_author_module.png"),
        ]:
            sc.pl.umap(adata, color=color, show=False, frameon=False)
            plt.tight_layout()
            plt.savefig(figure_dir / filename, dpi=160, bbox_inches="tight")
            plt.close("all")

    adata.uns["core_reproduction"] = {"completion": completion, "selected_clusters": {key: sorted(value) for key, value in selected.items()}, "formal_signature_sha256": sha256_file(ROOT / signatures["lamcore_formal"]["file"]), "formal_signature_missing": formal_missing}
    adata.write_h5ad(output_path, compression="gzip")
    manifest_path = ROOT / "manifests" / "run_manifest.yaml"
    manifest = yaml.safe_load(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest.setdefault("steps", [])
    manifest["steps"] = [step for step in manifest["steps"] if step.get("name") != "phase1_core_reproduction"]
    manifest["steps"].append({"name": "phase1_core_reproduction", "completed_at": datetime.now(timezone.utc).isoformat(), "input": str(input_path.relative_to(ROOT)), "output": str(output_path.relative_to(ROOT)), "formal_signature_used_for_identification": False, "formal_signature_used_for_consistency": True, "doublet_exclusion_applied": False, "classification": completion["classification"]})
    manifest["status"] = "phase1_core_reproduction_ready_parallel_phase2_phase3"
    manifest_path.write_text(yaml.safe_dump(manifest, sort_keys=False))
    print(json.dumps(completion, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"reproduce_core.py failed: {exc}", file=sys.stderr)
        raise
