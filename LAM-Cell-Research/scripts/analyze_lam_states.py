"""Score the formal LAMCORE signature and exploratory cell-state programs."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".cache" / "matplotlib"))

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import yaml
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_gene_names(adata: ad.AnnData, genes: list[str]) -> tuple[list[str], list[str]]:
    """Resolve requested symbols against the raw matrix, case-insensitively."""
    raw_names = [str(name) for name in adata.raw.var_names]
    exact = set(raw_names)
    upper_lookup: dict[str, str] = {}
    if "gene_symbol_upper" in adata.raw.var:
        for actual, upper in zip(raw_names, adata.raw.var["gene_symbol_upper"].astype(str)):
            upper_lookup.setdefault(upper.upper(), actual)
    available: list[str] = []
    missing: list[str] = []
    for gene in genes:
        requested = str(gene)
        actual = requested if requested in exact else upper_lookup.get(requested.upper())
        if actual is None:
            missing.append(requested)
        elif actual not in available:
            available.append(actual)
    return available, sorted(set(missing))


def score_gene_set(adata: ad.AnnData, genes: list[str], score_name: str) -> dict:
    available, missing = resolve_gene_names(adata, genes)
    if len(available) < 2:
        adata.obs[score_name] = np.nan
        return {
            "requested_genes": len(genes),
            "available": available,
            "missing": missing,
            "status": "insufficient_genes",
        }
    try:
        ctrl_size = min(50, max(1, len(adata.raw.var_names) // 20))
        sc.tl.score_genes(
            adata,
            gene_list=available,
            score_name=score_name,
            use_raw=True,
            ctrl_size=ctrl_size,
            random_state=0,
        )
        status = "completed"
    except Exception as exc:
        adata.obs[score_name] = np.nan
        status = f"failed: {type(exc).__name__}: {exc}"
    return {
        "requested_genes": len(genes),
        "available": available,
        "missing": missing,
        "status": status,
    }


def detected_gene_count(adata: ad.AnnData, genes: list[str]) -> tuple[pd.Series, dict]:
    available, missing = resolve_gene_names(adata, genes)
    if not available:
        return pd.Series(0, index=adata.obs_names, dtype="int64"), {
            "available": [],
            "missing": missing,
            "status": "no_genes_available",
        }
    matrix = adata.raw[:, available].X
    detected = (matrix > 0).sum(axis=1)
    values = np.asarray(detected).ravel().astype(int)
    return pd.Series(values, index=adata.obs_names), {
        "available": available,
        "missing": missing,
        "status": "completed",
    }


def load_formal_signature(spec: dict) -> tuple[list[str], dict]:
    path = ROOT / spec["file"]
    if not path.exists():
        raise FileNotFoundError(f"Formal LAMCORE signature file is missing: {path}")
    table = pd.read_csv(path)
    table.columns = table.columns.astype(str).str.strip()
    column = spec.get("gene_column", "Gene")
    if column not in table.columns:
        raise ValueError(f"Signature file does not contain the configured gene column: {column}")
    genes = table[column].dropna().astype(str).str.strip()
    genes = list(dict.fromkeys(gene for gene in genes if gene))
    expected = int(spec.get("expected_genes", len(genes)))
    if len(genes) != expected:
        raise ValueError(f"Expected {expected} formal signature genes, found {len(genes)} in {path}")
    metadata = {
        "source": spec.get("source"),
        "source_url": spec.get("source_url"),
        "file": str(path.relative_to(ROOT)),
        "sha256": sha256_file(path),
        "rows": int(len(table)),
        "unique_genes": int(len(genes)),
        "gene_column": column,
    }
    return genes, metadata


def sparse_sum(matrix, axis=0) -> np.ndarray:
    values = matrix.sum(axis=axis)
    return np.asarray(values).ravel()


def write_umap(adata: ad.AnnData, color: str, output: Path) -> None:
    if color not in adata.obs:
        return
    sc.pl.umap(adata, color=color, show=False, frameon=False)
    plt.tight_layout()
    plt.savefig(output, dpi=160, bbox_inches="tight")
    plt.close("all")


def discover_candidates(
    adata: ad.AnnData,
    score_col: str,
    detected_col: str,
    cluster_col: str,
    cluster_quantile: float,
    cell_quantile: float,
    min_detected: int,
) -> tuple[pd.Series, pd.DataFrame, str, float]:
    cluster_summary = (
        adata.obs.groupby(cluster_col, observed=True)
        .agg(
            mean_score=(score_col, "mean"),
            mean_genes_detected=(detected_col, "mean"),
            cells=(score_col, "size"),
        )
        .reset_index()
    )
    cluster_cutoff = float(cluster_summary["mean_score"].quantile(cluster_quantile))
    eligible = cluster_summary.loc[
        (cluster_summary["mean_score"] >= cluster_cutoff)
        & (cluster_summary["mean_genes_detected"] >= min_detected),
        cluster_col,
    ].astype(str)
    candidate = adata.obs[cluster_col].astype(str).isin(set(eligible))
    method = "cluster_mean_score"
    if int(candidate.sum()) < 10:
        score_cutoff = float(adata.obs[score_col].quantile(cell_quantile))
        candidate = (adata.obs[score_col] >= score_cutoff) & (adata.obs[detected_col] >= min_detected)
        method = "cell_score_fallback"
    else:
        score_cutoff = cluster_cutoff
    return candidate.astype(bool), cluster_summary, method, score_cutoff


def summarize_candidate(
    adata: ad.AnnData,
    candidate_col: str,
    program_columns: list[str],
    min_cells: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    lam_obs = adata.obs[adata.obs["condition"].astype(str) == "LAM"]
    donor_rows: list[dict] = []
    for donor_id, donor in lam_obs.groupby("donor_id", observed=True):
        candidate = donor[donor[candidate_col].astype(bool)]
        other = donor[~donor[candidate_col].astype(bool)]
        base = {
            "donor_id": donor_id,
            "candidate_cells": int(len(candidate)),
            "other_cells": int(len(other)),
        }
        for column in program_columns:
            usable = len(candidate) >= min_cells and len(other) >= min_cells
            base[f"{column}_candidate_mean"] = float(candidate[column].mean()) if usable else np.nan
            base[f"{column}_other_mean"] = float(other[column].mean()) if usable else np.nan
            base[f"{column}_difference"] = (
                base[f"{column}_candidate_mean"] - base[f"{column}_other_mean"] if usable else np.nan
            )
            base[f"{column}_candidate_higher"] = (
                base[f"{column}_candidate_mean"] > base[f"{column}_other_mean"] if usable else np.nan
            )
        donor_rows.append(base)
    donor_summary = pd.DataFrame(donor_rows)
    concordance_rows: list[dict] = []
    lam_donor_count = int(lam_obs["donor_id"].nunique())
    for column in program_columns:
        difference_col = f"{column}_difference"
        higher_col = f"{column}_candidate_higher"
        valid = donor_summary[difference_col].notna() if not donor_summary.empty else pd.Series(dtype=bool)
        valid_rows = donor_summary.loc[valid] if not donor_summary.empty else donor_summary
        positive = int(valid_rows[higher_col].sum()) if not valid_rows.empty else 0
        concordance_rows.append(
            {
                "program": column,
                "n_lam_donors_tested": int(len(valid_rows)),
                "n_candidate_higher": positive,
                "concordance_fraction": positive / len(valid_rows) if len(valid_rows) else np.nan,
                "passes_observed_donor_rule": positive >= 3 and len(valid_rows) >= 3,
                "passes_3_of_4_rule": positive >= 3 and len(valid_rows) == lam_donor_count,
            }
        )
    return donor_summary, pd.DataFrame(concordance_rows)


def write_pseudobulk(adata: ad.AnnData, candidate_col: str, output: Path) -> None:
    counts = adata.layers["counts"]
    if not sp.issparse(counts):
        counts = sp.csr_matrix(counts)
    rows: list[np.ndarray] = []
    index: list[str] = []
    for donor_id, donor_indices in adata.obs.groupby("donor_id", observed=True).groups.items():
        donor_mask = adata.obs.index.isin(donor_indices)
        for label, label_mask in (("candidate", adata.obs[candidate_col].to_numpy()), ("other", ~adata.obs[candidate_col].to_numpy())):
            cell_mask = donor_mask & label_mask
            if not cell_mask.any():
                continue
            rows.append(sparse_sum(counts[cell_mask], axis=0))
            index.append(f"{donor_id}__{label}")
    if rows:
        pd.DataFrame(rows, index=index, columns=adata.var_names).to_csv(output)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/analysis.yaml")
    parser.add_argument("--signatures", default="config/signatures.yaml")
    args = parser.parse_args()
    config = yaml.safe_load((ROOT / args.config).read_text())
    signatures = yaml.safe_load((ROOT / args.signatures).read_text())
    seed = int(config["random_seed"])
    np.random.seed(seed)
    sc.settings.set_figure_params(dpi=120, facecolor="white")

    accession = config["accession"]
    input_path = ROOT / config["paths"]["processed"] / f"{accession}_preprocessed.h5ad"
    output_path = ROOT / config["paths"]["processed"] / f"{accession}_lam_states.h5ad"
    figure_dir = ROOT / config["paths"]["results"] / "figures"
    table_dir = ROOT / config["paths"]["results"] / "tables"
    report_dir = ROOT / config["paths"]["results"] / "report"
    for directory in (figure_dir, table_dir, report_dir):
        directory.mkdir(parents=True, exist_ok=True)
    if not input_path.exists():
        raise FileNotFoundError(f"Missing preprocessed AnnData: {input_path}")

    adata = ad.read_h5ad(input_path)
    if adata.raw is None:
        raise RuntimeError("AnnData is missing .raw; the preprocessing contract was not satisfied")

    formal_genes, formal_file_metadata = load_formal_signature(signatures["lamcore_formal"])
    fallback_genes = signatures["lamcore_fallback"]["genes"]
    signature_status: dict[str, dict] = {
        "lamcore_formal_file": formal_file_metadata,
        "lamcore_formal": score_gene_set(adata, formal_genes, "lamcore_score_formal"),
        "lamcore_fallback": score_gene_set(adata, fallback_genes, "lamcore_score_fallback"),
    }
    formal_detected, formal_detect_status = detected_gene_count(adata, formal_genes)
    fallback_detected, fallback_detect_status = detected_gene_count(adata, fallback_genes)
    signature_status["lamcore_formal_detection"] = formal_detect_status
    signature_status["lamcore_fallback_detection"] = fallback_detect_status
    adata.obs["lamcore_genes_detected_formal"] = formal_detected
    adata.obs["lamcore_genes_detected_fallback"] = fallback_detected

    for name, genes in signatures["programs"].items():
        column = f"program_{name}"
        signature_status[name] = score_gene_set(adata, genes, column)

    primary_cluster = "leiden_r0_8"
    if primary_cluster not in adata.obs:
        cluster_keys = [key for key in adata.obs.columns if key.startswith("leiden_")]
        if not cluster_keys:
            raise RuntimeError("No Leiden cluster column is available for LAMCORE candidate discovery")
        primary_cluster = cluster_keys[0]
    analysis_cfg = config["analysis"]
    formal_candidate, formal_cluster_summary, formal_method, formal_cutoff = discover_candidates(
        adata,
        "lamcore_score_formal",
        "lamcore_genes_detected_formal",
        primary_cluster,
        float(analysis_cfg["lamcore_cluster_quantile"]),
        float(analysis_cfg["lamcore_cell_quantile_fallback"]),
        int(analysis_cfg["min_formal_signature_genes_detected"]),
    )
    fallback_candidate, fallback_cluster_summary, fallback_method, fallback_cutoff = discover_candidates(
        adata,
        "lamcore_score_fallback",
        "lamcore_genes_detected_fallback",
        primary_cluster,
        float(analysis_cfg["lamcore_cluster_quantile"]),
        float(analysis_cfg["lamcore_cell_quantile_fallback"]),
        int(analysis_cfg["min_fallback_marker_genes_detected"]),
    )
    adata.obs["lamcore_candidate_formal"] = formal_candidate
    adata.obs["lamcore_candidate_fallback"] = fallback_candidate
    adata.obs["lamcore_candidate"] = formal_candidate
    adata.obs["lamcore_score"] = adata.obs["lamcore_score_formal"]
    adata.obs["lamcore_marker_detected"] = adata.obs["lamcore_genes_detected_formal"]
    adata.obs["lamcore_label"] = np.where(formal_candidate, "candidate", "other")
    adata.obs["lamcore_label_fallback"] = np.where(fallback_candidate, "candidate", "other")

    formal_cluster_summary.to_csv(table_dir / "lamcore_cluster_summary_formal.csv", index=False)
    fallback_cluster_summary.to_csv(table_dir / "lamcore_cluster_summary_fallback.csv", index=False)
    score_columns = [
        "sample_id",
        "donor_id",
        "condition",
        "assay",
        "lamcore_score_formal",
        "lamcore_genes_detected_formal",
        "lamcore_candidate_formal",
        "lamcore_score_fallback",
        "lamcore_genes_detected_fallback",
        "lamcore_candidate_fallback",
    ]
    adata.obs[score_columns].to_csv(table_dir / "lamcore_cell_scores.csv")

    formal_candidate_count = int(formal_candidate.sum())
    if formal_candidate_count >= 10 and adata.obs["lamcore_label"].nunique() == 2:
        sc.tl.rank_genes_groups(
            adata,
            groupby="lamcore_label",
            groups=["candidate"],
            reference="other",
            method="wilcoxon",
            use_raw=True,
            key_added="lamcore_markers_formal",
        )
        markers = sc.get.rank_genes_groups_df(adata, group="candidate", key="lamcore_markers_formal")
        markers.to_csv(table_dir / "lamcore_candidate_markers_formal_exploratory.csv", index=False)
    fallback_candidate_count = int(fallback_candidate.sum())
    if fallback_candidate_count >= 10 and adata.obs["lamcore_label_fallback"].nunique() == 2:
        sc.tl.rank_genes_groups(
            adata,
            groupby="lamcore_label_fallback",
            groups=["candidate"],
            reference="other",
            method="wilcoxon",
            use_raw=True,
            key_added="lamcore_markers_fallback",
        )
        markers = sc.get.rank_genes_groups_df(adata, group="candidate", key="lamcore_markers_fallback")
        markers.to_csv(table_dir / "lamcore_candidate_markers_fallback_exploratory.csv", index=False)

    program_columns = [f"program_{name}" for name in signatures["programs"]]
    formal_donor_summary, formal_concordance = summarize_candidate(
        adata, "lamcore_candidate_formal", program_columns, int(analysis_cfg["min_cells_per_donor_for_summary"])
    )
    fallback_donor_summary, fallback_concordance = summarize_candidate(
        adata, "lamcore_candidate_fallback", program_columns, int(analysis_cfg["min_cells_per_donor_for_summary"])
    )
    formal_donor_summary.to_csv(table_dir / "donor_state_summary.csv", index=False)
    formal_donor_summary.to_csv(table_dir / "donor_state_summary_formal.csv", index=False)
    fallback_donor_summary.to_csv(table_dir / "donor_state_summary_fallback.csv", index=False)
    formal_concordance.to_csv(table_dir / "program_concordance.csv", index=False)
    formal_concordance.to_csv(table_dir / "program_concordance_formal.csv", index=False)
    fallback_concordance.to_csv(table_dir / "program_concordance_fallback.csv", index=False)

    lam_obs = adata.obs[adata.obs["condition"].astype(str) == "LAM"]
    lam_donors = sorted(lam_obs["donor_id"].unique().tolist())
    leave_one_out_rows: list[dict] = []
    for excluded in lam_donors:
        remaining = lam_obs[lam_obs["donor_id"] != excluded]
        row = {"excluded_donor": excluded, "remaining_cells": int(len(remaining))}
        for column in ["lamcore_score_formal", "lamcore_score_fallback", *program_columns]:
            row[f"mean_{column}_formal_candidates"] = float(remaining.loc[remaining["lamcore_candidate_formal"], column].mean())
            row[f"mean_{column}_fallback_candidates"] = float(remaining.loc[remaining["lamcore_candidate_fallback"], column].mean())
        leave_one_out_rows.append(row)
    pd.DataFrame(leave_one_out_rows).to_csv(table_dir / "leave_one_donor_out.csv", index=False)

    write_pseudobulk(adata, "lamcore_candidate_formal", table_dir / "donor_pseudobulk_counts.csv")
    write_pseudobulk(adata, "lamcore_candidate_fallback", table_dir / "donor_pseudobulk_counts_fallback.csv")

    candidate_comparison = (
        adata.obs.groupby(["condition", "donor_id", "sample_id"], observed=True)
        .agg(
            formal_candidates=("lamcore_candidate_formal", "sum"),
            formal_cells=("lamcore_candidate_formal", "size"),
            fallback_candidates=("lamcore_candidate_fallback", "sum"),
        )
        .reset_index()
    )
    candidate_comparison["formal_fraction"] = candidate_comparison["formal_candidates"] / candidate_comparison["formal_cells"]
    candidate_comparison["fallback_fraction"] = candidate_comparison["fallback_candidates"] / candidate_comparison["formal_cells"]
    candidate_comparison.to_csv(table_dir / "lamcore_candidate_comparison.csv", index=False)
    overlap = int((formal_candidate & fallback_candidate).sum())
    comparison_metadata = {
        "formal_candidate_cells": formal_candidate_count,
        "fallback_candidate_cells": fallback_candidate_count,
        "overlap_cells": overlap,
        "formal_fraction_of_all_cells": formal_candidate_count / adata.n_obs,
        "fallback_fraction_of_all_cells": fallback_candidate_count / adata.n_obs,
        "overlap_fraction_of_formal": overlap / formal_candidate_count if formal_candidate_count else np.nan,
        "overlap_fraction_of_fallback": overlap / fallback_candidate_count if fallback_candidate_count else np.nan,
    }
    (table_dir / "lamcore_candidate_comparison.json").write_text(json.dumps(comparison_metadata, indent=2))

    stability_rows: list[dict] = []
    stability_labels: dict[int, np.ndarray] = {}
    for stability_seed in config["preprocess"]["stability_seeds"]:
        key = f"leiden_stability_{stability_seed}"
        sc.tl.leiden(
            adata,
            resolution=float(config["preprocess"]["primary_leiden_resolution"]),
            random_state=int(stability_seed),
            key_added=key,
        )
        labels = adata.obs[key].astype(str).to_numpy()
        stability_labels[int(stability_seed)] = labels
        sample_size = min(5000, adata.n_obs)
        rng = np.random.default_rng(seed + int(stability_seed))
        indices = rng.choice(adata.n_obs, size=sample_size, replace=False)
        stability_rows.append(
            {
                "seed": int(stability_seed),
                "n_clusters": int(pd.Series(labels).nunique()),
                "silhouette": float(silhouette_score(adata.obsm["X_pca"][indices], labels[indices])),
            }
        )
    for left_seed, right_seed in itertools.combinations(sorted(stability_labels), 2):
        stability_rows.append(
            {
                "seed": f"{left_seed}_vs_{right_seed}",
                "n_clusters": np.nan,
                "silhouette": np.nan,
                "ari": float(adjusted_rand_score(stability_labels[left_seed], stability_labels[right_seed])),
                "nmi": float(normalized_mutual_info_score(stability_labels[left_seed], stability_labels[right_seed])),
            }
        )
    pd.DataFrame(stability_rows).to_csv(table_dir / "clustering_stability.csv", index=False)

    for color, filename in [
        ("lamcore_score_formal", "umap_lamcore_score_formal.png"),
        ("lamcore_score_fallback", "umap_lamcore_score_fallback.png"),
        ("lamcore_label", "umap_lamcore_label_formal.png"),
        ("lamcore_label_fallback", "umap_lamcore_label_fallback.png"),
        ("program_contractile", "umap_contractile.png"),
        ("program_ecM_remodeling", "umap_ecm_remodeling.png"),
        ("program_proliferation", "umap_proliferation.png"),
        ("program_stress_hypoxia", "umap_stress_hypoxia.png"),
    ]:
        write_umap(adata, color, figure_dir / filename)

    adata.uns["lamcore"] = {
        "formal_signature": formal_file_metadata,
        "fallback_signature": {
            "source": signatures["lamcore_fallback"].get("source"),
            "genes": fallback_genes,
            "available": signature_status["lamcore_fallback"]["available"],
        },
        "formal_candidate_method": formal_method,
        "fallback_candidate_method": fallback_method,
        "formal_candidate_cutoff": formal_cutoff,
        "fallback_candidate_cutoff": fallback_cutoff,
        "primary_cluster": primary_cluster,
        "signature_status": signature_status,
        "score_method": "scanpy.tl.score_genes on normalized log1p expression with control genes; not an exact reimplementation of the original ssGSEA/classification threshold",
        "candidate_comparison": comparison_metadata,
    }
    adata.uns["stability"] = {
        "lam_donors": lam_donors,
        "formal_passes_any_3_of_4_program_rule": bool(formal_concordance["passes_3_of_4_rule"].any()),
        "fallback_passes_any_3_of_4_program_rule": bool(fallback_concordance["passes_3_of_4_rule"].any()),
    }
    adata.write_h5ad(output_path, compression="gzip")

    program_names_zh = {
        "program_contractile": "收缩/平滑肌样",
        "program_ecM_remodeling": "ECM重塑",
        "program_proliferation": "增殖",
        "program_stress_hypoxia": "应激/缺氧",
        "program_inflammatory_response": "炎症/免疫响应",
    }
    concordance_zh = formal_concordance.assign(
        program=formal_concordance["program"].map(program_names_zh).fillna(formal_concordance["program"])
    ).rename(
        columns={
            "program": "程序",
            "n_lam_donors_tested": "可比较LAM供体数",
            "n_candidate_higher": "候选细胞较高的供体数",
            "concordance_fraction": "方向一致比例",
            "passes_observed_donor_rule": "通过可观察供体规则",
            "passes_3_of_4_rule": "通过严格3/4规则",
        }
    )
    formal_by_donor = candidate_comparison[candidate_comparison["condition"].astype(str) == "LAM"]
    formal_usable_donors = int((formal_by_donor["formal_candidates"] >= int(analysis_cfg["min_cells_per_donor_for_summary"])).sum())
    fallback_usable_donors = int((formal_by_donor["fallback_candidates"] >= int(analysis_cfg["min_cells_per_donor_for_summary"])).sum())
    lam2_formal = int(candidate_comparison.loc[candidate_comparison["donor_id"].astype(str) == "LAM2", "formal_candidates"].sum())
    lam2_fallback = int(candidate_comparison.loc[candidate_comparison["donor_id"].astype(str) == "LAM2", "fallback_candidates"].sum())
    donor1_formal = int(candidate_comparison.loc[candidate_comparison["donor_id"].astype(str) == "Donor1", "formal_candidates"].sum())
    formal_passed_programs = formal_concordance.loc[
        formal_concordance["passes_3_of_4_rule"], "program"
    ].map(program_names_zh).tolist()
    formal_failed_programs = formal_concordance.loc[
        ~formal_concordance["passes_3_of_4_rule"], "program"
    ].map(program_names_zh).tolist()
    program_overlap_rows = []
    formal_gene_upper = {gene.upper() for gene in formal_genes}
    for name, genes in signatures["programs"].items():
        overlap_genes = sorted({gene.upper() for gene in genes} & formal_gene_upper)
        program_overlap_rows.append(
            {
                "program": f"program_{name}",
                "program_genes": len(genes),
                "overlap_with_formal_lamcore_signature": len(overlap_genes),
                "overlap_genes": ";".join(overlap_genes),
            }
        )
    pd.DataFrame(program_overlap_rows).to_csv(table_dir / "lamcore_signature_program_overlap.csv", index=False)
    report_lines_en = [
        "# LAM Cell-State Research Report",
        "",
        f"- Accession: `{accession}`",
        f"- Cells analyzed: **{adata.n_obs:,}**",
        f"- Genes retained in the counts contract: **{adata.n_vars:,}**",
        f"- Formal LAMCORE signature genes: **{len(formal_genes):,}**; available in this matrix: **{len(signature_status['lamcore_formal']['available']):,}**",
        f"- Formal LAMCORE candidate method: `{formal_method}`",
        f"- Formal LAMCORE candidate cells: **{formal_candidate_count:,}**",
        f"- Exploratory four-gene candidate cells: **{fallback_candidate_count:,}**",
        f"- Overlap between the two candidate sets: **{overlap:,}** cells",
        "",
        "## What was reproduced",
        "",
        "The official LAM Cell Atlas table was downloaded and its 777 unique LAMCORE genes were scored in the GSE135851 cells. Scores were calculated with Scanpy's control-gene module-score implementation on normalized log1p expression. This reproduces the published gene set and provides a transparent re-scoring, but it is not a byte-for-byte replication of the original study's cell labels or scoring threshold.",
        "The four-gene panel (ACTA2, PMEL, FIGF, MLANA) is retained only as an exploratory comparison, not as a substitute for the formal signature.",
        "",
        "## Interpretation guardrails",
        "",
        "This is an exploratory analysis. Cell-level tests are not treated as donor-level evidence.",
        f"A cross-donor state is considered stable only when the pre-specified donor-wise rule is met in all four LAM donors. The formal candidate set has usable candidate counts in {formal_usable_donors}/4 LAM donors; the four-gene comparison has usable counts in {fallback_usable_donors}/4.",
        "Because candidate cells are selected using the formal LAMCORE signature, higher scores for signature-overlapping programs are not independent validation. The program/signature overlap is recorded in `lamcore_signature_program_overlap.csv`.",
        "",
        "## Formal-signature donor-level program concordance",
        "",
        "```text",
        formal_concordance.to_string(index=False),
        "```",
        "",
        "## Formal signature versus exploratory panel",
        "",
        f"LAM2 has {lam2_formal} formal-signature candidates versus {lam2_fallback} four-gene candidates. Donor1 has {donor1_formal} formal-signature candidates. These counts show why the candidate definition must be reported explicitly and why neither set alone is treated as a cell-type label.",
        "",
        "## Robustness checks",
        "",
        "Threshold sensitivity, loose/primary/strict QC sensitivity, and assay-stratified summaries are generated by `scripts/robustness_tests.py`. These checks assess sensitivity of the operational candidate definition; they do not replace validation in an independent cohort.",
        "",
        "## Conclusion",
        "",
        f"Under the formal 777-gene candidate rule, the following programs meet the pre-specified directional 3/4 donor rule: {', '.join(formal_passed_programs) if formal_passed_programs else 'none'}. Programs not meeting it: {', '.join(formal_failed_programs) if formal_failed_programs else 'none'}. This is evidence for reproducible expression programs under this operational definition, not proof of a new LAMCORE subtype. The candidate fraction varies substantially by donor, LAM2 contributes {lam2_formal} formal candidates, and the score/candidate definition is not independent of signature-overlapping programs; independent validation is still required.",
        "The original Guo et al. study reported no LAMCORE cells in LAM2 under its own classification procedure. Our 28 LAM2 candidates therefore must not be read as a contradiction: this run uses the published gene set but a new, transparent Scanpy score and cluster threshold rather than the original cell-level classification.",
        "",
        "## Data and run outputs",
        "",
        f"- Formal signature file: `data/raw/reference/LAM_core_signature_genes.csv` (SHA-256 `{formal_file_metadata['sha256']}`)",
        f"- Processed object: `{output_path.relative_to(ROOT)}`",
        f"- Tables: `{table_dir.relative_to(ROOT)}`",
        f"- Figures: `{figure_dir.relative_to(ROOT)}`",
        f"- Official source: {signatures['lamcore_formal'].get('source_url')}",
    ]
    report_lines_zh = [
        "# LAM 细胞状态研究报告",
        "",
        f"- 数据集：`{accession}`",
        f"- 纳入分析的细胞数：**{adata.n_obs:,}**",
        f"- counts 数据契约中的基因数：**{adata.n_vars:,}**",
        f"- 正式 LAMCORE signature 基因数：**{len(formal_genes):,}**；在本矩阵中可匹配：**{len(signature_status['lamcore_formal']['available']):,}**",
        f"- 正式 LAMCORE 候选识别方法：`{formal_method}`",
        f"- 正式 LAMCORE 候选细胞数：**{formal_candidate_count:,}**",
        f"- 四基因探索版候选细胞数：**{fallback_candidate_count:,}**",
        f"- 两种候选集合的重叠细胞数：**{overlap:,}**",
        "",
        "## 这次真正复现了什么",
        "",
        "本轮从 LAM Cell Atlas 官方页面下载了 777 个不重复的 LAMCORE 基因，并在 GSE135851 细胞上用完整基因表重新计算分数。分数采用 Scanpy 的 control-gene module score，输入是标准化并 log1p 后的表达矩阵。因此，本轮复现了正式基因集合和透明的再评分过程，但不是对原论文细胞标签或分类阈值的逐字节复制。",
        "ACTA2、PMEL、FIGF、MLANA 四基因只保留作探索性对照，不能替代正式 signature。",
        "",
        "## 结果解释边界",
        "",
        "本分析属于探索性研究；单细胞层面的统计不作为供体级证据。",
        f"只有四位 LAM 供体都满足预先设定的供体级规则，才能把某个状态称为跨供体稳定。正式 signature 候选细胞数达到最低要求的供体为 {formal_usable_donors}/4；四基因对照为 {fallback_usable_donors}/4。由于候选细胞本身是用正式 LAMCORE signature 选出来的，与 signature 重叠的程序分数升高不构成独立验证；重叠关系记录在 `lamcore_signature_program_overlap.csv`。",
        "",
        "## 正式 signature 的供体级表达程序一致性",
        "",
        "```text",
        concordance_zh.to_string(index=False),
        "```",
        "",
        "## 正式 signature 与探索版对照",
        "",
        f"LAM2 中正式 signature 候选细胞为 {lam2_formal} 个，四基因探索版为 {lam2_fallback} 个；Donor1 中正式 signature 候选细胞为 {donor1_formal} 个。这说明候选定义必须明确报告，任何一种候选集合都不能直接当作细胞类型标签。",
        "",
        "## 稳健性检查",
        "",
        "阈值敏感性、宽松/主分析/严格 QC 敏感性以及按 assay 分层的汇总由 `scripts/robustness_tests.py` 生成。这些检查用于评估当前候选定义是否对参数敏感，不能替代独立队列验证。",
        "",
        "## 结论",
        "",
        f"按照正式 777 基因候选规则，以下程序满足预设的供体方向 3/4 规则：{'、'.join(formal_passed_programs) if formal_passed_programs else '无'}；未满足的程序：{'、'.join(formal_failed_programs) if formal_failed_programs else '无'}。这说明在当前操作性定义下，部分表达程序具有可重复方向，但不能据此证明存在新的 LAMCORE 亚型。候选比例在供体之间差异很大，LAM2 只有 {lam2_formal} 个正式候选，而且候选定义与部分程序存在 signature 基因重叠，因此仍需独立数据验证。",
        "Guo 等人的原始研究使用其自己的分类流程，在 LAM2 中报告没有 LAMCORE 细胞。因此，本轮出现 28 个 LAM2 候选并不等于推翻原结果：本轮使用了发表的基因集合，但采用了新的、透明的 Scanpy 分数和 cluster 阈值，并没有复刻原始细胞级分类阈值。",
        "",
        "## 数据与运行输出",
        "",
        f"- 正式 signature 文件：`data/raw/reference/LAM_core_signature_genes.csv`（SHA-256 `{formal_file_metadata['sha256']}`）",
        f"- 处理后的 AnnData：`{output_path.relative_to(ROOT)}`",
        f"- 结果表：`{table_dir.relative_to(ROOT)}`",
        f"- 图形：`{figure_dir.relative_to(ROOT)}`",
        f"- 官方来源：{signatures['lamcore_formal'].get('source_url')}",
    ]
    (report_dir / "LAM_state_report_en.md").write_text("\n".join(report_lines_en) + "\n")
    (report_dir / "LAM_state_report_zh.md").write_text("\n".join(report_lines_zh) + "\n")
    (report_dir / "LAM_state_report.md").write_text("\n".join(report_lines_en) + "\n")

    run_manifest_path = ROOT / config["paths"]["manifests"] / "run_manifest.yaml"
    run_manifest = yaml.safe_load(run_manifest_path.read_text()) if run_manifest_path.exists() else {}
    run_manifest["steps"] = [step for step in run_manifest.get("steps", []) if step.get("name") != "analyze_lam_states"]
    run_manifest["steps"].append(
        {
            "name": "analyze_lam_states",
            "completed_at": datetime.now(timezone.utc).isoformat(),
            "input": str(input_path.relative_to(ROOT)),
            "output": str(output_path.relative_to(ROOT)),
            "formal_signature_file": formal_file_metadata["file"],
            "formal_signature_sha256": formal_file_metadata["sha256"],
            "formal_candidate_method": formal_method,
            "fallback_candidate_method": fallback_method,
        }
    )
    run_manifest["status"] = "analysis_completed_formal_signature"
    run_manifest_path.write_text(yaml.safe_dump(run_manifest, sort_keys=False))
    print(
        json.dumps(
            {
                "output": str(output_path),
                "formal_candidate_method": formal_method,
                "formal_candidate_cells": formal_candidate_count,
                "fallback_candidate_method": fallback_method,
                "fallback_candidate_cells": fallback_candidate_count,
                "formal_signature_genes": len(formal_genes),
                "lam_donors": lam_donors,
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"analyze_lam_states.py failed: {exc}", file=sys.stderr)
        raise
