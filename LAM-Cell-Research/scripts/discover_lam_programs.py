#!/usr/bin/env python3
"""Discover candidate LAM programs without closing the search space first.

This script is deliberately exploratory.  It maintains a conservative
high-confidence pool, a broad marker-supported pool, and an unrestricted LAM
pool as a guardrail against removing weak/noncanonical states before discovery.
Pooled NMF proposes programs; donor-wise NMF and meta-program matching provide
the first stability check.  Known programs are compared after discovery and
are never regressed out of the primary matrix.
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
import zlib
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy.optimize import linear_sum_assignment
from scipy.sparse import issparse
from sklearn.decomposition import NMF
from sklearn.linear_model import LinearRegression
from sklearn.metrics import adjusted_rand_score

warnings.filterwarnings("ignore", category=FutureWarning)


MARKER_COLUMNS = [
    "marker_expr_PMEL",
    "marker_expr_MLANA",
    "marker_expr_MITF",
    "marker_expr_ACTA2",
    "marker_expr_ESR1",
    "marker_expr_FIGF",
    "marker_expr_VEGFD",
    "marker_expr_CTSK",
]
MARKER_NAMES = [x.replace("marker_expr_", "") for x in MARKER_COLUMNS]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--config", default="config/program_discovery.yaml")
    p.add_argument("--input", default=None)
    p.add_argument("--output-dir", default=None)
    p.add_argument("--fast", action="store_true", help="small deterministic smoke run")
    p.add_argument("--skip-external", action="store_true")
    return p.parse_args()


def load_yaml(path: Path) -> dict:
    with path.open() as handle:
        return yaml.safe_load(handle) or {}


def as_dense(x) -> np.ndarray:
    if issparse(x):
        return x.toarray().astype(np.float32, copy=False)
    return np.asarray(x, dtype=np.float32)


def gene_index(adata: ad.AnnData) -> dict[str, int]:
    symbols = adata.var.get("gene_symbol_upper", pd.Series(adata.var_names, index=adata.var_names))
    result: dict[str, int] = {}
    for i, symbol in enumerate(symbols.astype(str).str.upper()):
        result.setdefault(symbol, i)
    for i, name in enumerate(adata.var_names.astype(str)):
        result.setdefault(name.upper(), i)
    return result


def resolve_genes(names: Iterable[str], lookup: dict[str, int]) -> list[int]:
    return [lookup[g.upper()] for g in names if g and g.upper() in lookup]


def stable_factor_similarity(w1: np.ndarray, w2: np.ndarray) -> float:
    a = w1 / (np.linalg.norm(w1, axis=0, keepdims=True) + 1e-8)
    b = w2 / (np.linalg.norm(w2, axis=0, keepdims=True) + 1e-8)
    sim = np.clip(a.T @ b, -1.0, 1.0)
    rows, cols = linear_sum_assignment(-sim)
    return float(sim[rows, cols].mean()) if len(rows) else 0.0


def sample_rows(obs: pd.DataFrame, rows: np.ndarray, max_n: int, seed: int) -> np.ndarray:
    if len(rows) <= max_n:
        return np.sort(np.asarray(rows, dtype=int))
    rng = np.random.default_rng(seed)
    return np.sort(rng.choice(rows, size=max_n, replace=False).astype(int))


def build_candidate_pools(obs: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=obs.index.astype(str))
    out["cell_id"] = obs.index.astype(str)
    out["condition"] = obs.get("condition", "unknown").astype(str).values
    out["donor_id"] = obs.get("donor_id", "unknown").astype(str).values
    out["assay"] = obs.get("assay", "unknown").astype(str).values
    lam = out["condition"].eq("LAM").to_numpy()

    values = []
    for c in MARKER_COLUMNS:
        if c in obs:
            values.append(pd.to_numeric(obs[c], errors="coerce").fillna(0).to_numpy() > 0)
    marker_count = np.column_stack(values).sum(axis=1) if values else np.zeros(len(obs), dtype=int)
    combo = obs.get("known_marker_combo_ge2", pd.Series(False, index=obs.index)).astype(bool).to_numpy()
    author = obs.get("lamcore_candidate_author_style", pd.Series(False, index=obs.index)).astype(bool).to_numpy()
    formal = obs.get("lamcore_candidate_formal", pd.Series(False, index=obs.index)).astype(bool).to_numpy()
    known_detected = pd.to_numeric(obs.get("known_marker_genes_detected", 0), errors="coerce").fillna(0).to_numpy()

    high = lam & (author | formal | (combo & (marker_count >= 2)))
    broad = lam & ((marker_count >= 1) | (known_detected >= 1) | combo)
    unrestricted = lam.copy()
    out["classical_marker_count"] = marker_count
    out["known_marker_genes_detected"] = known_detected
    out["source_author_style"] = author
    out["source_formal_signature"] = formal
    out["pool_high_confidence"] = high
    out["pool_broad_lam_like"] = broad
    out["pool_unrestricted_lam"] = unrestricted
    out["identity_status"] = np.select(
        [high, broad, unrestricted],
        ["high_confidence_lamcore_like", "provisional_broad_lam_like", "unrestricted_lam_guardrail"],
        default="not_in_lam_discovery_pool",
    )
    out["candidate_reason"] = np.select(
        [author | formal, combo & (marker_count >= 2), marker_count >= 1, known_detected >= 1, unrestricted],
        ["author_or_formal_candidate", "known_marker_combo_and_marker_support", "weak_classical_marker", "known_marker_detection", "LAM_condition_guardrail"],
        default="none",
    )
    return out.reset_index(drop=True)


def select_feature_genes(adata: ad.AnnData, lookup: dict[str, int], known: dict, max_genes: int) -> tuple[list[int], list[str]]:
    var = adata.var.copy()
    symbols = var.get("gene_symbol_upper", pd.Series(adata.var_names, index=adata.var_names)).astype(str).str.upper()
    exclude = symbols.str.match(r"^(MT-|RPS|RPL)")
    hvg = var.get("highly_variable", pd.Series(False, index=var.index)).astype(bool).to_numpy()
    if hvg.sum() < 100:
        hvg = ~exclude.to_numpy()
    else:
        hvg = hvg & ~exclude.to_numpy()
    candidate = np.flatnonzero(hvg)
    if len(candidate) > max_genes:
        dispersion_series = var.get("dispersions_norm", var.get("dispersions"))
        disp = None
        if dispersion_series is not None:
            disp = pd.to_numeric(dispersion_series, errors="coerce").fillna(-np.inf).to_numpy()
            if np.all(~np.isfinite(disp)) or np.nanmax(disp) == 0:
                disp = None
        if disp is None:
            # External matrices often do not carry HVG annotations. Estimate
            # variance from a deterministic cell sample without densifying the
            # complete matrix, so feature selection remains expression-driven.
            rng = np.random.default_rng(20260822)
            n_sample = min(2000, adata.n_obs)
            sample_rows = np.sort(rng.choice(np.arange(adata.n_obs), size=n_sample, replace=False))
            block = adata.X[sample_rows, candidate]
            if issparse(block):
                mean = np.asarray(block.mean(axis=0)).ravel()
                mean_sq = np.asarray(block.power(2).mean(axis=0)).ravel()
            else:
                block = np.asarray(block, dtype=np.float32)
                mean = block.mean(axis=0)
                mean_sq = (block * block).mean(axis=0)
            disp = np.zeros(len(var), dtype=float)
            disp[candidate] = np.maximum(mean_sq - mean * mean, 0.0)
        candidate = candidate[np.argsort(disp[candidate])[-max_genes:]]
    wanted = set()
    for item in known.get("programs", []):
        wanted.update(g.upper() for g in item.get("genes", []))
    wanted_idx = resolve_genes(wanted, lookup)
    selected = list(dict.fromkeys([*candidate.tolist(), *wanted_idx]))
    if len(selected) > max_genes:
        # Keep all resolved comparison genes where possible, then fill with HVGs.
        selected = list(dict.fromkeys(wanted_idx + candidate.tolist()))[:max_genes]
    names = [str(symbols.iloc[i]) for i in selected]
    return selected, names


def read_matrix(adata: ad.AnnData, rows: np.ndarray, cols: list[int]) -> np.ndarray:
    if len(rows) == 0 or len(cols) == 0:
        return np.zeros((len(rows), len(cols)), dtype=np.float32)
    return as_dense(adata.X[rows, cols])


def run_nmf_replicates(
    x: np.ndarray,
    genes: list[str],
    ranks: list[int],
    seeds: list[int],
    max_iter: int,
    tol: float,
    label: str,
) -> tuple[pd.DataFrame, dict[tuple[int, int], tuple[np.ndarray, np.ndarray]]]:
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    x[x < 0] = 0
    keep = np.isfinite(x).all(axis=0) & (x.var(axis=0) > 1e-10)
    x = x[:, keep]
    genes = [g for g, k in zip(genes, keep) if k]
    if x.shape[1] < 5 or x.shape[0] < 5:
        return pd.DataFrame(), {}
    x_norm = float(np.linalg.norm(x)) + 1e-8
    rows = []
    models: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]] = {}
    for rank in ranks:
        if rank >= min(x.shape):
            continue
        for seed in seeds:
            model = NMF(
                n_components=rank,
                # Use randomized initialization so the configured seed sweep
                # measures real factor stability rather than repeating the
                # same deterministic NNDSVDA start.
                init="random",
                solver="cd",
                beta_loss="frobenius",
                max_iter=max_iter,
                tol=tol,
                random_state=seed,
            )
            w = model.fit_transform(x)
            h = model.components_
            models[(rank, seed)] = (w, h)
            rows.append(
                {
                    "pool": label,
                    "rank": rank,
                    "seed": seed,
                    "n_cells": x.shape[0],
                    "n_genes": x.shape[1],
                    "reconstruction_error": float(model.reconstruction_err_),
                    "normalized_reconstruction_error": float(model.reconstruction_err_ / x_norm),
                    "converged_or_max_iter": int(model.n_iter_),
                }
            )
    summary = pd.DataFrame(rows)
    if summary.empty:
        return summary, models
    stabilities = []
    for rank in summary["rank"].unique():
        runs = [(s, models[(rank, s)][0]) for s in seeds if (rank, s) in models]
        values = [stable_factor_similarity(a, b) for i, (_, a) in enumerate(runs) for _, b in runs[i + 1 :]]
        stabilities.append((rank, float(np.mean(values)) if values else 0.0))
    stability_map = dict(stabilities)
    summary["factor_stability"] = summary["rank"].map(stability_map).fillna(0.0)
    summary["selection_score"] = 0.7 * summary["factor_stability"] + 0.3 / (1.0 + summary["normalized_reconstruction_error"])
    return summary, models


def top_genes(h: np.ndarray, genes: list[str], top_n: int) -> list[tuple[str, float, int]]:
    order = np.argsort(h)[::-1][:top_n]
    return [(genes[i], float(h[i]), rank + 1) for rank, i in enumerate(order)]


def best_model(summary: pd.DataFrame, models: dict[tuple[int, int], tuple[np.ndarray, np.ndarray]]) -> tuple[tuple[int, int], np.ndarray, np.ndarray]:
    row = summary.sort_values(["selection_score", "factor_stability", "normalized_reconstruction_error"], ascending=[False, False, True]).iloc[0]
    key = (int(row["rank"]), int(row["seed"]))
    w, h = models[key]
    return key, w, h


def known_comparisons(
    pool: str,
    w: np.ndarray,
    top_gene_lists: dict[int, list[str]],
    x: np.ndarray,
    feature_genes: list[str],
    known: dict,
) -> pd.DataFrame:
    rows = []
    gene_to_idx = {g.upper(): i for i, g in enumerate(feature_genes)}
    known_scores = {}
    for item in known.get("programs", []):
        name = item["program_name"]
        genes = [g.upper() for g in item.get("genes", []) if g.upper() in gene_to_idx]
        if not genes:
            continue
        idx = [gene_to_idx[g] for g in genes]
        known_scores[name] = x[:, idx].mean(axis=1)
        for program_id, top in top_gene_lists.items():
            top_set = {g.upper() for g in top}
            overlap = sorted(top_set.intersection(genes))
            union = top_set.union(genes)
            corr = np.corrcoef(w[:, program_id], known_scores[name])[0, 1] if len(overlap) or np.std(known_scores[name]) > 0 else np.nan
            rows.append(
                {
                    "pool": pool,
                    "candidate_program": f"program_{program_id + 1}",
                    "known_program": name,
                    "evidence_scope": item.get("evidence_scope"),
                    "evidence_level": item.get("evidence_level"),
                    "known_type": item.get("type"),
                    "overlap_n": len(overlap),
                    "overlap_genes": ",".join(overlap),
                    "jaccard_top_genes": len(overlap) / len(union) if union else 0.0,
                    "score_correlation": float(corr) if np.isfinite(corr) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def core3_scores(adata: ad.AnnData, lam_rows: np.ndarray, lookup: dict[str, int], known: dict) -> pd.DataFrame:
    model_cfg = known.get("core3_model", {})
    ident_genes = next((p.get("genes", []) for p in known.get("programs", []) if p.get("program_name") == model_cfg.get("identity_program")), [])
    trans_genes = model_cfg.get("translation_program", [])
    ident_idx = resolve_genes(ident_genes, lookup)
    trans_idx = resolve_genes(trans_genes, lookup)
    all_idx = list(dict.fromkeys(ident_idx + trans_idx))
    if not all_idx:
        return pd.DataFrame()
    x = read_matrix(adata, lam_rows, all_idx)
    ident = x[:, : len(ident_idx)].mean(axis=1) if ident_idx else np.zeros(len(lam_rows))
    trans_start = len(ident_idx)
    trans = x[:, trans_start:].mean(axis=1) if trans_idx else np.zeros(len(lam_rows))
    obs = adata.obs.iloc[lam_rows].copy()
    depth = np.log1p(pd.to_numeric(obs.get("n_genes_by_counts", 0), errors="coerce").fillna(0).to_numpy())
    counts = np.log1p(pd.to_numeric(obs.get("total_counts", 0), errors="coerce").fillna(0).to_numpy())
    assay = pd.get_dummies(obs.get("assay", "unknown").astype(str), dtype=float).to_numpy()
    design = np.column_stack([depth, counts, assay])
    residual = ident - LinearRegression().fit(design, ident).predict(design) if len(ident) > design.shape[1] + 2 else ident
    ident_q = np.quantile(ident, 0.75) if len(ident) else 0
    low_q = np.quantile(residual, 0.25) if len(residual) else 0
    trans_q = np.quantile(trans, 0.75) if len(trans) else 0
    out = pd.DataFrame(
        {
            "cell_id": adata.obs.index.astype(str)[lam_rows],
            "donor_id": obs.get("donor_id", "unknown").astype(str).to_numpy(),
            "assay": obs.get("assay", "unknown").astype(str).to_numpy(),
            "core3_identity_score": ident,
            "depth_adjusted_identity_activity": residual,
            "translation_enrichment_score": trans,
            "n_genes_by_counts": obs.get("n_genes_by_counts", 0).to_numpy(),
            "total_counts": obs.get("total_counts", 0).to_numpy(),
        }
    )
    out["core3_like"] = (out["core3_identity_score"] >= ident_q) & (out["depth_adjusted_identity_activity"] <= low_q) & (out["translation_enrichment_score"] >= trans_q)
    out["core3_model_note"] = "identity + depth-adjusted low activity + translation; low complexity alone is insufficient"
    return out


def residual_sensitivity(pool: str, w: np.ndarray, known_scores: pd.DataFrame, obs: pd.DataFrame) -> pd.DataFrame:
    if known_scores.empty or w.size == 0:
        return pd.DataFrame()
    covars = known_scores.to_numpy(dtype=float)
    assay = pd.get_dummies(obs.get("assay", "unknown").astype(str), dtype=float).to_numpy()
    depth = np.column_stack([
        np.log1p(pd.to_numeric(obs.get("n_genes_by_counts", 0), errors="coerce").fillna(0).to_numpy()),
        np.log1p(pd.to_numeric(obs.get("total_counts", 0), errors="coerce").fillna(0).to_numpy()),
        assay,
    ])
    design = np.column_stack([covars, depth])
    rows = []
    for j in range(w.shape[1]):
        target = w[:, j]
        if len(target) <= design.shape[1] + 2:
            continue
        model = LinearRegression().fit(design, target)
        rows.append({"pool": pool, "candidate_program": f"program_{j + 1}", "known_program_explained_r2": float(model.score(design, target)), "residual_sd": float(np.std(target - model.predict(design)))})
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    root = Path(__file__).resolve().parents[1]
    cfg = load_yaml(root / args.config)
    known = load_yaml(root / cfg["known_programs"])
    input_path = root / (args.input or cfg["input_h5ad"])
    outdir = root / (args.output_dir or cfg["output_dir"])
    outdir.mkdir(parents=True, exist_ok=True)
    seed = int(cfg.get("seed", 20260822))
    nmf_cfg = cfg.get("nmf", {})
    ranks = [2, 3, 4, 5, 6]
    seeds = [0, 1, 2, 3, 4]
    max_genes = int(nmf_cfg.get("max_genes", 2000))
    max_pool = int(nmf_cfg.get("max_cells_pooled", 6000))
    max_donor = int(nmf_cfg.get("max_cells_per_donor", 1500))
    max_iter = int(nmf_cfg.get("max_iter", 500))
    tol = float(nmf_cfg.get("tol", 1e-4))
    if args.fast:
        ranks, seeds, max_genes, max_pool, max_donor, max_iter = [2, 3], [0, 1], 800, 1200, 600, 200

    adata = ad.read_h5ad(input_path, backed="r")
    obs = adata.obs.copy()
    pools = build_candidate_pools(obs)
    pools.to_csv(outdir / "candidate_pool_labels.csv", index=False)
    lookup = gene_index(adata)
    feature_cols, feature_names = select_feature_genes(adata, lookup, known, max_genes)

    all_summary = []
    all_gene_rows = []
    all_score_rows = []
    all_compare = []
    all_residual = []
    donor_gene_rows = []
    donor_match_rows = []
    donor_level_rows = []
    pool_models = {}
    pool_masks = {"high_confidence": pools["pool_high_confidence"].to_numpy(), "broad_lam_like": pools["pool_broad_lam_like"].to_numpy(), "unrestricted_lam": pools["pool_unrestricted_lam"].to_numpy()}

    for pool, mask in pool_masks.items():
        rows = np.flatnonzero(mask)
        rows = sample_rows(obs, rows, max_pool, seed)
        x = read_matrix(adata, rows, feature_cols)
        summary, models = run_nmf_replicates(x, feature_names, ranks, seeds, max_iter, tol, pool)
        if summary.empty:
            continue
        summary["n_source_cells_before_pool_sampling"] = int(mask.sum())
        summary["n_source_donors_before_pool_sampling"] = int(obs.iloc[np.flatnonzero(mask)]["donor_id"].nunique())
        all_summary.append(summary)
        key, w, h = best_model(summary, models)
        pool_models[pool] = {"rows": rows, "x": x, "w": w, "h": h, "genes": feature_names, "key": key}
        top_lists = {}
        for j in range(h.shape[0]):
            top = top_genes(h[j], feature_names, int(nmf_cfg.get("top_genes_per_program", 50)))
            top_lists[j] = [g for g, _, _ in top]
            for gene, weight, rank_pos in top:
                all_gene_rows.append({"pool": pool, "candidate_program": f"program_{j + 1}", "gene": gene, "weight": weight, "rank_position": rank_pos, "selected_rank": key[0], "selected_seed": key[1]})
            for i, row_id in enumerate(rows):
                all_score_rows.append({"pool": pool, "cell_id": str(obs.index[row_id]), "donor_id": str(obs.iloc[row_id].get("donor_id", "unknown")), "assay": str(obs.iloc[row_id].get("assay", "unknown")), "candidate_program": f"program_{j + 1}", "score": float(w[i, j])})
        comp = known_comparisons(pool, w, top_lists, x, feature_names, known)
        if not comp.empty:
            all_compare.append(comp)
        score_frame = pd.DataFrame({p: w[:, j] for j, p in enumerate([f"program_{i + 1}" for i in range(w.shape[1])])})
        known_score_dict = {}
        feature_map = {g.upper(): i for i, g in enumerate(feature_names)}
        for item in known.get("programs", []):
            idx = [feature_map[g.upper()] for g in item.get("genes", []) if g.upper() in feature_map]
            if idx:
                known_score_dict[item["program_name"]] = x[:, idx].mean(axis=1)
        residual = residual_sensitivity(pool, w, pd.DataFrame(known_score_dict), obs.iloc[rows].reset_index(drop=True))
        if not residual.empty:
            all_residual.append(residual)

        selected_obs = obs.iloc[rows].copy()
        for donor, donor_frame in selected_obs.groupby("donor_id", observed=True):
            donor_rows = donor_frame.index
            donor_pos = np.array([obs.index.get_loc(x) for x in donor_rows], dtype=int)
            donor_seed = seed + int(zlib.crc32(str(donor).encode("utf-8")) % 10000)
            donor_pos = sample_rows(obs, donor_pos, max_donor, donor_seed)
            donor_x = read_matrix(adata, donor_pos, feature_cols)
            if len(donor_pos) < int(nmf_cfg.get("min_cells_per_donor_for_nmf", 40)):
                donor_level_rows.append({"pool": pool, "donor_id": str(donor), "n_cells": len(donor_pos), "analysis": "pseudobulk_only_insufficient_cells_for_independent_nmf"})
                continue
            dsummary, dmodels = run_nmf_replicates(donor_x, feature_names, ranks, seeds[: min(3, len(seeds))], max_iter, tol, f"{pool}__{donor}")
            if dsummary.empty:
                continue
            dkey, dw, dh = best_model(dsummary, dmodels)
            for j in range(dh.shape[0]):
                top = top_genes(dh[j], feature_names, int(nmf_cfg.get("top_genes_per_program", 50)))
                for gene, weight, rank_pos in top:
                    donor_gene_rows.append({"pool": pool, "donor_id": str(donor), "candidate_program": f"program_{j + 1}", "gene": gene, "weight": weight, "rank_position": rank_pos, "selected_rank": dkey[0], "selected_seed": dkey[1], "n_cells": len(donor_pos)})
            pooled_top = [set(v) for v in top_lists.values()]
            for j in range(dh.shape[0]):
                donor_top = {g for g, _, _ in top_genes(dh[j], feature_names, int(nmf_cfg.get("top_genes_per_program", 50)))}
                overlaps = [len(donor_top & p) / len(donor_top | p) if donor_top | p else 0 for p in pooled_top]
                best_p = int(np.argmax(overlaps)) if overlaps else -1
                donor_match_rows.append({"pool": pool, "donor_id": str(donor), "donor_program": f"program_{j + 1}", "donor_rank": dkey[0], "best_pooled_program": f"program_{best_p + 1}" if best_p >= 0 else "", "top_gene_jaccard": float(max(overlaps)) if overlaps else 0.0, "independently_discovered": bool(max(overlaps, default=0.0) >= 0.15), "n_cells": len(donor_pos)})
            donor_level_rows.append({"pool": pool, "donor_id": str(donor), "n_cells": len(donor_pos), "analysis": "independent_nmf", "selected_rank": dkey[0], "selected_seed": dkey[1], "n_programs": dh.shape[0]})

        # Donor-level scores for the pooled model, kept separate from cell-level evidence.
        selected_scores = pd.DataFrame(w, columns=[f"program_{j + 1}" for j in range(w.shape[1])])
        selected_scores["donor_id"] = obs.iloc[rows]["donor_id"].astype(str).to_numpy()
        selected_scores["pool"] = pool
        for donor, frame in selected_scores.groupby("donor_id", sort=True):
            row = {"pool": pool, "donor_id": donor, "n_cells": len(frame)}
            row.update({c: float(frame[c].mean()) for c in selected_scores.columns if c.startswith("program_")})
            donor_level_rows.append(row)

    pd.concat(all_summary, ignore_index=True).to_csv(outdir / "pooled_nmf_summary.csv", index=False) if all_summary else pd.DataFrame().to_csv(outdir / "pooled_nmf_summary.csv", index=False)
    pd.DataFrame(all_gene_rows).to_csv(outdir / "pooled_program_genes.csv", index=False)
    pd.DataFrame(all_score_rows).to_csv(outdir / "pooled_program_scores.csv", index=False)
    pd.DataFrame(donor_gene_rows).to_csv(outdir / "donor_wise_program_genes.csv", index=False)
    pd.DataFrame(donor_match_rows).to_csv(outdir / "donor_meta_program_matches.csv", index=False)
    pd.DataFrame(donor_level_rows).to_csv(outdir / "donor_level_program_scores.csv", index=False)
    if donor_match_rows:
        match_df = pd.DataFrame(donor_match_rows)
        meta = (
            match_df[match_df["independently_discovered"]]
            .groupby(["pool", "best_pooled_program"], as_index=False)
            .agg(
                independently_discovered_donors=("donor_id", "nunique"),
                mean_top_gene_jaccard=("top_gene_jaccard", "mean"),
                min_top_gene_jaccard=("top_gene_jaccard", "min"),
                max_top_gene_jaccard=("top_gene_jaccard", "max"),
            )
        )
        meta.to_csv(outdir / "meta_program_summary.csv", index=False)
    else:
        pd.DataFrame().to_csv(outdir / "meta_program_summary.csv", index=False)
    pd.concat(all_compare, ignore_index=True).to_csv(outdir / "known_program_comparisons.csv", index=False) if all_compare else pd.DataFrame().to_csv(outdir / "known_program_comparisons.csv", index=False)
    pd.concat(all_residual, ignore_index=True).to_csv(outdir / "known_program_residual_sensitivity.csv", index=False) if all_residual else pd.DataFrame().to_csv(outdir / "known_program_residual_sensitivity.csv", index=False)

    lam_rows = np.flatnonzero(obs.get("condition", "").astype(str).eq("LAM").to_numpy())
    core3 = core3_scores(adata, lam_rows, lookup, known)
    core3.to_csv(outdir / "core3_structured_scores.csv", index=False)
    pools.loc[pools["pool_broad_lam_like"] & ~pools["pool_high_confidence"], "identity_status"] = "uncertain_identity_hypothesis_pending_validation"
    pools.to_csv(outdir / "candidate_pool_labels.csv", index=False)

    external_status = []
    for accession, relpath in cfg.get("external_inputs", {}).items():
        path = root / relpath
        external_status.append({"accession": accession, "path": str(path.relative_to(root)), "available": path.exists(), "status": "ready_for_external_validation" if path.exists() else "not_downloaded_or_not_prepared"})
    pd.DataFrame(external_status).to_csv(outdir / "external_data_status.csv", index=False)

    summary = {
        "run_mode": "fast_smoke" if args.fast else "configured_full",
        "input": str(input_path.relative_to(root)),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "pool_counts": {k: int(v.sum()) for k, v in pool_masks.items()},
        "pool_donors": {k: int(obs.iloc[np.flatnonzero(v)]["donor_id"].nunique()) for k, v in pool_masks.items()},
        "feature_genes": len(feature_names),
        "pooled_models": {k: {"selected_rank": int(v["key"][0]), "selected_seed": int(v["key"][1]), "n_cells": int(len(v["rows"]))} for k, v in pool_models.items()},
        "external_data": external_status,
        "interpretation": "Candidate programs are exploratory until independent PatientID and orthogonal evidence are available.",
    }
    (outdir / "program_discovery_run_manifest.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False))
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
