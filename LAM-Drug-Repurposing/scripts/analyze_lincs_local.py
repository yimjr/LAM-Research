"""Local LINCS/CMap connectivity analysis for the LAM residual programs.

This module deliberately keeps the CMap-like calculations explicit and
auditable.  It reads the GEO Level 5 GCTX matrices, ranks the complete BING
space for every perturbation signature, computes the weighted KS ES/WTCS,
then applies the published mean-scaling NCS procedure within dataset/query/
cell/perturbagen-type groups.

The two GEO releases are analysed independently.  Their comparison is called
cross-phase/cross-release recurrence, not independent biological replication.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator

import h5py
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from common import CANDIDATE_RESULTS, LINCS_CANDIDATE_INTERMEDIATE, ROOT, write_json


DEFAULT_MATRIX = {
    "GSE92742": {
        "gctx": "data/processed/LINCS/gctx/GSE92742_Level5.gctx",
        "sig_info": "data/raw/LINCS/GSE92742/GSE92742_Broad_LINCS_sig_info.txt.gz",
        "gene_info": "data/raw/LINCS/GSE92742/GSE92742_Broad_LINCS_gene_info.txt.gz",
        "pert_info": "data/raw/LINCS/GSE92742/GSE92742_Broad_LINCS_pert_info.txt.gz",
    },
    "GSE70138": {
        "gctx": "data/processed/LINCS/gctx/GSE70138_Level5.gctx",
        "sig_info": "data/raw/LINCS/GSE70138/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
        "gene_info": "data/raw/LINCS/GSE70138/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz",
        "pert_info": "data/raw/LINCS/GSE70138/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
    },
}

TARGET_SIZES = (50, 100, 150)
NON_PERT_TYPES = {"ctl_vehicle", "ctl_vector", "ctl_untrt", "ctl_vehicle.cns", "ctl_vector.cns", "ctl_untrt.cns"}
CORE_MTOR_NAMES = {"sirolimus", "rapamycin", "everolimus", "temsirolimus", "torin-1", "torin-2"}
CORE_MTOR_GENES = {"mtor", "rptor", "rheb"}
EXTENDED_PI3K_NAMES = {"gdc-0941", "nvp-bez235", "pi-103"}
EXTENDED_PI3K_GENES = {"akt1", "akt2", "akt3", "lamtor3"}


@dataclass(frozen=True)
class Query:
    query_id: str
    contrast: str
    canonical_hash: str
    up_symbols: tuple[str, ...]
    down_symbols: tuple[str, ...]
    up_indices: tuple[int, ...]
    down_indices: tuple[int, ...]
    up_signed_scores: tuple[float, ...]
    down_signed_scores: tuple[float, ...]
    requested_target_sizes: tuple[int, ...]
    actual_target_size: int
    exploratory: bool
    unmapped_up: tuple[str, ...]
    unmapped_down: tuple[str, ...]


def _decode(values: Iterable[object]) -> list[str]:
    result = []
    for value in values:
        if isinstance(value, bytes):
            value = value.decode("utf-8", errors="replace")
        result.append(str(value))
    return result


def _norm_text(value: object) -> str:
    return re.sub(r"\s+", " ", str(value).strip().lower())


def _normalise_perturbation_key(pert_type: str, pert_iname: str) -> tuple[str, str]:
    """Return a cross-release key and a human-readable perturbation class."""
    name = _norm_text(pert_iname)
    if pert_type == "trt_cp":
        return f"compound::{name}", "compound"
    if pert_type.startswith("trt_"):
        return f"genetic::{name}", "genetic"
    return f"other::{pert_type}::{name}", "other"


def _stable_hash(up: Iterable[str], down: Iterable[str]) -> str:
    payload = "up:" + ",".join(sorted(up)) + "|down:" + ",".join(sorted(down))
    return hashlib.sha256(payload.encode()).hexdigest()[:16]


def load_metadata(dataset: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cfg = DEFAULT_MATRIX[dataset]
    sig = pd.read_csv(ROOT / cfg["sig_info"], sep="\t", compression="gzip", low_memory=False)
    gene = pd.read_csv(ROOT / cfg["gene_info"], sep="\t", compression="gzip", low_memory=False)
    pert = pd.read_csv(ROOT / cfg["pert_info"], sep="\t", compression="gzip", low_memory=False)
    sig["sig_id"] = sig["sig_id"].astype(str)
    sig["cell_id"] = sig["cell_id"].astype(str)
    sig["pert_type"] = sig["pert_type"].astype(str)
    sig["pert_iname"] = sig["pert_iname"].astype(str)
    sig["pert_id"] = sig["pert_id"].astype(str)
    sig["perturbation_key"], sig["perturbation_class"] = zip(
        *[_normalise_perturbation_key(t, n) for t, n in zip(sig["pert_type"], sig["pert_iname"])]
    )
    if "pert_idose" in sig.columns:
        sig["dose"] = sig["pert_idose"].astype(str)
    elif "pert_dose" in sig.columns:
        sig["dose"] = sig["pert_dose"].astype(str) + " " + sig.get("pert_dose_unit", "").astype(str)
    else:
        sig["dose"] = ""
    sig["time"] = sig.get("pert_itime", sig.get("pert_time", "")).astype(str)
    pert = pert.drop_duplicates("pert_id")
    return sig, gene, pert


def build_queries(signature_path: Path, gene: pd.DataFrame) -> list[Query]:
    signature = pd.read_csv(signature_path)
    required = {"contrast", "gene", "direction", "signed_score"}
    missing = required - set(signature.columns)
    if missing:
        raise ValueError(f"{signature_path} missing columns: {sorted(missing)}")

    symbol_to_index: dict[str, int] = {}
    bing = gene.loc[gene["pr_is_bing"].astype(int).eq(1)].copy()
    for idx, symbol in enumerate(bing["pr_gene_symbol"].astype(str).str.upper()):
        symbol_to_index.setdefault(symbol, idx)

    queries: dict[str, Query] = {}
    for contrast, group in signature.groupby("contrast", sort=True):
        group = group.copy()
        group["gene"] = group["gene"].astype(str).str.upper().str.strip()
        group = group.drop_duplicates(["direction", "gene"])
        exploratory = str(contrast) == "environment_dependent_escape"
        selected_by_size: dict[int, tuple[list[str], list[str], list[float], list[float]]] = {}
        for target in TARGET_SIZES:
            up = group.loc[group.direction.eq("up")].sort_values(["moderated_q", "signed_score"], na_position="last", ascending=[True, False]).head(target)
            down = group.loc[group.direction.eq("down")].sort_values(["moderated_q", "signed_score"], na_position="last", ascending=[True, True]).head(target)
            selected_by_size[target] = (
                up["gene"].tolist(), down["gene"].tolist(), up["signed_score"].astype(float).tolist(), down["signed_score"].astype(float).tolist()
            )

        grouped_sizes: dict[tuple[tuple[str, ...], tuple[str, ...]], list[int]] = defaultdict(list)
        for target, (up, down, ups, downs) in selected_by_size.items():
            grouped_sizes[(tuple(up), tuple(down))].append(target)
        for (up, down), targets in grouped_sizes.items():
            up_indices = tuple(symbol_to_index[g] for g in up if g in symbol_to_index)
            down_indices = tuple(symbol_to_index[g] for g in down if g in symbol_to_index)
            mapped_up = tuple(g for g in up if g in symbol_to_index)
            mapped_down = tuple(g for g in down if g in symbol_to_index)
            raw_up_scores = selected_by_size[targets[0]][2]
            raw_down_scores = selected_by_size[targets[0]][3]
            up_score_lookup = dict(zip(up, raw_up_scores))
            down_score_lookup = dict(zip(down, raw_down_scores))
            mapped_up_scores = tuple(float(up_score_lookup[g]) for g in mapped_up)
            mapped_down_scores = tuple(float(down_score_lookup[g]) for g in mapped_down)
            canonical_hash = _stable_hash(mapped_up, mapped_down)
            query_id = f"{contrast}__q{len(mapped_up)}x{len(mapped_down)}__{canonical_hash}"
            queries[canonical_hash] = Query(
                query_id=query_id,
                contrast=str(contrast),
                canonical_hash=canonical_hash,
                up_symbols=tuple(mapped_up),
                down_symbols=tuple(mapped_down),
                up_indices=up_indices,
                down_indices=down_indices,
                up_signed_scores=mapped_up_scores,
                down_signed_scores=mapped_down_scores,
                requested_target_sizes=tuple(targets),
                actual_target_size=max(len(mapped_up), len(mapped_down)),
                exploratory=exploratory,
                unmapped_up=tuple(g for g in up if g not in symbol_to_index),
                unmapped_down=tuple(g for g in down if g not in symbol_to_index),
            )
    return sorted(queries.values(), key=lambda q: (q.contrast, q.actual_target_size, q.query_id))


def weighted_ks_from_positions(positions: np.ndarray, weights: np.ndarray, n_genes: int) -> float:
    """CMap/GSEA weighted signed KS ES from zero-based hit positions.

    Positions must identify unique genes in a complete ranked universe.  The
    hit weights are non-negative Level 5 magnitudes; p=1 is the CMap weighted
    KS convention used here.  The implementation evaluates both sides of each
    jump so the sign is not lost at a miss-to-hit boundary.
    """
    if len(positions) == 0 or n_genes <= len(positions):
        return 0.0
    order = np.argsort(positions, kind="mergesort")
    pos = positions[order].astype(np.int64, copy=False)
    w = np.abs(weights[order].astype(float, copy=False))
    total = float(w.sum())
    if not np.isfinite(total) or total <= 0:
        w = np.ones(len(pos), dtype=float)
        total = float(len(pos))
    cum = np.cumsum(w) / total
    hit_before = np.arange(len(pos), dtype=float)
    miss_after = (pos.astype(float) - hit_before) / float(n_genes - len(pos))
    diff_before = np.concatenate(([0.0], cum[:-1])) - miss_after
    diff_after = cum - miss_after
    candidates = np.concatenate((diff_before, diff_after))
    return float(candidates[np.argmax(np.abs(candidates))])


def wtcs(es_up: float, es_down: float) -> float:
    if es_up == 0 or es_down == 0 or es_up * es_down >= 0:
        return 0.0
    return float((es_up - es_down) / 2.0)


def weighted_pearson(x: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y) & np.isfinite(weights) & (weights >= 0)
    x, y, weights = x[mask], y[mask], weights[mask]
    if len(x) < 2 or weights.sum() <= 0:
        return np.nan
    x_bar = np.average(x, weights=weights)
    y_bar = np.average(y, weights=weights)
    dx, dy = x - x_bar, y - y_bar
    denominator = math.sqrt(float(np.sum(weights * dx * dx) * np.sum(weights * dy * dy)))
    return float(np.sum(weights * dx * dy) / denominator) if denominator > 0 else np.nan


def score_rows(values: np.ndarray, queries: list[Query], bing_count: int) -> list[dict[str, np.ndarray]]:
    """Score a matrix chunk shaped (n_signatures, n_BING_genes)."""
    n_rows = values.shape[0]
    # Sorting every complete BING row is intentional: query-only ranking would
    # not be comparable to CMap's rank-space enrichment.
    order = np.argsort(values, axis=1, kind="quicksort")[:, ::-1]
    ranks = np.empty_like(order, dtype=np.int32)
    rank_numbers = np.arange(bing_count, dtype=np.int32)[None, :]
    row_numbers = np.arange(n_rows, dtype=np.int32)[:, None]
    ranks[row_numbers, order] = rank_numbers
    output: list[dict[str, np.ndarray]] = []
    for query in queries:
        up_idx = np.asarray(query.up_indices, dtype=np.int64)
        down_idx = np.asarray(query.down_indices, dtype=np.int64)
        up_score = np.asarray(query.up_signed_scores, dtype=float)
        down_score = np.asarray(query.down_signed_scores, dtype=float)
        es_up = np.zeros(n_rows, dtype=np.float32)
        es_down = np.zeros(n_rows, dtype=np.float32)
        wcs = np.zeros(n_rows, dtype=np.float32)
        corr = np.full(n_rows, np.nan, dtype=np.float32)
        for i in range(n_rows):
            if len(up_idx):
                es_up[i] = weighted_ks_from_positions(ranks[i, up_idx], values[i, up_idx], bing_count)
            if len(down_idx):
                es_down[i] = weighted_ks_from_positions(ranks[i, down_idx], values[i, down_idx], bing_count)
            wcs[i] = wtcs(float(es_up[i]), float(es_down[i]))
            x = np.concatenate((up_score, down_score))
            y = np.concatenate((values[i, up_idx], values[i, down_idx]))
            corr[i] = weighted_pearson(x, y, np.abs(x))
        output.append({"ES_up": es_up, "ES_down": es_down, "WTCS": wcs, "weighted_reversal_correlation": -corr})
    return output


def _metadata_for_matrix(sig: pd.DataFrame, row_ids: list[str]) -> pd.DataFrame:
    lookup = sig.set_index("sig_id")
    missing = [x for x in row_ids if x not in lookup.index]
    if missing:
        raise ValueError(f"GCTX contains {len(missing)} signature IDs missing from sig_info; first={missing[:3]}")
    return lookup.loc[row_ids].reset_index()


def _ncs_denominators(raw_path: Path) -> pd.DataFrame:
    sums: dict[tuple[str, str, str, str], list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0, 0.0])
    parquet = pq.ParquetFile(raw_path)
    columns = ["dataset", "query_id", "cell_id", "pert_type", "WTCS"]
    for batch in parquet.iter_batches(columns=columns, batch_size=200_000):
        frame = batch.to_pandas()
        for key, group in frame.groupby(["dataset", "query_id", "cell_id", "pert_type"], sort=False):
            vals = group["WTCS"].to_numpy(float)
            pos = vals[vals > 0]
            neg = vals[vals < 0]
            item = sums[key]
            item[0] += float(pos.sum())
            item[1] += float(len(pos))
            item[2] += float(neg.sum())
            item[3] += float(len(neg))
    rows = []
    for key, (pos_sum, pos_n, neg_sum, neg_n) in sums.items():
        rows.append({
            "dataset": key[0], "query_id": key[1], "cell_id": key[2], "pert_type": key[3],
            "mu_pos": pos_sum / pos_n if pos_n else np.nan,
            "mu_neg_signed": neg_sum / neg_n if neg_n else np.nan,
            "mu_neg_abs": abs(neg_sum / neg_n) if neg_n else np.nan,
            "n_positive_WTCS": int(pos_n), "n_negative_WTCS": int(neg_n),
        })
    return pd.DataFrame(rows)


def _add_ncs(frame: pd.DataFrame, denominators: pd.DataFrame) -> pd.DataFrame:
    keys = ["dataset", "query_id", "cell_id", "pert_type"]
    frame = frame.merge(denominators, on=keys, how="left", validate="many_to_one")
    # Zero WTCS is the explicit CMap null connectivity value and therefore
    # receives NCS=0 even when a group has no positive/negative denominator.
    frame["NCS"] = 0.0
    positive = frame.WTCS > 0
    negative = frame.WTCS < 0
    frame.loc[positive, "NCS"] = frame.loc[positive, "WTCS"] / frame.loc[positive, "mu_pos"]
    frame.loc[negative, "NCS"] = frame.loc[negative, "WTCS"] / frame.loc[negative, "mu_neg_abs"]
    frame["reversal_WTCS"] = -frame["WTCS"]
    frame["reversal_NCS"] = -frame["NCS"]
    return frame


def run_dataset(dataset: str, queries: list[Query], signature_path: Path, chunk_size: int = 256) -> tuple[Path, Path, dict]:
    cfg = DEFAULT_MATRIX[dataset]
    gctx_path = ROOT / cfg["gctx"]
    if not gctx_path.exists():
        raise FileNotFoundError(f"Missing decompressed GCTX: {gctx_path}")
    sig, gene, pert = load_metadata(dataset)
    bing = gene.loc[gene.pr_is_bing.astype(int).eq(1)].copy()
    with h5py.File(gctx_path, "r") as h5:
        matrix = h5["0/DATA/0/matrix"]
        gene_ids = _decode(h5["0/META/ROW/id"][:])
        sig_ids = _decode(h5["0/META/COL/id"][:])
        # These GEO GCTX files store metadata in the conventional
        # ROW=gene/COL=signature orientation, while the matrix is physically
        # laid out as signature x gene.  Check both dimensions explicitly.
        if matrix.shape != (len(sig_ids), len(gene_ids)):
            raise ValueError(f"GCTX shape/id mismatch for {dataset}: {matrix.shape}, genes={len(gene_ids)}, signatures={len(sig_ids)}")
        col_lookup = {str(gene_id): i for i, gene_id in enumerate(gene_ids)}
        bing_gene_ids = [str(x) for x in bing.pr_gene_id]
        missing_bing = [x for x in bing_gene_ids if x not in col_lookup]
        if missing_bing:
            raise ValueError(f"{dataset}: {len(missing_bing)} BING gene IDs absent from GCTX")
        bing_columns = np.asarray([col_lookup[x] for x in bing_gene_ids], dtype=np.int64)
        if len(np.unique(bing_columns)) != len(bing_columns):
            raise ValueError(f"{dataset}: duplicate BING columns detected")

        row_meta = _metadata_for_matrix(sig, sig_ids)
        valid = ~row_meta.pert_type.isin(NON_PERT_TYPES)
        valid_rows = np.flatnonzero(valid.to_numpy())
        row_meta = row_meta.iloc[valid_rows].reset_index(drop=True)
        out_dir = CANDIDATE_RESULTS
        tmp_dir = LINCS_CANDIDATE_INTERMEDIATE
        out_dir.mkdir(parents=True, exist_ok=True)
        tmp_dir.mkdir(parents=True, exist_ok=True)
        raw_path = tmp_dir / f"{dataset}_raw_WTCS.parquet"
        final_path = out_dir / f"{dataset}_LINCS_signature_WTCS.parquet"
        writer = None
        processed = 0
        try:
            for start in range(0, len(valid_rows), chunk_size):
                row_pos = valid_rows[start : start + chunk_size]
                values = np.asarray(matrix[row_pos, :][:, bing_columns], dtype=np.float32)
                scored = score_rows(values, queries, len(bing_columns))
                meta_chunk = row_meta.iloc[start : start + len(row_pos)].reset_index(drop=True)
                for query, scores in zip(queries, scored):
                    frame = meta_chunk.copy()
                    frame.insert(0, "dataset", dataset)
                    frame.insert(1, "query_id", query.query_id)
                    frame.insert(2, "contrast", query.contrast)
                    frame.insert(3, "canonical_query_hash", query.canonical_hash)
                    frame.insert(4, "target_size", query.actual_target_size)
                    frame.insert(5, "applicable_target_sizes", ",".join(map(str, query.requested_target_sizes)))
                    frame.insert(6, "n_up_mapped", len(query.up_indices))
                    frame.insert(7, "n_down_mapped", len(query.down_indices))
                    frame.insert(8, "n_up_unmapped", len(query.unmapped_up))
                    frame.insert(9, "n_down_unmapped", len(query.unmapped_down))
                    for col, values_out in scores.items():
                        frame[col] = values_out
                    frame["exploratory"] = query.exploratory
                    # Several requested sizes can intentionally map to the
                    # same actual gene set.  That is a stable, deduplicated
                    # query; the applicable sizes are audited separately.
                    frame["query_size_stable"] = True
                    frame["query_gene_count"] = len(query.up_indices) + len(query.down_indices)
                    table = pa.Table.from_pandas(frame, preserve_index=False)
                    if writer is None:
                        writer = pq.ParquetWriter(raw_path, table.schema, compression="zstd")
                    writer.write_table(table)
                    processed += len(frame)
                if start == 0 or (start // chunk_size) % 25 == 0:
                    print(f"[{dataset}] scored {start + len(row_pos):,}/{len(valid_rows):,} signatures; rows={processed:,}", flush=True)
        finally:
            if writer is not None:
                writer.close()
        denominators = _ncs_denominators(raw_path)
        denom_path = tmp_dir / f"{dataset}_NCS_denominators.csv"
        denominators.to_csv(denom_path, index=False)
        final_writer = None
        try:
            for batch in pq.ParquetFile(raw_path).iter_batches(batch_size=200_000):
                frame = _add_ncs(batch.to_pandas(), denominators)
                table = pa.Table.from_pandas(frame, preserve_index=False)
                if final_writer is None:
                    final_writer = pq.ParquetWriter(final_path, table.schema, compression="zstd")
                final_writer.write_table(table)
        finally:
            if final_writer is not None:
                final_writer.close()
    summary = {
        "dataset": dataset,
        "gctx": str(gctx_path),
        "n_matrix_rows": len(sig_ids),
        "n_matrix_genes": len(gene_ids),
        "n_bing_genes": len(bing_columns),
        "n_scored_signatures": len(valid_rows),
        "n_queries": len(queries),
        "n_output_rows": processed,
        "raw_intermediate": str(raw_path),
        "output": str(final_path),
    }
    return final_path, denom_path, summary


def _aggregate_frame(frame: pd.DataFrame, group_cols: list[str], prefix: str = "") -> pd.DataFrame:
    frame = frame.copy()
    frame["is_reversed"] = frame["reversal_WTCS"] > 0
    grouped = frame.groupby(group_cols, sort=False, dropna=False)
    result = grouped.agg(
        median_reversal_WTCS=("reversal_WTCS", "median"),
        median_reversal_NCS=("reversal_NCS", "median"),
        max_reversal_NCS=("reversal_NCS", "max"),
        fraction_reversed=("is_reversed", "mean"),
        n_signatures=("WTCS", "size"),
        query_size_stable=("query_size_stable", "all"),
        exploratory=("exploratory", "any"),
    ).reset_index()
    quantiles = (
        grouped["reversal_NCS"].quantile([0.25, 0.75]).unstack(level=-1)
        .rename(columns={0.25: "q25_reversal_NCS", 0.75: "q75_reversal_NCS"})
        .reset_index()
    )
    result = result.merge(quantiles, on=group_cols, how="left", validate="one_to_one")
    corr = grouped["weighted_reversal_correlation"].median().rename("median_weighted_reversal_correlation").reset_index()
    result = result.merge(corr, on=group_cols, how="left", validate="one_to_one")
    if prefix:
        rename = {column: f"{prefix}{column}" for column in result.columns if column not in group_cols}
        result = result.rename(columns=rename)
    return result


def aggregate_results(paths: list[Path]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frames = []
    for path in paths:
        keep = [
            "dataset", "query_id", "contrast", "canonical_query_hash", "target_size", "applicable_target_sizes",
            "pert_id", "pert_iname", "pert_type", "perturbation_key", "perturbation_class", "cell_id", "dose", "time",
            "ES_up", "ES_down", "WTCS", "reversal_WTCS", "NCS", "reversal_NCS", "weighted_reversal_correlation",
            "query_size_stable", "exploratory",
        ]
        frame = pd.read_parquet(path, columns=keep)
        frames.append(frame[keep])
    all_scores = pd.concat(frames, ignore_index=True)
    # Aggregate at normalized perturbation level, not BRD/shRNA ID level.
    # This lets cross-release matching recognize the same compound under
    # different BRD IDs and the same genetic target under shRNA versus xpr.
    group_cols = ["dataset", "canonical_query_hash", "query_id", "contrast", "perturbation_key", "perturbation_class", "pert_iname"]
    context_cols = group_cols + ["pert_type", "cell_id", "dose", "time"]
    context = _aggregate_frame(all_scores, context_cols)
    all_scores["is_reversed"] = all_scores["reversal_WTCS"] > 0
    grouped = all_scores.groupby(group_cols, sort=False, dropna=False)
    perturbation = grouped.agg(
        median_reversal_WTCS=("reversal_WTCS", "median"),
        median_reversal_NCS=("reversal_NCS", "median"),
        max_reversal_NCS=("reversal_NCS", "max"),
        fraction_reversed=("is_reversed", "mean"),
        n_signatures=("WTCS", "size"),
        n_cells=("cell_id", "nunique"),
        n_query_sizes=("target_size", "nunique"),
        query_size_stable=("query_size_stable", "all"),
        exploratory=("exploratory", "any"),
        applicable_target_sizes=("applicable_target_sizes", "first"),
        pert_type=("pert_type", lambda value: "|".join(sorted(set(value.astype(str))))),
    ).reset_index()
    quantiles = (
        grouped["reversal_NCS"].quantile([0.25, 0.75]).unstack(level=-1)
        .rename(columns={0.25: "q25_reversal_NCS", 0.75: "q75_reversal_NCS"})
        .reset_index()
    )
    perturbation = perturbation.merge(quantiles, on=group_cols, how="left", validate="one_to_one")
    best = (
        all_scores.assign(_best_score=all_scores["reversal_NCS"].fillna(-np.inf))
        .sort_values("_best_score", ascending=False)
        .drop_duplicates(group_cols)
    )
    best_context = best[group_cols + ["cell_id", "dose", "time"]].copy()
    best_context["best_context"] = (
        best_context["cell_id"].astype(str) + "|" + best_context["dose"].astype(str) + "|" + best_context["time"].astype(str)
    )
    perturbation = perturbation.merge(best_context[group_cols + ["best_context"]], on=group_cols, how="left", validate="one_to_one")
    recurrence = cross_release_recurrence(perturbation)
    return context, perturbation, recurrence


def _direction(row: pd.Series) -> str:
    value = row.get("median_reversal_NCS", np.nan)
    if not np.isfinite(value) or abs(float(value)) < 0.25:
        return "weak_or_neutral"
    return "reversal" if value > 0 else "mimic"


def cross_release_recurrence(perturbation: pd.DataFrame) -> pd.DataFrame:
    keys = ["canonical_query_hash", "perturbation_key", "perturbation_class", "pert_iname"]
    left = perturbation.loc[perturbation.dataset.eq("GSE92742")].copy()
    right = perturbation.loc[perturbation.dataset.eq("GSE70138")].copy()
    merged = left.merge(right, on=keys, how="outer", suffixes=("_GSE92742", "_GSE70138"), indicator=True)
    rows = []
    for _, row in merged.iterrows():
        available = row["_merge"] == "both"
        if not available:
            status = "replication_not_available"
        else:
            a = _direction(row.rename({"median_reversal_NCS_GSE92742": "median_reversal_NCS"}))
            b = _direction(row.rename({"median_reversal_NCS_GSE70138": "median_reversal_NCS"}))
            if {a, b} == {"reversal", "mimic"}:
                status = "replicated_discordant"
            elif a == b == "reversal" or a == b == "mimic":
                status = "replicated_concordant"
            else:
                status = "replication_available_but_weak"
        record = {key: row.get(key) for key in keys}
        record.update({
            "perturbation_class": row.get("perturbation_class"),
            "measured_in_GSE92742": bool(row["_merge"] in {"both", "left_only"}),
            "measured_in_GSE70138": bool(row["_merge"] in {"both", "right_only"}),
            "measured_in_other_dataset": available,
            "cross_phase_status": status,
            "median_reversal_NCS_GSE92742": row.get("median_reversal_NCS_GSE92742", np.nan),
            "median_reversal_NCS_GSE70138": row.get("median_reversal_NCS_GSE70138", np.nan),
            "fraction_reversed_GSE92742": row.get("fraction_reversed_GSE92742", np.nan),
            "fraction_reversed_GSE70138": row.get("fraction_reversed_GSE70138", np.nan),
            "n_cells_GSE92742": row.get("n_cells_GSE92742", 0),
            "n_cells_GSE70138": row.get("n_cells_GSE70138", 0),
            "pert_type_GSE92742": row.get("pert_type_GSE92742", ""),
            "pert_type_GSE70138": row.get("pert_type_GSE70138", ""),
        })
        rows.append(record)
    return pd.DataFrame(rows)


def positive_control_validation(perturbation: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, row in perturbation.iterrows():
        name = _norm_text(row.get("pert_iname", ""))
        gene = name.replace(" ", "")
        if name in CORE_MTOR_NAMES or gene in CORE_MTOR_GENES:
            panel = "core_mTOR"
        elif name in EXTENDED_PI3K_NAMES or gene in EXTENDED_PI3K_GENES:
            panel = "extended_PI3K_AKT"
        else:
            continue
        rows.append({
            "dataset": row["dataset"], "query_id": row["query_id"], "contrast": row["contrast"],
            "perturbation_key": row["perturbation_key"], "pert_iname": row["pert_iname"],
            "pert_type": row["pert_type"], "panel": panel,
            "median_reversal_WTCS": row["median_reversal_WTCS"],
            "median_reversal_NCS": row["median_reversal_NCS"],
            "max_reversal_NCS": row["max_reversal_NCS"], "n_cells": row["n_cells"],
            "sanity_check_only": True,
            "interpretation": "directional sanity check; not a hard pipeline acceptance criterion",
        })
    return pd.DataFrame(rows)


def self_test() -> None:
    # Small hand-verifiable rank universe with two query lists.
    assert math.isclose(weighted_ks_from_positions(np.array([0]), np.array([1.0]), 4), 1.0)
    assert math.isclose(weighted_ks_from_positions(np.array([3]), np.array([1.0]), 4), -1.0)
    assert wtcs(0.5, -0.5) == 0.5
    assert wtcs(0.5, 0.5) == 0.0
    assert wtcs(-0.5, -0.5) == 0.0
    assert wtcs(0.5, -0.5) == -wtcs(-0.5, 0.5)
    x = np.array([1.0, 2.0, -1.0])
    y = np.array([-1.0, -2.0, 1.0])
    assert weighted_pearson(x, y, np.abs(x)) < -0.99
    assert _normalise_perturbation_key("trt_cp", "sirolimus")[1] == "compound"
    assert _normalise_perturbation_key("trt_sh", "MTOR")[1] == "genetic"
    print("self_test: PASS")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--signature", default="results/signatures/GSE179044_cmap_query_signatures.csv")
    parser.add_argument("--datasets", nargs="+", choices=sorted(DEFAULT_MATRIX), default=sorted(DEFAULT_MATRIX))
    parser.add_argument("--chunk-size", type=int, default=256)
    parser.add_argument("--aggregate-only", action="store_true", help="Aggregate existing dataset Parquet outputs without rescoring GCTX")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    signature_path = ROOT / args.signature
    outputs = []
    manifests = []
    previous_manifest: dict = {}
    if args.aggregate_only:
        manifest_path = ROOT / "manifests/lincs_analysis.json"
        if manifest_path.exists():
            previous_manifest = json.loads(manifest_path.read_text())
        # Rebuild the lightweight query audit even when scoring is skipped.
        # This prevents an aggregate-only run from erasing the query manifest.
        query_gene = load_metadata(args.datasets[0])[1]
        queries = build_queries(signature_path, query_gene)
        for dataset in args.datasets:
            final_path = CANDIDATE_RESULTS / f"{dataset}_LINCS_signature_WTCS.parquet"
            if not final_path.exists():
                raise FileNotFoundError(f"Missing existing dataset output for --aggregate-only: {final_path}")
            outputs.append(final_path)
    else:
        # Both releases use the same BING feature order, but queries are built
        # from each release's local annotation and checked independently.
        query_gene = load_metadata(args.datasets[0])[1]
        queries = build_queries(signature_path, query_gene)
        if not queries:
            raise ValueError("No usable query signatures were generated")
        for dataset in args.datasets:
            final_path, denom_path, summary = run_dataset(dataset, queries, signature_path, args.chunk_size)
            outputs.append(final_path)
            manifests.append(summary | {"NCS_denominators": str(denom_path)})
    context, perturbation, recurrence = aggregate_results(outputs)
    out_dir = CANDIDATE_RESULTS
    out_dir.mkdir(parents=True, exist_ok=True)
    perturbation.to_csv(out_dir / "LINCS_perturbation_summary.csv.gz", index=False)
    context.to_csv(out_dir / "LINCS_context_summary.csv.gz", index=False)
    recurrence.to_csv(out_dir / "LINCS_cross_dataset_recurrence.csv", index=False)
    positive_control_validation(perturbation).to_csv(out_dir / "LINCS_positive_control_validation.csv", index=False)
    # Candidate ranking is deliberately descriptive only.  It does not create
    # Tier 1 candidates before human-state, toxicity, target and exposure filters.
    recurrence_keys = ["canonical_query_hash", "perturbation_key", "perturbation_class", "pert_iname"]
    ranking = perturbation.merge(
        recurrence[recurrence_keys + ["cross_phase_status", "measured_in_other_dataset"]],
        on=recurrence_keys,
        how="left",
        validate="many_to_one",
    )
    ranking["tier"] = np.where(ranking["median_reversal_NCS"].fillna(0) > 0.25, "exploratory_reversal", "weak_or_neutral")
    ranking["tier1_eligible"] = False
    ranking.to_csv(out_dir / "LINCS_candidate_ranking.csv", index=False)
    manifest_queries = [
        {"query_id": q.query_id, "contrast": q.contrast, "canonical_query_hash": q.canonical_hash,
         "requested_target_sizes": q.requested_target_sizes, "n_up": len(q.up_indices), "n_down": len(q.down_indices),
         "unmapped_up": q.unmapped_up, "unmapped_down": q.unmapped_down, "exploratory": q.exploratory}
        for q in queries
    ]
    manifest_summaries = previous_manifest.get("matrix_summaries", []) if args.aggregate_only else manifests
    if args.aggregate_only and not manifest_summaries:
        manifest_summaries = []
        for dataset, output in zip(args.datasets, outputs):
            cfg = DEFAULT_MATRIX[dataset]
            sig, gene, _ = load_metadata(dataset)
            with h5py.File(ROOT / cfg["gctx"], "r") as h5:
                matrix_shape = list(h5["0/DATA/0/matrix"].shape)
            manifest_summaries.append({
                "dataset": dataset,
                "gctx": str(ROOT / cfg["gctx"]),
                "n_matrix_rows": matrix_shape[0],
                "n_matrix_genes": matrix_shape[1],
                "n_bing_genes": int(gene["pr_is_bing"].astype(int).sum()),
                "n_scored_signatures": int((~sig["pert_type"].isin(NON_PERT_TYPES)).sum()),
                "n_queries": len(queries),
                "n_output_rows": int(pq.ParquetFile(output).metadata.num_rows),
                "output": str(output),
            })
    write_json(ROOT / "manifests/lincs_analysis.json", {
        "status": "completed",
        "datasets": args.datasets,
        "analysis_name": "local_LINCS_CMap_residual_programs",
        "recurrence_label": "cross-phase/cross-release recurrence; not independent biological replication",
        "signature_file": str(signature_path),
        "n_queries": len(manifest_queries),
        "queries": manifest_queries,
        "matrix_summaries": manifest_summaries,
        "methods": {
            "feature_space": "complete BING (10,174 genes)",
            "es": "weighted KS on full BING rank space, hit weight abs(Level5 score), exponent 1",
            "wtcs": "(ES_up - ES_down)/2 only when ES_up*ES_down < 0; otherwise zero",
            "ncs": "published mean-scaling within dataset/query/cell/pert_type, negative side sign preserved",
            "weighted_correlation": "negative weighted Pearson with weight abs(disease signed_score)",
            "tau": "not calculated; no fixed Touchstone query compendium",
        },
        "hard_acceptance": ["synthetic ES", "up/down swap", "BING mapping", "same-sign WTCS zero", "NCS denominator/sign", "full BING coverage", "duplicate query-size dedup"],
        "biological_sanity_checks": ["core mTOR/rapamycin panel", "extended PI3K/AKT panel"],
        "tier1_generated": 0,
        "note": "Cross-release absence is not a failed replication; discordance is recorded only when both releases contain the perturbation.",
    })
    print(json.dumps({
        "status": "completed", "n_queries": len(manifest_queries), "n_perturbation_rows": len(perturbation),
        "n_recurrence_rows": len(recurrence), "tier1_generated": 0,
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
