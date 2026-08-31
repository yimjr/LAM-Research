#!/usr/bin/env python3
"""Step 10: independent patient-aware DE and program annotation per state."""

from __future__ import annotations

import argparse
import gc
import os
from pathlib import Path

os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/lam-state-mpl")

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from statsmodels.stats.multitest import multipletests

from state_modeling_utils import PROJECT_ROOT, load_config, resolve_shared, validate_integer_counts, write_json


def matrix_row_sum(matrix, positions: np.ndarray) -> np.ndarray:
    values = matrix[positions].sum(axis=0)
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel()


def aggregate_counts(matrix, positions: np.ndarray) -> np.ndarray:
    return matrix_row_sum(matrix, positions).astype(np.int64)


def mean_expression(matrix, positions: np.ndarray, gene_positions: list[int]) -> np.ndarray:
    if len(positions) == 0 or len(gene_positions) == 0:
        return np.zeros(len(gene_positions), dtype=float)
    values = matrix[positions][:, gene_positions].mean(axis=0)
    if sparse.issparse(values):
        values = values.toarray()
    return np.asarray(values).ravel().astype(float)


def program_definitions(config: dict) -> tuple[list[dict], str]:
    path = resolve_shared(config, config["annotation_files"]["shared"]["known_programs"])
    if path is None:
        return [], "not_available"
    payload = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    programs = []
    for item in payload.get("programs", []):
        genes = ["VEGFD" if str(g).upper() == "FIGF" else str(g).upper() for g in (item.get("genes") or [])]
        if genes:
            programs.append({
                "program_name": str(item.get("program_name", "unnamed")),
                "genes": list(dict.fromkeys(genes)),
                "evidence_scope": item.get("evidence_scope", ""),
                "evidence_level": item.get("evidence_level", ""),
                "source_study": item.get("source_study", ""),
            })
    return programs, str(path)


def load_high_prepared(config: dict, consensus: ad.AnnData) -> ad.AnnData:
    path = PROJECT_ROOT / config["outputs"]["prepared_h5ad"]
    prepared = ad.read_h5ad(path, backed="r")
    if "counts" not in prepared.layers:
        raise RuntimeError("BLOCKED_INPUT: prepared AnnData has no raw counts layer")
    # The state artifact and prepared artifact share the canonical AnnData
    # index.  Keep an analysis_cell_id fallback for older upstream artifacts.
    target_ids = consensus.obs.index.astype(str)
    source_ids = prepared.obs.index.astype(str)
    positions = pd.Series(np.arange(len(source_ids)), index=source_ids).reindex(target_ids)
    if positions.isna().any() and "analysis_cell_id" in prepared.obs and "analysis_cell_id" in consensus.obs:
        lookup = pd.Series(np.arange(len(prepared)), index=prepared.obs["analysis_cell_id"].astype(str))
        positions = lookup.reindex(consensus.obs["analysis_cell_id"].astype(str).to_numpy())
    if positions.isna().any():
        raise RuntimeError("BLOCKED_INPUT: consensus cells cannot be mapped to prepared AnnData")
    high = prepared[positions.astype(int).to_numpy()].to_memory()
    high.obs_names = consensus.obs_names.copy()
    high.obs = consensus.obs.copy()
    audit = validate_integer_counts(high)
    if not audit["valid"]:
        raise RuntimeError(f"BLOCKED_INPUT: prepared raw counts failed integer-value validation: {audit}")
    high.uns["step10_counts_audit"] = audit
    return high


def build_pseudobulk(high: ad.AnnData, state: str, min_state: int, min_rest: int) -> tuple[pd.DataFrame, list[dict], dict[str, dict[str, np.ndarray]]]:
    labels = high.obs["consensus_state"].astype(str).to_numpy()
    patients = high.obs["patient_id"].astype(str).to_numpy()
    state_group = f"State_{state}"
    matrix = high.layers["counts"]
    rows = []
    metadata = []
    patient_vectors: dict[str, dict[str, np.ndarray]] = {}
    for patient in sorted(pd.unique(patients)):
        patient_pos = np.flatnonzero(patients == patient)
        state_pos = patient_pos[labels[patient_pos] == state]
        rest_pos = patient_pos[labels[patient_pos] != state]
        if len(state_pos) < min_state or len(rest_pos) < min_rest:
            continue
        state_counts = aggregate_counts(matrix, state_pos)
        rest_counts = aggregate_counts(matrix, rest_pos)
        for group, positions, counts in [(state_group, state_pos, state_counts), ("Rest_of_LAM", rest_pos, rest_counts)]:
            sample_id = f"{patient}__{group}"
            rows.append(np.asarray(counts, dtype=np.int64))
            metadata.append({
                "state_id": state,
                "sample_id": sample_id,
                "patient_id": patient,
                "group": group,
                "cells": int(len(positions)),
                "total_counts": int(counts.sum()),
            })
        patient_vectors[patient] = {state_group: state_counts, "Rest_of_LAM": rest_counts}
    if not rows:
        return pd.DataFrame(), metadata, patient_vectors
    counts = pd.DataFrame(np.vstack(rows), columns=high.var_names.astype(str))
    meta = pd.DataFrame(metadata)
    return pd.concat([meta.reset_index(drop=True), counts], axis=1), metadata, patient_vectors


def patient_direction(patient_vectors: dict[str, dict[str, np.ndarray]], state_group: str, total_genes: int) -> np.ndarray:
    effects = []
    for values in patient_vectors.values():
        state_counts = values[state_group].astype(float)
        rest_counts = values["Rest_of_LAM"].astype(float)
        state_cpm = state_counts / max(1.0, state_counts.sum()) * 1e6
        rest_cpm = rest_counts / max(1.0, rest_counts.sum()) * 1e6
        effects.append(np.log2((state_cpm + 0.5) / (rest_cpm + 0.5)))
    if not effects:
        return np.full(total_genes, np.nan)
    return np.vstack(effects)


def run_one_state(
    state: str,
    pb: pd.DataFrame,
    patient_vectors: dict[str, dict[str, np.ndarray]],
    gene_names: list[str],
    config: dict,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    state_group = f"State_{state}"
    n_patients = int(pb["patient_id"].nunique()) if not pb.empty else 0
    if pb.empty:
        return pd.DataFrame(), pd.DataFrame(), {"state_id": state, "status": "descriptive_only", "n_patients": 0, "reason": "no eligible patient has both groups"}
    count_data = pb[gene_names].copy()
    detected = (count_data > 0).sum(axis=0)
    keep = (count_data.sum(axis=0) >= int(config["step10"]["min_gene_count"])) & (detected >= int(config["step10"]["min_gene_samples"]))
    count_data = count_data.loc[:, keep]
    directions = patient_direction(patient_vectors, state_group, len(gene_names))
    direction_map = dict(zip(gene_names, directions.T))
    if n_patients < int(config["step10"]["min_de_patients"]) or count_data.shape[1] == 0:
        marker = pd.DataFrame({"state_id": [state], "status": ["descriptive_only"], "n_patients": [n_patients]})
        return pd.DataFrame(), marker, {
            "state_id": state,
            "status": "descriptive_only",
            "n_patients": n_patients,
            "n_genes_tested": int(count_data.shape[1]),
            "reason": "fewer than min_de_patients or no genes passed count filter",
        }
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    metadata = pb[["patient_id", "group"]].copy()
    metadata.index = pb["sample_id"].astype(str)
    counts = count_data.copy()
    counts.index = metadata.index
    counts = counts.astype(np.int64)
    status = {
        "state_id": state,
        "status": "formal_de_attempted",
        "n_patients": n_patients,
        "n_genes_tested": int(counts.shape[1]),
        "design": "~ patient_id + group",
        "contrast": f"group: {state_group} vs Rest_of_LAM",
    }
    try:
        dds = DeseqDataSet(
            counts=counts,
            metadata=metadata,
            design_factors=["patient_id", "group"],
            n_cpus=1,
            min_replicates=max(1, min(3, n_patients)),
            quiet=True,
        )
        dds.deseq2()
        stats = DeseqStats(dds, contrast=["group", state_group, "Rest_of_LAM"], alpha=float(config["step10"]["fdr_alpha"]), n_cpus=1, quiet=True)
        stats.summary()
        result = stats.results_df.reset_index().rename(columns={"index": "gene"})
    except Exception as exc:
        status.update({"status": "formal_de_failed", "error": f"{type(exc).__name__}: {exc}"})
        return pd.DataFrame(), pd.DataFrame([status]), status
    result["state_id"] = state
    result["contrast"] = f"{state_group} vs Rest_of_LAM"
    result["n_patients"] = n_patients
    effect_sign = np.sign(pd.to_numeric(result["log2FoldChange"], errors="coerce").to_numpy())
    concordance = []
    for gene, sign in zip(result["gene"].astype(str), effect_sign):
        values = np.asarray(direction_map.get(gene, np.asarray([np.nan])), dtype=float)
        values = values[np.isfinite(values)]
        concordance.append(float(np.mean(np.sign(values) == sign)) if len(values) and sign != 0 else np.nan)
    result["patient_direction_concordance"] = concordance
    result["direction"] = np.where(result["log2FoldChange"] > 0, "up", np.where(result["log2FoldChange"] < 0, "down", "flat"))
    result = result[["state_id", "gene", "contrast", "n_patients", "log2FoldChange", "pvalue", "padj", "direction", "patient_direction_concordance"]]
    marker = result[result["padj"].le(float(config["step10"]["fdr_alpha"])) & result["log2FoldChange"].abs().ge(float(config["step10"]["min_abs_log2fc"]))].copy()
    marker["status"] = "marker_candidate"
    return result, marker, status


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/state_modeling.yaml")
    args = parser.parse_args()
    config = load_config(args.config)
    out = PROJECT_ROOT / config["outputs"]["step10_dir"]
    out.mkdir(parents=True, exist_ok=True)
    consensus = ad.read_h5ad(PROJECT_ROOT / config["outputs"]["consensus_h5ad"])
    high = load_high_prepared(config, consensus)
    gene_names = high.var_names.astype(str).tolist()
    states = sorted(high.obs["consensus_state"].astype(str).unique(), key=lambda value: (len(value), value))

    pb_path = out / "state_pseudobulk_counts.csv"
    de_frames = []
    marker_frames = []
    program_rows = []
    statuses = []
    first_pb = True
    programs, program_source = program_definitions(config)
    normalized_matrix = high.X
    labels = high.obs["consensus_state"].astype(str).to_numpy()
    for state in states:
        pb, metadata, patient_vectors = build_pseudobulk(
            high,
            state,
            int(config["step10"]["min_state_cells_per_patient"]),
            int(config["step10"]["min_rest_cells_per_patient"]),
        )
        if not pb.empty:
            pb.to_csv(pb_path, mode="w" if first_pb else "a", header=first_pb, index=False)
            first_pb = False
        de, markers, status = run_one_state(state, pb, patient_vectors, gene_names, config)
        statuses.append(status)
        if not de.empty:
            de_frames.append(de)
        if not markers.empty:
            marker_frames.append(markers)
        state_positions = np.flatnonzero(labels == state)
        rest_positions = np.flatnonzero(labels != state)
        for program in programs:
            positions = [gene_names.index(gene) for gene in program["genes"] if gene in gene_names]
            state_mean = float(mean_expression(normalized_matrix, state_positions, positions).mean()) if positions else np.nan
            rest_mean = float(mean_expression(normalized_matrix, rest_positions, positions).mean()) if positions and len(rest_positions) else np.nan
            program_rows.append({
                "state_id": state,
                "program_name": program["program_name"],
                "genes_present": int(len(positions)),
                "genes_requested": int(len(program["genes"])),
                "state_mean_expression": state_mean,
                "rest_mean_expression": rest_mean,
                "delta_state_minus_rest": state_mean - rest_mean if np.isfinite(state_mean) and np.isfinite(rest_mean) else np.nan,
                "evidence_scope": program["evidence_scope"],
                "evidence_level": program["evidence_level"],
                "source_study": program["source_study"],
            })
        del pb, patient_vectors
        gc.collect()

    de_all = pd.concat(de_frames, ignore_index=True) if de_frames else pd.DataFrame(columns=["state_id", "gene", "contrast", "n_patients", "log2FoldChange", "pvalue", "padj", "direction", "patient_direction_concordance"])
    if not de_all.empty:
        pvalues = pd.to_numeric(de_all["pvalue"], errors="coerce")
        valid = pvalues.notna().to_numpy()
        de_all["global_padj"] = np.nan
        if valid.any():
            de_all.loc[valid, "global_padj"] = multipletests(pvalues[valid], method="fdr_bh")[1]
    de_all.to_csv(out / "state_de_results.csv", index=False)
    (pd.concat(marker_frames, ignore_index=True) if marker_frames else pd.DataFrame()).to_csv(out / "state_markers.csv", index=False)
    pd.DataFrame(program_rows).to_csv(out / "state_program_scores.csv", index=False)
    pathway_rows = [{"state_id": state, "status": "not_available", "source": "no local state-specific pathway result in current scope"} for state in states]
    regulon_rows = [{"state_id": state, "status": "not_available", "source": "no local state-specific regulon result in current scope"} for state in states]
    pd.DataFrame(pathway_rows).to_csv(out / "state_pathway_enrichment.csv", index=False)
    pd.DataFrame(regulon_rows).to_csv(out / "state_regulon_summary.csv", index=False)
    write_json(out / "step10_manifest.json", {
        "n_states": len(states),
        "state_independent_models": True,
        "design": "~ patient_id + group",
        "contrast_template": "State_k vs Rest_of_LAM",
        "min_de_patients": int(config["step10"]["min_de_patients"]),
        "program_source": program_source,
        "scvi_training_called": False,
        "state_status": statuses,
    })
    report_lines = [
        "# Stage 10 biology annotation",
        "",
        f"- States analyzed independently: {len(states)}",
        "- Each eligible state uses patient × group pseudobulk and `~ patient_id + group`.",
        "- Each state has its own DE FDR; `global_padj` is supplemental and does not replace state-level FDR.",
        "- No all-state multi-class or pooled binary DE model was constructed.",
        "- Pathway/regulon outputs are explicit `not_available` placeholders unless a local state-specific result exists.",
    ]
    (PROJECT_ROOT / "reports/stage10_biology_annotation.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    print(f"Step 10 complete: {len(states)} state-specific models attempted")


if __name__ == "__main__":
    main()
