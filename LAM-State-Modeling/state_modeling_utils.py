"""Shared, deliberately small utilities for the numbered LAM state scripts."""

from __future__ import annotations

import json
import os
import re
from pathlib import Path
from typing import Any, Iterable

import anndata as ad
import numpy as np
import pandas as pd
import yaml
from scipy import sparse


# The project venv intentionally lives outside the repository.  Point Numba's
# optional cache to a writable Linux path instead of the read-only default
# home cache used by this managed WSL session.
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/lam-state-numba-cache")

PROJECT_ROOT = Path(__file__).resolve().parent


def load_config(path: str | Path = "config/state_modeling.yaml") -> dict[str, Any]:
    path = Path(path)
    if not path.is_absolute():
        path = PROJECT_ROOT / path
    with path.open(encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def ensure_dir(path: str | Path) -> Path:
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def resolve_path(candidates: Iterable[str], roots: Iterable[str | Path] | None = None) -> Path | None:
    """Resolve the first existing candidate without copying or mutating inputs."""
    roots = list(roots or ["."])
    for root in roots:
        root_path = Path(root)
        if not root_path.is_absolute():
            root_path = PROJECT_ROOT / root_path
        for candidate in candidates:
            candidate_path = Path(candidate)
            path = candidate_path if candidate_path.is_absolute() else root_path / candidate_path
            if path.exists():
                return path.resolve()
    return None


def resolve_dataset_h5ad(config: dict[str, Any], dataset: str) -> Path | None:
    spec = config["datasets"][dataset]
    return resolve_path(spec.get("h5ad_candidates", []), config.get("input_roots", []))


def resolve_shared(config: dict[str, Any], relative_path: str) -> Path | None:
    return resolve_path([relative_path], config.get("input_roots", []))


def annotation_directory(config: dict[str, Any], dataset: str) -> Path | None:
    rel = config["datasets"][dataset].get("annotation_dir")
    if not rel:
        return None
    return resolve_path([rel], config.get("input_roots", []))


def json_default(value: Any) -> Any:
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value)!r}")


def write_json(path: str | Path, payload: Any) -> None:
    path = ensure_dir(Path(path).parent) / Path(path).name
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=json_default) + "\n", encoding="utf-8")


def as_bool(value: Any) -> bool:
    if pd.isna(value):
        return False
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(value)
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def bool_series(values: pd.Series, index: pd.Index | None = None) -> pd.Series:
    result = values.map(as_bool)
    if index is not None:
        result.index = index
    return result.astype(bool)


def safe_column(value: Any) -> str:
    text = re.sub(r"[^0-9A-Za-z_]+", "_", str(value).strip()).strip("_")
    return text or "unnamed"


def matrix_data(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.data)
    return np.asarray(matrix)


def validate_integer_counts(adata: ad.AnnData) -> dict[str, Any]:
    """Validate numeric integer-valued counts; integer dtype is not required."""
    source = "layers[counts]" if "counts" in adata.layers else "X"
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    values = matrix_data(matrix)
    finite = bool(np.isfinite(values).all())
    nonnegative = bool((values >= 0).all()) if values.size else True
    integer_valued = bool(np.isclose(values, np.round(values), rtol=0, atol=1e-6).all()) if values.size else True
    valid = finite and nonnegative and integer_valued
    return {
        "source": source,
        "dtype": str(getattr(matrix, "dtype", "unknown")),
        "n_values_checked": int(values.size),
        "finite": finite,
        "nonnegative": nonnegative,
        "integer_valued": integer_valued,
        "valid": valid,
    }


def ensure_counts_layer(adata: ad.AnnData) -> tuple[dict[str, Any], bool]:
    """Ensure counts exists, accepting integer-valued float X with an audit flag."""
    audit = validate_integer_counts(adata)
    copied_from_x = "counts" not in adata.layers
    if not audit["valid"]:
        raise ValueError(f"Invalid raw counts: {audit}")
    if copied_from_x:
        adata.layers["counts"] = adata.X.copy()
        audit["warning"] = "layers[counts] was absent; integer-valued X was copied as counts"
    return audit, copied_from_x


def _aggregate_matrix(matrix: Any, groups: list[list[int]]) -> Any:
    if len(groups) == matrix.shape[1] and all(len(group) == 1 for group in groups):
        return matrix.copy()
    if sparse.issparse(matrix):
        # A column-by-column hstack is prohibitively slow for 30k–60k gene
        # matrices.  A sparse membership matrix performs the same summation
        # in one sparse matrix multiplication.
        n_old = matrix.shape[1]
        group_ids = np.empty(n_old, dtype=np.int64)
        for new_index, group in enumerate(groups):
            group_ids[np.asarray(group, dtype=np.int64)] = new_index
        membership = sparse.csr_matrix(
            (np.ones(n_old, dtype=np.float32), (np.arange(n_old), group_ids)),
            shape=(n_old, len(groups)),
        )
        return (matrix.tocsr() @ membership).tocsr()
    columns = []
    for group in groups:
        columns.append(np.asarray(matrix[:, group]).sum(axis=1, keepdims=True))
    return np.column_stack(columns)


def canonicalize_gene_aliases(adata: ad.AnnData, aliases: dict[str, str] | None = None) -> tuple[ad.AnnData, dict[str, Any]]:
    """Canonicalize gene symbols and sum duplicate columns, including FIGF→VEGFD."""
    aliases = {str(k).upper(): str(v).upper() for k, v in (aliases or {"FIGF": "VEGFD"}).items()}
    original_names = adata.var_names.astype(str).tolist()
    symbol_col = "gene_symbol" if "gene_symbol" in adata.var else None
    symbols = adata.var[symbol_col].astype(str).tolist() if symbol_col else original_names
    canonical = [aliases.get(symbol.strip().upper(), symbol.strip()) for symbol in symbols]
    canonical = [name if name else original_names[i] for i, name in enumerate(canonical)]

    groups: dict[str, list[int]] = {}
    for i, name in enumerate(canonical):
        groups.setdefault(name, []).append(i)
    ordered_names = list(groups)
    group_indices = list(groups.values())
    alias_groups = {name: [original_names[i] for i in indices] for name, indices in groups.items() if len(indices) > 1 or any(original_names[i].upper() in aliases for i in indices)}

    old_x = adata.X
    old_layers = {name: value for name, value in adata.layers.items()}
    new_x = _aggregate_matrix(old_x, group_indices)
    new_layers = {
        layer_name: _aggregate_matrix(matrix, group_indices)
        for layer_name, matrix in old_layers.items()
        if getattr(matrix, "shape", None) == (adata.n_obs, len(original_names))
    }
    var = adata.var.copy().iloc[[indices[0] for indices in group_indices]].copy()
    var.index = pd.Index(ordered_names, dtype=str)
    var["gene_symbol"] = ordered_names
    var["gene_symbol_upper"] = [name.upper() for name in ordered_names]
    var["state_model_gene_alias"] = [";".join(alias_groups.get(name, [])) if name in alias_groups else "" for name in ordered_names]
    # Change var and X together.  Assigning X before var would fail AnnData's
    # shape validation when an alias merge reduces the number of genes.
    obsm = {key: value.copy() for key, value in adata.obsm.items()}
    uns = dict(adata.uns)
    adata = ad.AnnData(X=new_x, obs=adata.obs.copy(), var=var, layers=new_layers, uns=uns, obsm=obsm)
    # The caller needs the replacement object when aliases were merged.  Keep
    # the return contract backward-compatible by returning it in the audit.
    audit = {
        "aliases": aliases,
        "n_genes_before": len(original_names),
        "n_genes_after": len(ordered_names),
        "merged_groups": alias_groups,
        "requires_recomputed_log_x": any(
            len(group) > 1 and target in {str(name).upper() for name in group}
            for target, group in alias_groups.items()
        ),
    }
    adata.uns["state_model_gene_alias_audit"] = audit
    return adata, audit


def row_sums(matrix: Any) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(float)


def row_nonzero(matrix: Any) -> np.ndarray:
    if sparse.issparse(matrix):
        return np.asarray(matrix.getnnz(axis=1)).ravel().astype(int)
    return np.asarray((np.asarray(matrix) != 0).sum(axis=1)).ravel().astype(int)


def calculate_count_qc(adata: ad.AnnData) -> None:
    counts = adata.layers["counts"]
    total = row_sums(counts)
    n_genes = row_nonzero(counts)
    symbols = pd.Series(
        adata.var.get("gene_symbol_upper", pd.Series(adata.var_names, index=adata.var_names)),
        index=adata.var_names,
    ).astype(str).str.upper()
    mt_mask = symbols.str.match(r"^MT-").to_numpy(dtype=bool)
    mt_total = row_sums(counts[:, mt_mask]) if mt_mask.any() else np.zeros(adata.n_obs, dtype=float)
    adata.obs["total_counts"] = total
    adata.obs["n_genes_by_counts"] = n_genes
    adata.obs["total_counts_mt"] = mt_total
    adata.obs["pct_counts_mt"] = np.divide(mt_total * 100.0, total, out=np.zeros_like(mt_total), where=total > 0)
    adata.var["mt"] = mt_mask


def assay_family(assay: Any) -> str:
    text = str(assay).lower()
    return "snrna" if "sn" in text else "scrna"


def threshold_for_assay(assay: Any, qc_config: dict[str, Any]) -> float:
    values = qc_config.get("mt_pct_by_assay", {})
    return float(values.get(assay_family(assay), values.get("default", 20.0)))


def registry_records(registry: dict[str, Any]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for item in registry.get("donors", []):
        item = dict(item)
        item.setdefault("sample_id", item.get("donor_id"))
        item.setdefault("patient_id", item.get("donor_id"))
        item.setdefault("specimen_id", item.get("sample_id"))
        records.append(item)
    for accession, rows in registry.get("sample_mappings", {}).items():
        for row in rows or []:
            item = dict(row)
            item["accession"] = accession
            records.append(item)
    return records


def _key(value: Any) -> str:
    return "" if value is None or (isinstance(value, float) and np.isnan(value)) else str(value).strip()


def build_registry_lookup(records: list[dict[str, Any]], accession: str) -> dict[str, dict[str, Any]]:
    lookup: dict[str, dict[str, Any]] = {}
    for record in records:
        record_accession = _key(record.get("accession"))
        if record_accession and accession not in record_accession and record_accession not in accession:
            continue
        for field in ["sample_id", "source_sample", "specimen_id", "donor_id", "gsm"]:
            key = _key(record.get(field))
            if key:
                lookup.setdefault(key, record)
    return lookup


def apply_registry_mapping(adata: ad.AnnData, accession: str, registry: dict[str, Any]) -> dict[str, Any]:
    lookup = build_registry_lookup(registry_records(registry), accession)
    obs = adata.obs.copy()
    existing = {field: obs[field].astype(str).tolist() if field in obs else [""] * len(obs) for field in ["patient_id", "donor_id"]}
    mapped_rows: list[dict[str, Any] | None] = []
    statuses: list[str] = []
    conflicts: list[dict[str, Any]] = []
    for pos, (cell_id, row) in enumerate(obs.iterrows()):
        candidates = []
        for field in ["sample_id", "specimen_id", "source_sample", "donor_id"]:
            if field in row:
                value = _key(row[field])
                if value:
                    candidates.append(value)
        cell_text = str(cell_id)
        if ":" in cell_text:
            candidates.append(cell_text.split(":", 1)[0])
        record = next((lookup.get(candidate) for candidate in candidates if lookup.get(candidate)), None)
        mapped_rows.append(record)
        if record is None:
            statuses.append("unresolved")
            continue
        statuses.append("registry")
        for field in ["patient_id", "donor_id"]:
            old = existing[field][pos]
            new = _key(record.get(field))
            if old and old not in {"nan", "None", "unknown"} and new and old != new:
                conflicts.append({"cell_id": str(cell_id), "field": field, "upstream": old, "registry": new})

    for field in ["patient_id", "donor_id", "specimen_id", "assay", "tissue", "condition", "independence_group"]:
        values = []
        for record, old in zip(mapped_rows, obs.get(field, pd.Series([""] * len(obs), index=obs.index)).tolist()):
            values.append(record.get(field, old) if record else old)
        obs[field] = pd.Series(values, index=obs.index).astype(str)
    obs["mapping_status"] = statuses
    obs["mapping_source"] = np.where(np.asarray(statuses) == "registry", "donor_registry.yaml", "unresolved")
    adata.obs = obs
    return {
        "accession": accession,
        "n_cells": int(len(obs)),
        "n_registry_mapped": int(sum(status == "registry" for status in statuses)),
        "n_unresolved": int(sum(status == "unresolved" for status in statuses)),
        "n_conflicts": int(len(conflicts)),
        "conflicts": conflicts[:50],
    }


def join_cell_table(adata: ad.AnnData, table: pd.DataFrame, prefix: str, table_name: str) -> int:
    if "cell_id" not in table.columns:
        raise ValueError(f"{table_name} has no cell_id column")
    table = table.copy()
    table["cell_id"] = table["cell_id"].astype(str)
    if table["cell_id"].duplicated().any():
        value_cols = [col for col in table.columns if col != "cell_id"]
        aggregations: dict[str, Any] = {}
        for col in value_cols:
            aggregations[col] = "mean" if pd.api.types.is_numeric_dtype(table[col]) else "first"
        table = table.groupby("cell_id", as_index=False).agg(aggregations)
    table = table.set_index("cell_id")
    source_ids = pd.Index(adata.obs["source_cell_id"].astype(str))
    matched = int(source_ids.isin(table.index).sum())
    aligned = table.reindex(source_ids)
    aligned.index = adata.obs.index
    for col in aligned.columns:
        name = f"{prefix}_{safe_column(col)}"
        values = aligned[col]
        if pd.api.types.is_bool_dtype(values):
            adata.obs[name] = values.fillna(False).astype(bool).to_numpy()
        elif pd.api.types.is_numeric_dtype(values):
            adata.obs[name] = pd.to_numeric(values, errors="coerce").to_numpy()
        else:
            adata.obs[name] = values.fillna("").astype(str).to_numpy()
    return matched


def join_key_table(
    adata: ad.AnnData,
    table: pd.DataFrame,
    obs_key: str,
    table_key: str,
    prefix: str,
    table_name: str,
    exclude_columns: set[str] | None = None,
) -> int:
    """Attach a donor/sample-level upstream table without using it for modeling."""
    if obs_key not in adata.obs or table_key not in table.columns:
        raise ValueError(f"{table_name} cannot join {obs_key} <- {table_key}")
    table = table.copy()
    table[table_key] = table[table_key].astype(str)
    exclude_columns = set(exclude_columns or ())
    value_cols = [col for col in table.columns if col != table_key and col not in exclude_columns]
    if table[table_key].duplicated().any():
        aggregations: dict[str, Any] = {}
        for col in value_cols:
            if pd.api.types.is_bool_dtype(table[col]):
                aggregations[col] = "max"
            elif pd.api.types.is_numeric_dtype(table[col]):
                aggregations[col] = "mean"
            else:
                aggregations[col] = "first"
        table = table.groupby(table_key, as_index=False).agg(aggregations)
    table = table.set_index(table_key)
    source_keys = pd.Index(adata.obs[obs_key].astype(str))
    matched = int(source_keys.isin(table.index).sum())
    aligned = table.reindex(source_keys)
    aligned.index = adata.obs.index
    for col in value_cols:
        if col not in aligned:
            continue
        name = f"{prefix}_{safe_column(col)}"
        values = aligned[col]
        if pd.api.types.is_bool_dtype(values):
            adata.obs[name] = values.fillna(False).astype(bool).to_numpy()
        elif pd.api.types.is_numeric_dtype(values):
            adata.obs[name] = pd.to_numeric(values, errors="coerce").to_numpy()
        else:
            adata.obs[name] = values.fillna("").astype(str).to_numpy()
    return matched


def attach_candidate_annotation(adata: ad.AnnData, table_path: Path, dataset: str) -> dict[str, Any]:
    table = pd.read_csv(table_path)
    required = {"cell_id", "pool_high_confidence", "pool_broad_lam_like", "pool_unrestricted_lam"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"{table_path} missing required candidate columns: {missing}")
    table["cell_id"] = table["cell_id"].astype(str)
    table = table.drop_duplicates("cell_id", keep="last").set_index("cell_id")
    source_ids = pd.Index(adata.obs["source_cell_id"].astype(str))
    aligned = table.reindex(source_ids)
    aligned.index = adata.obs.index
    pool_names = ["high_confidence", "broad_lam_like", "unrestricted_lam"]
    for pool in pool_names:
        adata.obs[f"upstream_pool_{pool}"] = aligned[f"pool_{pool}"].map(as_bool).fillna(False).astype(bool).to_numpy()
    for col in table.columns:
        if col.startswith("pool_"):
            continue
        name = f"upstream_{dataset}_candidate_{safe_column(col)}"
        values = aligned[col]
        if pd.api.types.is_numeric_dtype(table[col]):
            adata.obs[name] = pd.to_numeric(values, errors="coerce").to_numpy()
        else:
            adata.obs[name] = values.fillna("").astype(str).to_numpy()
    matched = np.asarray(source_ids.isin(table.index), dtype=bool)
    adata.obs["upstream_candidate_annotation_matched"] = matched
    adata.obs["lam_candidate"] = adata.obs["upstream_pool_high_confidence"].astype(bool)
    adata.obs["boundary"] = adata.obs["upstream_pool_broad_lam_like"].astype(bool) & ~adata.obs["lam_candidate"].astype(bool)
    unrestricted_only = adata.obs["upstream_pool_unrestricted_lam"].astype(bool) & ~adata.obs["upstream_pool_broad_lam_like"].astype(bool)
    adata.obs["analysis_role"] = np.select(
        [adata.obs["lam_candidate"].to_numpy(), adata.obs["boundary"].to_numpy(), unrestricted_only.to_numpy()],
        ["primary_candidate", "boundary", "unrestricted_audit_only"],
        default="excluded_context",
    )
    return {
        "path": str(table_path),
        "n_rows": int(len(table)),
        "n_matched": int(matched.sum()),
        "n_unmatched": int((~matched).sum()),
        "pool_counts": {pool: int(adata.obs[f"upstream_pool_{pool}"].sum()) for pool in pool_names},
    }


def attach_upstream_cell_tables(adata: ad.AnnData, directory: Path, dataset: str) -> list[dict[str, Any]]:
    results: list[dict[str, Any]] = []
    processed: set[str] = set()
    for name in ["core3_structured_scores.csv"]:
        path = directory / name
        if not path.exists():
            continue
        table = pd.read_csv(path)
        matched = join_cell_table(adata, table, f"upstream_{dataset}_core3", name)
        results.append({"path": str(path), "kind": "cell_level", "n_rows": int(len(table)), "n_matched": matched})
        processed.add(name)

    program_path = directory / "pooled_program_scores.csv"
    if program_path.exists():
        table = pd.read_csv(program_path)
        if {"cell_id", "pool", "candidate_program", "score"}.issubset(table.columns):
            table["column"] = table["pool"].map(safe_column) + "_" + table["candidate_program"].map(safe_column)
            wide = table.pivot_table(index="cell_id", columns="column", values="score", aggfunc="mean").reset_index()
            matched = join_cell_table(adata, wide, f"upstream_{dataset}_program", "pooled_program_scores")
            results.append({"path": str(program_path), "kind": "cell_level_pivot", "n_rows": int(len(table)), "n_matched": matched})
        processed.add(program_path.name)

    # Donor-level outputs cannot be used as cell-level training labels, but
    # they are inherited for post-hoc interpretation through the canonical
    # donor_id key.  Keep the gene-level columns out of obs to avoid creating
    # enormous repeated strings; their source files remain listed below.
    donor_tables = {
        "donor_level_program_scores.csv": "donor_level",
        "donor_meta_program_matches.csv": "donor_meta",
        "donor_wise_program_genes.csv": "donor_wise",
    }
    for name, kind in donor_tables.items():
        path = directory / name
        if not path.exists():
            continue
        table = pd.read_csv(path)
        item = {"path": str(path), "kind": kind, "n_rows": int(len(table)), "n_matched": 0}
        if "donor_id" in table.columns and "donor_id" in adata.obs:
            exclude = {"gene", "weight", "rank_position"} if kind == "donor_wise" else set()
            item["n_matched"] = join_key_table(
                adata,
                table,
                "donor_id",
                "donor_id",
                f"upstream_{dataset}_{kind}",
                name,
                exclude_columns=exclude,
            )
        results.append(item)
        processed.add(name)

    # Record every other CSV (including meta-program summaries and sensitivity
    # tables) so the full upstream result inventory is preserved even when it
    # has no cell/donor join key.
    for path in sorted(directory.glob("*.csv")):
        if path.name in processed or path.name == "candidate_pool_labels.csv":
            continue
        try:
            header = pd.read_csv(path, nrows=0)
            n_rows = int(sum(1 for _ in path.open(encoding="utf-8"))) - 1
        except Exception as exc:
            results.append({"path": str(path), "kind": "unreadable", "error": f"{type(exc).__name__}: {exc}"})
            continue
        results.append({
            "path": str(path),
            "kind": "audit_only",
            "n_rows": max(0, n_rows),
            "columns": [str(col) for col in header.columns],
            "join_key": next((key for key in ["cell_id", "donor_id", "sample_id", "specimen_id"] if key in header.columns), None),
        })
    return results


def discover_annotation_files(config: dict[str, Any], dataset: str) -> list[dict[str, Any]]:
    directory = annotation_directory(config, dataset)
    if directory is None or not directory.exists():
        return []
    files = sorted(directory.glob("*.csv"))
    result = []
    for path in files:
        try:
            columns = [str(col) for col in pd.read_csv(path, nrows=0).columns]
        except Exception as exc:
            columns = []
            error = f"{type(exc).__name__}: {exc}"
        else:
            error = None
        result.append({
            "dataset": dataset,
            "path": str(path),
            "name": path.name,
            "available": True,
            "columns": columns,
            "error": error,
        })
    return result


def recreate_log_normalized_x(adata: ad.AnnData, target_sum: float = 10000.0) -> None:
    import scanpy as sc

    adata.X = adata.layers["counts"].copy()
    # Inherited AnnData often carries the upstream ``uns['log1p']`` marker.
    # Remove it before rebuilding X so scanpy does not treat raw counts as
    # already log-transformed.
    adata.uns.pop("log1p", None)
    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)
    sc.pp.log1p(adata)


def model_mask(adata: ad.AnnData, include_normal: bool = False) -> np.ndarray:
    roles = adata.obs["analysis_role"].astype(str)
    allowed = {"primary_candidate", "boundary"}
    if include_normal:
        allowed.add("normal_reference")
    return roles.isin(allowed).to_numpy()


def latent_leiden_labels(
    latent: np.ndarray,
    obs_index: pd.Index,
    n_neighbors: int,
    resolution: float,
    seed: int,
    key_added: str = "leiden",
) -> pd.Series:
    """Run a bounded Leiden graph on an already-computed latent matrix.

    The temporary AnnData has a one-column dummy X and uses only ``X_scVI``.
    In particular, this helper cannot accidentally consume counts, scaled X,
    PCA, or any upstream annotation.
    """
    import scanpy as sc

    latent = np.asarray(latent, dtype=np.float32)
    if latent.ndim != 2 or latent.shape[0] != len(obs_index):
        raise ValueError("latent and obs_index have incompatible shapes")
    if latent.shape[0] < 2:
        return pd.Series("0", index=obs_index, dtype=str)
    graph = ad.AnnData(
        X=np.zeros((latent.shape[0], 1), dtype=np.float32),
        obs=pd.DataFrame(index=pd.Index(obs_index).copy()),
    )
    graph.obsm["X_scVI"] = latent
    effective_neighbors = min(max(2, int(n_neighbors)), latent.shape[0] - 1)
    sc.pp.neighbors(
        graph,
        n_neighbors=effective_neighbors,
        use_rep="X_scVI",
        random_state=int(seed),
        key_added="neighbors_scvi_lam_only",
    )
    sc.tl.leiden(
        graph,
        resolution=float(resolution),
        key_added=key_added,
        neighbors_key="neighbors_scvi_lam_only",
        random_state=int(seed),
        flavor="igraph",
        directed=False,
    )
    return pd.Series(graph.obs[key_added].astype(str).to_numpy(), index=obs_index, dtype=str)


def partition_cluster_matches(labels_a: Iterable[Any], labels_b: Iterable[Any]) -> pd.DataFrame:
    """Return all non-zero cluster intersections, retaining split/merge matches."""
    a = pd.Series(np.asarray(list(labels_a), dtype=str), name="cluster_a")
    b = pd.Series(np.asarray(list(labels_b), dtype=str), name="cluster_b")
    if len(a) != len(b):
        raise ValueError("partitions must contain the same number of cells")
    frame = pd.DataFrame({"cluster_a": a, "cluster_b": b})
    intersections = frame.groupby(["cluster_a", "cluster_b"], observed=True).size().reset_index(name="intersection")
    if intersections.empty:
        return pd.DataFrame(columns=["cluster_a", "cluster_b", "size_a", "size_b", "intersection", "union", "jaccard", "containment_a", "containment_b"])
    sizes_a = frame.groupby("cluster_a", observed=True).size().rename("size_a")
    sizes_b = frame.groupby("cluster_b", observed=True).size().rename("size_b")
    result = intersections.join(sizes_a, on="cluster_a").join(sizes_b, on="cluster_b")
    result["union"] = result["size_a"] + result["size_b"] - result["intersection"]
    result["jaccard"] = result["intersection"] / result["union"].replace(0, np.nan)
    result["containment_a"] = result["intersection"] / result["size_a"].replace(0, np.nan)
    result["containment_b"] = result["intersection"] / result["size_b"].replace(0, np.nan)
    return result.reset_index(drop=True)


def pairwise_ari(labels_by_name: dict[str, Iterable[Any]]) -> pd.DataFrame:
    from sklearn.metrics import adjusted_rand_score

    names = list(labels_by_name)
    rows = []
    for i, left in enumerate(names):
        for right in names[i + 1:]:
            rows.append({
                "left": left,
                "right": right,
                "ari": float(adjusted_rand_score(
                    np.asarray(list(labels_by_name[left]), dtype=str),
                    np.asarray(list(labels_by_name[right]), dtype=str),
                )),
            })
    return pd.DataFrame(rows, columns=["left", "right", "ari"])


def mean_log1p_expression(adata: ad.AnnData, cell_mask: np.ndarray, gene_names: list[str]) -> dict[str, float]:
    """Compute mean normalized expression for a small gene set without densifying X."""
    if not gene_names or not cell_mask.any():
        return {}
    positions = [adata.var_names.get_loc(gene) for gene in gene_names if gene in adata.var_names]
    if not positions:
        return {}
    values = adata.X[cell_mask][:, positions]
    if sparse.issparse(values):
        values = values.toarray()
    values = np.asarray(values, dtype=float)
    return {adata.var_names[pos]: float(values[:, i].mean()) for i, pos in enumerate(positions)}
