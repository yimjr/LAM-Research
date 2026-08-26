"""Bridge immune-visibility genes to existing rapamycin retention results."""

from __future__ import annotations

import pandas as pd

from common import PROJECT_ROOT, ensure_output_path, load_signatures, load_source_manifest, project_relative, resolve_source, write_json


def main() -> None:
    manifest = load_source_manifest()
    signatures = load_signatures()
    retention_path = resolve_source(manifest["sirolimus"]["retention_gene"], manifest["source_root"])
    cross_path = resolve_source(manifest["sirolimus"]["retention_cross_environment"], manifest["source_root"])
    retention = pd.read_csv(retention_path)
    cross = pd.read_csv(cross_path) if cross_path.exists() else pd.DataFrame()
    gene_to_modules = {}
    for module, spec in signatures.items():
        for gene in spec.get("genes", []):
            gene_to_modules.setdefault(str(gene).upper(), []).append(module)
    retention["gene"] = retention["gene"].astype(str).str.upper()
    selected = retention[retention["gene"].isin(gene_to_modules)].copy()
    selected["modules"] = selected["gene"].map(lambda gene: ";".join(gene_to_modules[gene]))
    if not cross.empty:
        cross["gene"] = cross["gene"].astype(str).str.upper()
        cross = cross[cross["gene"].isin(gene_to_modules)].copy()
        cross["modules"] = cross["gene"].map(lambda gene: ";".join(gene_to_modules[gene]))
    output_dir = PROJECT_ROOT / "results" / "sirolimus_bridge"
    selected.to_csv(ensure_output_path(output_dir / "visibility_gene_retention_by_environment.csv"), index=False)
    if not cross.empty:
        cross.to_csv(ensure_output_path(output_dir / "visibility_cross_environment_retention.csv"), index=False)
    human_mapping = resolve_source(manifest["sirolimus"]["human_state_scores"], manifest["source_root"])
    module_mapping = resolve_source(manifest["sirolimus"]["human_module_scores"], manifest["source_root"])
    write_json(PROJECT_ROOT / "manifests" / "sirolimus_bridge_manifest.json", {
        "retention_source": project_relative(retention_path),
        "cross_environment_source": project_relative(cross_path),
        "visibility_genes_found": sorted(selected["gene"].unique().tolist()),
        "human_mapping_available": human_mapping.exists(),
        "human_module_mapping_available": module_mapping.exists(),
        "interpretation": "Existing perturbation retention supports candidate persistence only; it is not patient-level sirolimus treatment evidence.",
    })
    print(f"retention rows: {len(selected)}; cross-environment rows: {len(cross)}")


if __name__ == "__main__":
    main()
