"""Build complete file-level and dataset-level inventories for project data.

The inventory is derived from the actual ``data/`` tree, then enriched with
the existing download manifests and dataset configuration.  README files and
macOS metadata files are excluded; all actual data, source-document, archive,
annotation, and derived-input files are included.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd
import yaml

from common import ROOT, sha256_file


DATASET_META = {
    "GSE179044": {
        "role": "core_factorial_discovery",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE179044",
        "notes": "Complete WT/TSC2-null x vehicle/rapamycin x plastic/hydrogel factorial discovery dataset.",
    },
    "GSE27982": {
        "role": "external_tsc2_rapamycin_2x2",
        "organism": "Mus musculus",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE27982",
        "notes": "Low-serum mouse MEF external genotype-dependent rapamycin-response validation; not automatic escape proof.",
    },
    "GSE277844": {
        "role": "translation_residual_cross_model_analysis",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277844",
        "notes": "Human isogenic TSC2 WT/null neural progenitor-cell total and polysome counts; used for translation-program comparison with GSE179044 residuals, not direct LAM replication.",
    },
    "GSE16944": {
        "role": "historical_rapamycin_insensitive_support",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16944",
        "notes": "Historical LAM-like rapamycin-insensitive support; incomplete for formal residual/escape interaction.",
    },
    "GSE84476": {
        "role": "STAT3_mechanism_support",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84476",
        "notes": "STAT3 mechanism-support expression tables.",
    },
    "GSE104335": {
        "role": "SRPK2_mechanism_support",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE104335",
        "notes": "HTA-2.0 archive with CEL/CHP files; CHP gene-level extraction is derived locally.",
    },
    "GSE135851": {
        "role": "human_lamcore_validation",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135851",
        "notes": "Copied local snapshot with preliminary candidate/other labels, not formal LAMCORE/LAF taxonomy.",
    },
    "GSE302356": {
        "role": "human_lamcore_laf_spatial_validation",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE302356",
        "notes": "Multiple non-fully-paired modalities; LAM19 Xenium archive is staged but not in the current module-score pass.",
    },
    "GSE207322": {
        "role": "independent_btk_oriented_reference",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE207322",
        "notes": "TMD8 ibrutinib response with WT and BTK C481 mutants; independent pharmacological reference, not LAM.",
    },
    "GSE49414": {
        "role": "independent_ret_oriented_reference",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE49414",
        "notes": "TPC1 RPI-1 versus DMSO profile; RET-oriented but not RET-selective and not LAM.",
    },
    "GSE92742": {
        "role": "LINCS_cross_release_cmap_candidate_generation",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92742",
        "notes": "LINCS phase/release Level 5 GCTX and metadata; not independent biological replication.",
    },
    "GSE70138": {
        "role": "LINCS_cross_release_cmap_candidate_generation",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70138",
        "notes": "LINCS phase/release Level 5 GCTX and metadata; not independent biological replication.",
    },
    "GPL17518": {
        "role": "GSE49414_probe_annotation",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17518",
        "notes": "Complete Illumina HumanHT-12 V3 platform SOFT annotation.",
    },
    "GPL2895": {
        "role": "GSE16944_probe_annotation",
        "organism": "Homo sapiens",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL2895",
        "notes": "GEO platform annotation.",
    },
    "GPL339": {
        "role": "GSE27982_probe_annotation",
        "organism": "Mus musculus",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL339",
        "notes": "GEO platform annotation.",
    },
    "Gencode": {
        "role": "transcript_gene_annotation",
        "organism": "Homo sapiens",
        "source_url": "https://www.gencodegenes.org/human/",
        "notes": "GENCODE v24/v44 annotations used for transcript-to-gene support.",
    },
    "GO_Biological_Process_2023": {
        "role": "interpretive_enrichment_gene_sets",
        "organism": "Homo sapiens",
        "source_url": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Biological_Process_2023",
        "notes": "Downloaded GMT; used with 204-gene common analyzable background.",
    },
    "GO_Cellular_Component_2023": {
        "role": "interpretive_enrichment_gene_sets",
        "organism": "Homo sapiens",
        "source_url": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=GO_Cellular_Component_2023",
        "notes": "Downloaded GMT; used for gene-level annotation of the fixed 13-gene hydrogel translation-residual core.",
    },
    "Reactome_2022": {
        "role": "interpretive_enrichment_gene_sets",
        "organism": "Homo sapiens",
        "source_url": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=Reactome_2022",
        "notes": "Downloaded GMT; used with 204-gene common analyzable background.",
    },
    "MSigDB_Hallmark_2020": {
        "role": "interpretive_enrichment_gene_sets",
        "organism": "Homo sapiens",
        "source_url": "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020",
        "notes": "Downloaded GMT; used with 204-gene common analyzable background.",
    },
    "Mechanism_annotation": {
        "role": "derived_transcript_gene_map",
        "organism": "Homo sapiens",
        "source_url": "https://www.gencodegenes.org/human/",
        "notes": "Locally generated transcript-to-gene map from GENCODE v24.",
    },
}


def load_manifest_records() -> dict[str, dict]:
    records: dict[str, dict] = {}

    download = json.loads((ROOT / "manifests/download_manifest.json").read_text())
    for item in download.get("records", []):
        records[item["path"]] = {
            "source_accession": item["accession"],
            "checksum_type": "sha256",
            "checksum": item.get("sha256", ""),
            "checksum_status": "manifest_sha256",
            "manifest_status": item.get("status", "downloaded"),
        }

    lincs = json.loads((ROOT / "manifests/LINCS_download_manifest.json").read_text())
    for dataset_id, block in lincs.get("datasets", {}).items():
        for item in block.get("files", []):
            records[item["path"]] = {
                "source_accession": dataset_id,
                "checksum_type": "sha512",
                "checksum": item.get("sha512", ""),
                "checksum_status": "manifest_sha512",
                "manifest_status": "downloaded_verified_and_moved",
            }

    # Dataset-specific reference manifests store checksums for the current
    # external validation inputs.
    records["data/raw/GSE207322/GSE207322_TPM.txt.gz"] = {
        "source_accession": "GSE207322",
        "checksum_type": "sha256",
        "checksum": json.loads((ROOT / "manifests/GSE207322_BTK_reference.json").read_text())["matrix_sha256"],
        "checksum_status": "manifest_sha256",
        "manifest_status": "downloaded",
    }
    ret = json.loads((ROOT / "manifests/GSE49414_RET_reference.json").read_text())
    records["data/raw/GSE49414/GSE49414_series_matrix.txt.gz"] = {
        "source_accession": "GSE49414",
        "checksum_type": "sha256",
        "checksum": ret["matrix_sha256"],
        "checksum_status": "manifest_sha256",
        "manifest_status": "downloaded",
    }
    records["data/raw/GPL17518/GPL17518_family.soft.gz"] = {
        "source_accession": "GPL17518",
        "checksum_type": "sha256",
        "checksum": ret["annotation_sha256"],
        "checksum_status": "manifest_sha256",
        "manifest_status": "downloaded_verified_for_probe_mapping",
    }
    return records


def identify_dataset(relative_path: str) -> tuple[str, str]:
    parts = Path(relative_path).parts
    if parts[:2] == ("data", "raw") and len(parts) >= 3:
        group = parts[2]
        if group == "LINCS" and len(parts) >= 4:
            return parts[3], "raw"
        if group == "annotation":
            filename = Path(relative_path).name
            if "gencode.v24" in filename or "gencode.v44" in filename:
                return "Gencode", "raw"
            return "Annotation", "raw"
        if group.startswith("GPL"):
            return group, "raw"
        if group.startswith("GSE"):
            return group, "raw"
    if parts[:2] == ("data", "processed"):
        group = parts[2] if len(parts) >= 3 else "processed"
        if group == "gene_sets" and len(parts) >= 4:
            return Path(parts[3]).stem.replace(".gmt", ""), "processed"
        if group == "LINCS" and len(parts) >= 4 and parts[3] == "gctx":
            filename = parts[-1]
            if "GSE92742" in filename:
                return "GSE92742", "processed"
            if "GSE70138" in filename:
                return "GSE70138", "processed"
        if group == "mechanisms":
            return "Mechanism_annotation", "processed"
        if group.startswith("GSE"):
            return group, "processed"
    return "Unassigned", "unknown"


def classify_file(relative_path: str, dataset_id: str) -> tuple[str, str, str]:
    path = Path(relative_path)
    if "extracted_gene_expression" in path.parts or "unpacked" in path.parts:
        return "derived_extracted", "derived_local", "Extracted from a staged archive; source archive retained in the same dataset folder."
    if relative_path.startswith("data/processed/gctx/"):
        return "derived_processed", "derived_local", "Uncompressed/local working copy derived from the raw LINCS GCTX archive."
    if relative_path.startswith("data/processed/"):
        if dataset_id in {"GO_Biological_Process_2023", "Reactome_2022", "MSigDB_Hallmark_2020"}:
            return "external_download", "downloaded_external", "Downloaded GMT gene-set library."
        if dataset_id == "GSE135851":
            return "copied_input_snapshot", "copied_input", "Read-only local snapshot used for preliminary human mapping."
        if dataset_id == "Mechanism_annotation":
            return "derived_annotation", "derived_local", "Generated locally from retained GENCODE annotation."
        if path.name == "paper_state_marker_panels.csv":
            return "derived_marker_panel", "derived_from_paper", "Operational marker panel derived from the GSE302356 study report."
        return "derived_processed", "derived_local", "Derived analysis input."
    if path.name == "lam_niche_preprint.pdf" and path.stat().st_size < 100:
        return "invalid_placeholder", "invalid_or_incomplete", "Tiny error-response placeholder; not a valid PDF."
    if relative_path.startswith("data/raw/"):
        return "external_download", "downloaded_external", "Downloaded raw/archive/platform input."
    return "data_file", "unknown", "Data file under project data tree."


def source_url_for(dataset_id: str, relative_path: str) -> str:
    if dataset_id in DATASET_META:
        base = DATASET_META[dataset_id]["source_url"]
        if dataset_id == "GSE92742" and relative_path.endswith(".gctx.gz"):
            return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE927nnn/GSE92742/suppl/"
        if dataset_id == "GSE70138" and relative_path.endswith(".gctx.gz"):
            return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE701nnn/GSE70138/suppl/"
        return base
    return ""


def file_source_url_for(dataset_id: str, relative_path: str) -> str:
    filename = Path(relative_path).name
    if relative_path.startswith("data/processed/"):
        return source_url_for(dataset_id, relative_path)
    if dataset_id == "GSE179044" and filename == "GSE179044_raw_counts_matrix.csv.gz":
        return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179044/suppl/GSE179044_raw_counts_matrix.csv.gz"
    if dataset_id == "GSE27982" and filename == "GSE27982_series_matrix.txt.gz":
        return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27982/matrix/GSE27982_series_matrix.txt.gz"
    if dataset_id == "GSE277844" and filename == "GSE277844_raw_counts.txt.gz":
        return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277844/suppl/GSE277844_raw_counts.txt.gz"
    if dataset_id == "GSE16944" and filename == "GSE16944_series_matrix.txt.gz":
        return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE16nnn/GSE16944/matrix/GSE16944_series_matrix.txt.gz"
    if dataset_id == "GSE84476":
        return f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84476/suppl/{filename}"
    if dataset_id == "GSE104335" and filename == "GSE104335_RAW.tar":
        return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104335/suppl/GSE104335_RAW.tar"
    if dataset_id == "GSE302356":
        if filename == "GSE302356_family.xml.tgz":
            return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/GSE302356/miniml/GSE302356_family.xml.tgz"
        if filename.startswith("GSM") and filename.endswith(".tar.gz"):
            gsm = filename.split("_", 1)[0]
            return f"https://ftp.ncbi.nlm.nih.gov/geo/samples/{gsm[:-3]}nnn/{gsm}/suppl/{filename}"
    if dataset_id == "GSE207322":
        return f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE207nnn/GSE207322/{'suppl/' if filename.endswith('TPM.txt.gz') else 'matrix/'}{filename}"
    if dataset_id == "GSE49414" and filename == "GSE49414_series_matrix.txt.gz":
        return "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE49nnn/GSE49414/matrix/GSE49414_series_matrix.txt.gz"
    if dataset_id == "GSE92742":
        return f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE927nnn/GSE92742/suppl/{filename}"
    if dataset_id == "GSE70138":
        return f"https://ftp.ncbi.nlm.nih.gov/geo/series/GSE701nnn/GSE70138/suppl/{filename}"
    if dataset_id == "GPL17518":
        return "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL17nnn/GPL17518/soft/GPL17518_family.soft.gz" if filename.endswith(".soft.gz") else "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17518"
    if dataset_id == "GPL2895":
        return "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL2nnn/GPL2895/annot/GPL2895.annot.gz"
    if dataset_id == "GPL339":
        return "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL339/annot/GPL339.annot.gz"
    if dataset_id == "Gencode":
        if "v24" in filename:
            return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz"
        if "v44" in filename:
            return "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gtf.gz"
    return source_url_for(dataset_id, relative_path)


def build_file_inventory() -> pd.DataFrame:
    manifest_records = load_manifest_records()
    rows = []
    data_root = ROOT / "data"
    for path in sorted(data_root.rglob("*")):
        if not path.is_file() or path.name == ".DS_Store" or path.name == "README.md":
            continue
        relative = path.relative_to(ROOT).as_posix()
        dataset_id, tree_kind = identify_dataset(relative)
        category, provenance_status, verification_note = classify_file(relative, dataset_id)
        metadata = DATASET_META.get(dataset_id, {})
        record = manifest_records.get(relative, {})
        size = path.stat().st_size
        checksum_type = record.get("checksum_type", "")
        checksum = record.get("checksum", "")
        checksum_status = record.get("checksum_status", "")
        if not checksum and category not in {"invalid_placeholder"} and size <= 200 * 1024 * 1024:
            checksum_type = "sha256"
            checksum = sha256_file(path)
            checksum_status = "computed_sha256"
        elif not checksum and category == "invalid_placeholder":
            checksum_status = "not_applicable_invalid_placeholder"
        elif not checksum:
            checksum_status = "not_computed_large_file"
        if category == "invalid_placeholder":
            provenance_status = "invalid_or_incomplete"
        rows.append({
            "dataset_id": dataset_id,
            "tree_kind": tree_kind,
            "role": metadata.get("role", "supporting_or_derived_data"),
            "organism": metadata.get("organism", ""),
            "relative_path": relative,
            "filename": path.name,
            "file_category": category,
            "provenance_status": provenance_status,
            "size_bytes": size,
            "size_gib": round(size / (1024 ** 3), 6),
            "source_url": source_url_for(dataset_id, relative),
            "file_source_url": file_source_url_for(dataset_id, relative),
            "checksum_type": checksum_type,
            "checksum": checksum,
            "checksum_status": checksum_status,
            "manifest_status": record.get("manifest_status", ""),
            "verification_note": verification_note,
            "dataset_notes": metadata.get("notes", ""),
        })
    return pd.DataFrame(rows)


def build_dataset_inventory(files: pd.DataFrame) -> pd.DataFrame:
    rows = []
    observed = set(files["dataset_id"])
    for dataset_id in sorted(observed):
        group = files.loc[files["dataset_id"].eq(dataset_id)]
        metadata = DATASET_META.get(dataset_id, {})
        raw = group.loc[group["tree_kind"].eq("raw"), "size_bytes"].sum()
        processed = group.loc[group["tree_kind"].eq("processed"), "size_bytes"].sum()
        invalid = int(group["provenance_status"].eq("invalid_or_incomplete").sum())
        rows.append({
            "dataset_id": dataset_id,
            "role": metadata.get("role", "supporting_or_derived_data"),
            "organism": metadata.get("organism", ""),
            "source_url": metadata.get("source_url", ""),
            "n_files": int(len(group)),
            "n_external_download_files": int(group["provenance_status"].eq("downloaded_external").sum()),
            "n_derived_or_processed_files": int((~group["provenance_status"].eq("downloaded_external")).sum() - invalid),
            "n_invalid_or_incomplete_files": invalid,
            "total_size_bytes": int(group["size_bytes"].sum()),
            "total_size_gib": round(group["size_bytes"].sum() / (1024 ** 3), 6),
            "raw_size_bytes": int(raw),
            "processed_size_bytes": int(processed),
            "checksum_covered_files": int((group["checksum"].fillna("") != "").sum()),
            "notes": metadata.get("notes", ""),
        })
    return pd.DataFrame(rows)


def main() -> None:
    files = build_file_inventory()
    datasets = build_dataset_inventory(files)
    out = ROOT / "manifests"
    out.mkdir(parents=True, exist_ok=True)
    files.to_csv(out / "data_inventory.csv", index=False)
    datasets.to_csv(out / "dataset_inventory.csv", index=False)
    (out / "data_inventory_scope.json").write_text(json.dumps({
        "project_root": ".",
        "scope": "All actual files under data/ excluding README.md and .DS_Store; includes raw downloads, source documents, extracted archives, processed inputs, gene sets, and local derived data files.",
        "n_file_rows": int(len(files)),
        "n_dataset_rows": int(len(datasets)),
        "total_size_bytes": int(files["size_bytes"].sum()),
        "total_size_gib": round(files["size_bytes"].sum() / (1024 ** 3), 6),
        "invalid_or_incomplete_files": files.loc[files["provenance_status"].eq("invalid_or_incomplete"), "relative_path"].tolist(),
        "hash_policy": "Use existing manifest SHA256/SHA512 when available; compute SHA256 for files <=200 MiB; mark larger unrecorded files not_computed_large_file.",
    }, indent=2, ensure_ascii=False) + "\n")
    print({
        "n_file_rows": int(len(files)),
        "n_dataset_rows": int(len(datasets)),
        "total_size_gib": round(files["size_bytes"].sum() / (1024 ** 3), 3),
        "invalid_or_incomplete": files.loc[files["provenance_status"].eq("invalid_or_incomplete"), "relative_path"].tolist(),
    })


if __name__ == "__main__":
    main()
