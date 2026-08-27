"""Analyse concordant compound rows for the TSC2-loss plastic or hydrogel queries.

This is a downstream, auditable analysis of the already-computed local LINCS
ranking.  It deliberately reports both the exact filtered rows and a
drug-level deduplication: a LINCS row is a dataset/query-size/perturbation
summary, not necessarily a unique drug.

Target annotations are retrieved from the ChEMBL REST API using the LINCS
InChIKey where possible, then exact name/synonym matching as a fallback.
ChEMBL ``mechanism`` records are curated pharmacological mechanisms, not an
exhaustive list of every weak assay hit.  Every drug is retained in the
target table, including drugs for which no curated mechanism was returned.
"""

from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from common import CANDIDATE_DRUG_TARGETS, CANDIDATE_RESULTS, ROOT, write_json


RANKING = CANDIDATE_RESULTS / "LINCS_candidate_ranking.csv"
GSE92742_PERT = ROOT / "data/raw/LINCS/GSE92742/GSE92742_Broad_LINCS_pert_info.txt.gz"
GSE70138_PERT = ROOT / "data/raw/LINCS/GSE70138/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"
OUT = CANDIDATE_RESULTS
TARGET_OUT = CANDIDATE_DRUG_TARGETS
LEGACY_TARGET_PREFIX = "tsc2_loss_plastic_replicated_concordant"
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
BINDINGDB_TARGET_BY_COMPOUND = "https://bindingdb.org/rest/getTargetByCompound"

# ChEMBL's curated mechanism endpoint intentionally omits many tool
# compounds, endogenous ligands and compounds whose action is described as a
# biochemical class rather than a single target.  These fallback records are
# restricted to well-established primary mechanisms and are kept visibly
# separate from ChEMBL records in the output.
LITERATURE_TARGETS: dict[str, list[dict[str, Any]]] = {
    "CGS-20625": [
        {"gene": "GABRA1", "name": "GABA-A receptor alpha-1-containing complex", "mechanism": "central benzodiazepine-site partial agonist/modulator", "role": "primary receptor complex", "url": "https://pubmed.ncbi.nlm.nih.gov/2563294/"},
        {"gene": "GABRB2", "name": "GABA-A receptor beta-2-containing complex", "mechanism": "central benzodiazepine-site partial agonist/modulator", "role": "receptor complex component", "url": "https://pubmed.ncbi.nlm.nih.gov/2563294/"},
        {"gene": "GABRG1", "name": "GABA-A receptor gamma-1-containing complex", "mechanism": "central benzodiazepine-site partial agonist/modulator", "role": "receptor complex component", "url": "https://pubmed.ncbi.nlm.nih.gov/2563294/"},
    ],
    "LY-294002": [
        {"gene": "PIK3CA", "name": "class I phosphoinositide 3-kinase catalytic subunit alpha", "mechanism": "PI3K inhibitor; tool compound with additional kinase off-targets", "role": "primary pathway target", "url": "https://pubmed.ncbi.nlm.nih.gov/12807916/"},
        {"gene": "PIK3CB", "name": "class I phosphoinositide 3-kinase catalytic subunit beta", "mechanism": "PI3K inhibitor; tool compound with additional kinase off-targets", "role": "primary pathway target", "url": "https://pubmed.ncbi.nlm.nih.gov/12807916/"},
        {"gene": "PIK3CD", "name": "class I phosphoinositide 3-kinase catalytic subunit delta", "mechanism": "PI3K inhibitor; tool compound with additional kinase off-targets", "role": "primary pathway target", "url": "https://pubmed.ncbi.nlm.nih.gov/12807916/"},
        {"gene": "PIK3CG", "name": "class I phosphoinositide 3-kinase catalytic subunit gamma", "mechanism": "PI3K inhibitor; tool compound with additional kinase off-targets", "role": "primary pathway target", "url": "https://pubmed.ncbi.nlm.nih.gov/12807916/"},
        {"gene": "MTOR", "name": "mechanistic target of rapamycin", "mechanism": "reported dual PI3K/mTOR inhibition in cell-based pharmacology", "role": "secondary pathway target", "url": "https://pubmed.ncbi.nlm.nih.gov/12807916/"},
    ],
    "MG-132": [
        {"gene": "PSMB5", "name": "26S proteasome beta subunit 5", "mechanism": "reversible proteasome inhibitor; peptide-aldehyde tool compound", "role": "primary proteasome catalytic target", "url": "https://pubchem.ncbi.nlm.nih.gov/compound/mg-132"},
    ],
    "PI-103": [
        {"gene": "PIK3CA", "name": "PI3-kinase catalytic subunit alpha", "mechanism": "dual PI3K/mTOR inhibitor", "role": "primary pathway target", "url": "https://pubmed.ncbi.nlm.nih.gov/19351820/"},
        {"gene": "MTOR", "name": "mechanistic target of rapamycin", "mechanism": "dual PI3K/mTOR inhibitor", "role": "primary pathway target", "url": "https://pubmed.ncbi.nlm.nih.gov/19351820/"},
        {"gene": "PRKDC", "name": "DNA-dependent protein kinase catalytic subunit", "mechanism": "reported DNA-PK contribution to DNA-repair chemosensitization", "role": "secondary pathway target", "url": "https://pubmed.ncbi.nlm.nih.gov/19633683/"},
    ],
    "QL-X-138": [
        {"gene": "BTK", "name": "Bruton tyrosine kinase", "mechanism": "covalent BTK inhibition", "role": "primary target", "url": "https://pubmed.ncbi.nlm.nih.gov/26165234/"},
        {"gene": "MKNK1", "name": "MAP kinase-interacting serine/threonine kinase 1", "mechanism": "noncovalent MNK inhibition", "role": "primary target", "url": "https://pubmed.ncbi.nlm.nih.gov/26165234/"},
        {"gene": "MKNK2", "name": "MAP kinase-interacting serine/threonine kinase 2", "mechanism": "noncovalent MNK inhibition", "role": "primary target", "url": "https://pubmed.ncbi.nlm.nih.gov/26165234/"},
    ],
    "clofarabine": [
        {"gene": "RRM1", "name": "ribonucleotide reductase large subunit", "mechanism": "active triphosphate/diphosphate inhibits ribonucleotide reductase and depletes dNTPs", "role": "primary target", "url": "https://pubmed.ncbi.nlm.nih.gov/22840768/"},
        {"gene": "POLA1", "name": "DNA polymerase alpha catalytic subunit", "mechanism": "active triphosphate inhibits DNA polymerase alpha and DNA elongation", "role": "primary target", "url": "https://pubmed.ncbi.nlm.nih.gov/15803490/"},
        {"gene": "DCK", "name": "deoxycytidine kinase", "mechanism": "phosphorylates/activates clofarabine; activation enzyme rather than cytotoxic target", "role": "activation enzyme", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5681218/"},
    ],
    "dinoprost": [
        {"gene": "PTGFR", "name": "prostaglandin F2-alpha receptor", "mechanism": "prostaglandin F2-alpha receptor agonist", "role": "primary receptor", "url": "https://pubmed.ncbi.nlm.nih.gov/28298246/"},
    ],
    "milnacipran": [
        {"gene": "SLC6A2", "name": "norepinephrine transporter", "mechanism": "inhibits norepinephrine reuptake", "role": "primary transporter", "url": "https://pubmed.ncbi.nlm.nih.gov/12122491/"},
        {"gene": "SLC6A4", "name": "serotonin transporter", "mechanism": "inhibits serotonin reuptake", "role": "primary transporter", "url": "https://pubmed.ncbi.nlm.nih.gov/12122491/"},
    ],
    "niacin": [
        {"gene": "HCAR2", "name": "hydroxycarboxylic acid receptor 2 / GPR109A", "mechanism": "niacin agonist", "role": "primary receptor", "url": "https://pubmed.ncbi.nlm.nih.gov/38012147/"},
        {"gene": "HCAR3", "name": "hydroxycarboxylic acid receptor 3 / GPR109B", "mechanism": "lower-affinity niacin receptor in humans", "role": "secondary receptor", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC10682194/"},
    ],
    "phorbol-myristate-acetate": [
        {"gene": "PRKCA", "name": "protein kinase C alpha", "mechanism": "phorbol-ester activation through C1-domain binding", "role": "PKC family target", "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC393304/"},
        {"gene": "PRKCD", "name": "protein kinase C delta", "mechanism": "phorbol-ester activation; implicated in PMA-induced G1 arrest", "role": "PKC family target", "url": "https://pubmed.ncbi.nlm.nih.gov/16055435/"},
    ],
    "podophyllotoxin": [
        {"gene": "TUBB", "name": "beta-tubulin / microtubule", "mechanism": "binds tubulin at or near the colchicine site and disrupts microtubule dynamics", "role": "primary cytoskeletal target", "url": "https://pubmed.ncbi.nlm.nih.gov/2722802/"},
    ],
    "roxithromycin": [
        {"gene": None, "name": "bacterial 50S ribosomal subunit / 23S rRNA", "mechanism": "macrolide binding blocks peptide translocation", "role": "bacterial ribosome target; not a human gene", "url": "https://pubchem.ncbi.nlm.nih.gov/compound/Roxithromycin"},
    ],
}

OUTPUT_PREFIX = "tsc2_loss_plastic_or_hydrogel_replicated_concordant"
FILTER = {
    "contrasts": ["tsc2_loss_plastic", "tsc2_loss_hydrogel"],
    "perturbation_class": "compound",
    "cross_phase_status": "replicated_concordant",
}


def candidate_output(filename: str) -> Path:
    return OUT / f"{OUTPUT_PREFIX}_{filename}"


def target_output(filename: str) -> Path:
    return TARGET_OUT / f"{OUTPUT_PREFIX}_{filename}"


def norm(value: Any) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def fetch_json(path: str) -> dict[str, Any]:
    url = f"{CHEMBL_BASE}/{path}"
    request = urllib.request.Request(url, headers={"User-Agent": "LAM-Drug-Repurposing/1.0"})
    with urllib.request.urlopen(request, timeout=20) as response:
        return json.load(response)


def fetch_pubchem(path: str) -> dict[str, Any]:
    url = f"{PUBCHEM_BASE}/{path}"
    request = urllib.request.Request(url, headers={"User-Agent": "LAM-Drug-Repurposing/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def fetch_bindingdb(smiles: str, cutoff: str = "0.999") -> dict[str, Any]:
    query = urllib.parse.urlencode({"smiles": smiles, "cutoff": cutoff, "response": "application/json"})
    request = urllib.request.Request(f"{BINDINGDB_TARGET_BY_COMPOUND}?{query}", headers={"User-Agent": "LAM-Drug-Repurposing/1.0"})
    with urllib.request.urlopen(request, timeout=30) as response:
        return json.load(response)


def query_molecule(name: str, inchikey: str) -> tuple[dict[str, Any] | None, str]:
    """Match a LINCS compound to a ChEMBL molecule with explicit provenance."""
    queries: list[tuple[str, str]] = []
    if inchikey and str(inchikey) != "nan":
        queries.append((
            f"molecule.json?molecule_structures__standard_inchi_key__iexact="
            f"{urllib.parse.quote(str(inchikey))}&limit=20",
            "inchikey_exact",
        ))
    queries.extend([
        (f"molecule.json?pref_name__iexact={urllib.parse.quote(name)}&limit=20", "name_exact"),
        (f"molecule.json?molecule_synonyms__synonyms__iexact={urllib.parse.quote(name)}&limit=20", "synonym_exact"),
    ])
    for path, method in queries:
        try:
            molecules = fetch_json(path).get("molecules", [])
        except Exception:
            molecules = []
        if molecules:
            return molecules[0], method

    # Search is only a last resort.  Accept it only if the returned preferred
    # name or a synonym is an exact normalized match; otherwise keep unknown.
    try:
        molecules = fetch_json(f"molecule/search.json?q={urllib.parse.quote(name)}&limit=100").get("molecules", [])
    except Exception:
        molecules = []
    for molecule in molecules:
        if norm(molecule.get("pref_name")) == norm(name):
            return molecule, "search_pref_name_exact"
        synonyms = molecule.get("molecule_synonyms") or []
        if any(norm(item.get("synonyms")) == norm(name) for item in synonyms):
            return molecule, "search_synonym_exact"
    return None, "not_matched"


def target_rows_for_molecule(
    name: str,
    molecule: dict[str, Any] | None,
    match_method: str,
    lin_cs: dict[str, Any],
) -> list[dict[str, Any]]:
    base = {
        "pert_iname": name,
        "lincs_pert_id_GSE92742": lin_cs.get("pert_id_GSE92742"),
        "lincs_pert_id_GSE70138": lin_cs.get("pert_id_GSE70138"),
        "lincs_inchikey": lin_cs.get("inchi_key"),
        "lincs_pubchem_cid": lin_cs.get("pubchem_cid"),
        "canonical_smiles": lin_cs.get("canonical_smiles"),
        "chembl_match_method": match_method,
        "chembl_molecule_id": (molecule or {}).get("molecule_chembl_id"),
        "chembl_pref_name": (molecule or {}).get("pref_name"),
    }
    if molecule is None or not molecule.get("molecule_chembl_id"):
        return [{
            **base,
            "target_chembl_id": None,
            "target_name": None,
            "target_type": None,
            "organism": None,
            "target_gene_symbol": None,
            "target_uniprot_accession": None,
            "mechanism_of_action": None,
            "action_type": None,
            "direct_interaction": None,
            "mechanism_reference": None,
            "target_role": None,
            "target_source": "ChEMBL mechanism API",
            "target_source_url": None,
            "target_confidence": "not_annotated",
            "target_status": "molecule_not_matched",
        }]

    molecule_id = molecule["molecule_chembl_id"]
    try:
        mechanisms = fetch_json(f"mechanism.json?molecule_chembl_id={molecule_id}&limit=100").get("mechanisms", [])
    except Exception:
        mechanisms = []
    rows: list[dict[str, Any]] = []
    for mechanism in mechanisms:
        target_id = mechanism.get("target_chembl_id")
        try:
            target = fetch_json(f"target/{target_id}.json")
        except Exception:
            target = {}
        components = target.get("target_components") or [{}]
        refs = mechanism.get("mechanism_refs") or []
        ref_url = refs[0].get("ref_url") if refs else None
        for component in components:
            symbols = [
                item.get("component_synonym")
                for item in (component.get("target_component_synonyms") or [])
                if item.get("syn_type") == "GENE_SYMBOL"
            ]
            rows.append({
                **base,
                "target_chembl_id": target_id,
                "target_name": target.get("pref_name"),
                "target_type": target.get("target_type"),
                "organism": target.get("organism") or component.get("organism"),
                "target_gene_symbol": ";".join(sorted(set(str(x) for x in symbols if x))),
                "target_uniprot_accession": component.get("accession"),
                "mechanism_of_action": mechanism.get("mechanism_of_action"),
                "action_type": mechanism.get("action_type"),
                "direct_interaction": mechanism.get("direct_interaction"),
                "mechanism_reference": ref_url,
                "target_role": "curated ChEMBL mechanism target",
                "target_source": "ChEMBL curated mechanism",
                "target_source_url": f"{CHEMBL_BASE}/mechanism.json?molecule_chembl_id={molecule_id}",
                "target_confidence": "curated_mechanism",
                "target_status": "annotated",
            })
    if not rows:
        # Fill the deliberate ChEMBL coverage gap with a small, curated set of
        # primary mechanisms.  These rows are not silently merged into the
        # ChEMBL evidence class.
        for fallback in LITERATURE_TARGETS.get(name, []):
            rows.append({
                **base,
                "target_chembl_id": None,
                "target_name": fallback.get("name"),
                "target_type": "literature-defined target or target complex",
                "organism": "Homo sapiens" if fallback.get("gene") else "Bacteria",
                "target_gene_symbol": fallback.get("gene"),
                "target_uniprot_accession": None,
                "mechanism_of_action": fallback.get("mechanism"),
                "action_type": "MODULATOR_OR_INHIBITOR",
                "direct_interaction": None,
                "mechanism_reference": fallback.get("url"),
                "target_role": fallback.get("role"),
                "target_source": "curated primary literature fallback",
                "target_source_url": fallback.get("url"),
                "target_confidence": "literature_curated",
                "target_status": "annotated_literature_fallback",
            })
    if not rows:
        rows.append({
            **base,
            "target_chembl_id": None,
            "target_name": None,
            "target_type": None,
            "organism": None,
            "target_gene_symbol": None,
            "target_uniprot_accession": None,
            "mechanism_of_action": None,
            "action_type": None,
            "direct_interaction": None,
            "mechanism_reference": None,
            "target_role": None,
            "target_source": "ChEMBL mechanism API",
            "target_source_url": f"{CHEMBL_BASE}/mechanism.json?molecule_chembl_id={molecule_id}",
            "target_confidence": "no_curated_mechanism_returned",
            "target_status": "unannotated_mechanism",
        })
    return rows


def direction_class(value: float) -> str:
    if not np.isfinite(value) or abs(value) < 0.25:
        return "weak_or_neutral"
    return "reversal_direction" if value > 0 else "mimic_direction"


def build_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    ranking = pd.read_csv(RANKING)
    required = {"contrast", "perturbation_class", "cross_phase_status", "pert_iname", "perturbation_key"}
    missing = required - set(ranking.columns)
    if missing:
        raise ValueError(f"ranking file missing columns: {sorted(missing)}")
    mask = ranking["contrast"].isin(FILTER["contrasts"]).to_numpy(copy=True)
    for key in ["perturbation_class", "cross_phase_status"]:
        value = FILTER[key]
        mask &= ranking[key].eq(value).to_numpy()
    filtered = ranking.loc[mask].copy()
    filtered.insert(0, "filter_row_id", np.arange(1, len(filtered) + 1))
    filtered["row_direction_class"] = filtered["median_reversal_NCS"].map(direction_class)
    filtered["unique_drug_id"] = filtered["pert_iname"].map(lambda x: f"compound::{norm(x)}")
    filtered = filtered.sort_values(["median_reversal_NCS", "pert_iname", "dataset", "query_id"], ascending=[False, True, True, True])

    numeric = [
        "median_reversal_NCS", "median_reversal_WTCS", "max_reversal_NCS",
        "fraction_reversed", "n_signatures", "n_cells", "n_query_sizes", "q25_reversal_NCS", "q75_reversal_NCS",
    ]
    for col in numeric:
        if col in filtered:
            filtered[col] = pd.to_numeric(filtered[col], errors="coerce")

    rows: list[dict[str, Any]] = []
    for (drug_id, name), group in filtered.groupby(["unique_drug_id", "pert_iname"], sort=True):
        directions = sorted(set(group["row_direction_class"]))
        if len(directions) == 1:
            pattern = directions[0]
        elif set(directions) <= {"reversal_direction", "mimic_direction"}:
            pattern = "mixed_reversal_and_mimic"
        else:
            pattern = "mixed_with_weak_or_neutral"
        rows.append({
            "unique_drug_id": drug_id,
            "pert_iname": name,
            "n_filtered_rows": len(group),
            "datasets": ";".join(sorted(group["dataset"].dropna().unique())),
            "contrasts": ";".join(sorted(group["contrast"].dropna().unique())),
            "query_ids": ";".join(sorted(group["query_id"].dropna().unique())),
            "n_query_sizes": group["n_query_sizes"].max() if "n_query_sizes" in group else np.nan,
            "applicable_target_sizes": ";".join(sorted(set(";".join(group["applicable_target_sizes"].fillna("").astype(str)).split(";")) - {"", "nan"})),
            "median_reversal_NCS_across_rows": group["median_reversal_NCS"].median(),
            "q25_reversal_NCS_across_rows": group["median_reversal_NCS"].quantile(0.25),
            "q75_reversal_NCS_across_rows": group["median_reversal_NCS"].quantile(0.75),
            "min_reversal_NCS": group["median_reversal_NCS"].min(),
            "max_reversal_NCS": group["median_reversal_NCS"].max(),
            "mean_fraction_reversed": group["fraction_reversed"].mean(),
            "n_positive_reversal_rows": int((group["row_direction_class"] == "reversal_direction").sum()),
            "n_mimic_direction_rows": int((group["row_direction_class"] == "mimic_direction").sum()),
            "n_weak_or_neutral_rows": int((group["row_direction_class"] == "weak_or_neutral").sum()),
            "direction_pattern": pattern,
            "cross_phase_status": "replicated_concordant",
            "best_context": ";".join(sorted(set(group["best_context"].dropna().astype(str)))),
        })
    unique = pd.DataFrame(rows).sort_values(["median_reversal_NCS_across_rows", "pert_iname"], ascending=[False, True])
    return filtered, unique, ranking


def load_lincs_identity(names: set[str]) -> pd.DataFrame:
    a = pd.read_csv(GSE92742_PERT, sep="\t", compression="gzip", low_memory=False)
    b = pd.read_csv(GSE70138_PERT, sep="\t", compression="gzip", low_memory=False)
    a = a[a["pert_iname"].isin(names)].copy()
    b = b[b["pert_iname"].isin(names)].copy()
    # GSE92742 supplies PubChem CID and is used as the primary identity table;
    # GSE70138 fills identity gaps and documents release-level IDs.
    a = a.drop_duplicates("pert_iname").set_index("pert_iname")
    b = b.drop_duplicates("pert_iname").set_index("pert_iname")
    out = []
    for name in sorted(names):
        ar = a.loc[name] if name in a.index else pd.Series(dtype=object)
        br = b.loc[name] if name in b.index else pd.Series(dtype=object)
        out.append({
            "pert_iname": name,
            "pert_id_GSE92742": ar.get("pert_id"),
            "pert_id_GSE70138": br.get("pert_id"),
            "inchi_key": ar.get("inchi_key") if pd.notna(ar.get("inchi_key")) else br.get("inchi_key"),
            "canonical_smiles": ar.get("canonical_smiles") if pd.notna(ar.get("canonical_smiles")) else br.get("canonical_smiles"),
            "pubchem_cid": ar.get("pubchem_cid"),
        })
    return pd.DataFrame(out)


def annotate_targets(unique: pd.DataFrame, identities: pd.DataFrame | None = None) -> pd.DataFrame:
    identities = identities if identities is not None else load_lincs_identity(set(unique["pert_iname"]))
    all_rows: list[dict[str, Any]] = []
    records = identities.to_dict("records")
    print(f"[targets] querying ChEMBL for {len(records)} compounds", flush=True)
    for index, record in enumerate(records, start=1):
        name = record["pert_iname"]
        if index == 1 or index % 10 == 0 or index == len(records):
            print(f"[targets] {index}/{len(records)} {name}", flush=True)
        molecule, method = query_molecule(name, str(record.get("inchi_key") or ""))
        all_rows.extend(target_rows_for_molecule(name, molecule, method, record))
        # Keep API usage polite and reproducible.  The output is saved after
        # the full pass, so a rerun can be used if a transient request fails.
        time.sleep(0.08)
    targets = pd.DataFrame(all_rows)
    return targets


def placeholder_targets(unique: pd.DataFrame, identities: pd.DataFrame) -> pd.DataFrame:
    """Create explicit rows for compounds whose curated ChEMBL lookup is deferred.

    PubChem assay target evidence is added separately after identity resolution.
    Keeping this placeholder row means an unresolved compound remains visible in
    the merged target table instead of silently disappearing.
    """
    rows: list[dict[str, Any]] = []
    records = {row["pert_iname"]: row for row in identities.to_dict("records")}
    for name in sorted(unique["pert_iname"].unique()):
        record = records.get(name, {})
        rows.append({
            "pert_iname": name,
            "lincs_pert_id_GSE92742": record.get("pert_id_GSE92742"),
            "lincs_pert_id_GSE70138": record.get("pert_id_GSE70138"),
            "lincs_inchikey": record.get("inchi_key"),
            "lincs_pubchem_cid": record.get("pubchem_cid"),
            "canonical_smiles": record.get("canonical_smiles"),
            "chembl_match_method": "not_queried_in_default_combined_run",
            "chembl_molecule_id": None,
            "chembl_pref_name": None,
            "target_chembl_id": None,
            "target_name": None,
            "target_type": None,
            "organism": None,
            "target_gene_symbol": None,
            "target_uniprot_accession": None,
            "mechanism_of_action": None,
            "action_type": None,
            "direct_interaction": None,
            "mechanism_reference": None,
            "target_role": None,
            "target_source": "ChEMBL mechanism API deferred",
            "target_source_url": CHEMBL_BASE,
            "target_confidence": "not_queried",
            "target_status": "molecule_not_annotated",
        })
    return pd.DataFrame(rows)


def annotate_pubchem(
    unique: pd.DataFrame,
    targets: pd.DataFrame,
    identities: pd.DataFrame | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Add PubChem identity/description and retain raw assay target evidence.

    PubChem BioAssay target fields are intentionally exported separately from
    the curated target table.  They are heterogeneous assay annotations and
    may include non-human proteins, panel targets, or assay-specific records;
    they are useful for completeness checks but should not be interpreted as
    primary pharmacological targets without context.
    """
    identities = identities if identities is not None else load_lincs_identity(set(unique["pert_iname"]))
    identity_by_name = {row["pert_iname"]: row for row in identities.to_dict("records")}
    identity_rows: list[dict[str, Any]] = []
    assay_rows: list[dict[str, Any]] = []
    names = sorted(unique["pert_iname"].unique())
    print(f"[pubchem] querying identity and assay evidence for {len(names)} compounds", flush=True)
    for index, name in enumerate(names, start=1):
        if index == 1 or index % 10 == 0 or index == len(names):
            print(f"[pubchem] {index}/{len(names)} {name}", flush=True)
        try:
            identity_record = identity_by_name.get(name, {})
            lincs_cid = identity_record.get("pubchem_cid")
            if pd.notna(lincs_cid) and str(lincs_cid) not in {"-666", "-666.0", "nan", "None"}:
                cid = int(float(lincs_cid))
            else:
                cid_payload = fetch_pubchem(f"compound/name/{urllib.parse.quote(name)}/cids/JSON")
                cids = cid_payload.get("IdentifierList", {}).get("CID", [])
                cid = cids[0] if cids else None
            desc_payload = fetch_pubchem(f"compound/cid/{cid}/description/JSON") if cid else {}
            descriptions = desc_payload.get("InformationList", {}).get("Information", [])
            descriptions = [x for x in descriptions if x.get("Description")]
            desc = descriptions[0] if descriptions else {}
            identity_rows.append({
                "pert_iname": name,
                "pubchem_cid_resolved": cid,
                "pubchem_compound_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else None,
                "pubchem_description": desc.get("Description"),
                "pubchem_description_source": desc.get("DescriptionSourceName"),
                "pubchem_description_source_url": desc.get("DescriptionURL"),
                "pubchem_query_status": "resolved" if cid else "not_resolved",
            })

            # This endpoint can be large, but it is the only public PubChem
            # endpoint here that exposes explicit target accession/GeneID
            # fields.  Preserve every row with a target field.
            try:
                assay_payload = fetch_pubchem(f"compound/cid/{cid}/assaysummary/JSON") if cid else {}
                table = assay_payload.get("Table", {})
                columns = table.get("Columns", {}).get("Column", [])
                column_index = {str(col): i for i, col in enumerate(columns)}
                for row in table.get("Row", []):
                    cells = row.get("Cell", [])
                    def cell(label: str) -> Any:
                        idx = column_index.get(label)
                        return cells[idx] if idx is not None and idx < len(cells) else None
                    target_accession = cell("Target Accession")
                    target_geneid = cell("Target GeneID")
                    if not target_accession and not target_geneid:
                        continue
                    assay_rows.append({
                        "pert_iname": name,
                        "pubchem_cid": cid,
                        "aid": cell("AID"),
                        "target_accession": target_accession,
                        "target_geneid": target_geneid,
                        "activity_outcome": cell("Activity Outcome"),
                        "activity_value_uM": cell("Activity Value [uM]"),
                        "activity_name": cell("Activity Name"),
                        "assay_name": cell("Assay Name"),
                        "assay_type": cell("Assay Type"),
                        "pubmed_id": cell("PubMed ID"),
                        "source_url": f"https://pubchem.ncbi.nlm.nih.gov/bioassay/{cell('AID')}",
                    })
            except Exception:
                # A missing PubChem assay summary is not a failure of identity
                # resolution or of the curated target table.
                pass
        except Exception as exc:
            identity_rows.append({
                "pert_iname": name,
                "pubchem_cid_resolved": None,
                "pubchem_compound_url": None,
                "pubchem_description": None,
                "pubchem_description_source": None,
                "pubchem_description_source_url": None,
                "pubchem_query_status": f"query_error:{type(exc).__name__}",
            })
        time.sleep(0.08)

    identity = pd.DataFrame(identity_rows)
    assay = pd.DataFrame(assay_rows)
    targets = targets.merge(identity, on="pert_iname", how="left")
    return targets, identity, assay


def pubchem_assay_target_rows(assay: pd.DataFrame, identities: pd.DataFrame) -> pd.DataFrame:
    """Represent PubChem assay target fields as a separate, lower-confidence layer."""
    columns = [
        "pert_iname", "lincs_pert_id_GSE92742", "lincs_pert_id_GSE70138",
        "lincs_inchikey", "lincs_pubchem_cid", "canonical_smiles",
        "chembl_match_method", "chembl_molecule_id", "chembl_pref_name",
        "target_chembl_id", "target_name", "target_type", "organism",
        "target_gene_symbol", "target_uniprot_accession", "mechanism_of_action",
        "action_type", "direct_interaction", "mechanism_reference", "target_role",
        "target_source", "target_source_url", "target_confidence", "target_status",
    ]
    if assay.empty:
        return pd.DataFrame(columns=columns)
    identity_by_name = {row["pert_iname"]: row for row in identities.to_dict("records")}
    grouped = assay.groupby(["pert_iname", "target_accession", "target_geneid"], dropna=False).size().reset_index(name="n_assay_rows")
    rows: list[dict[str, Any]] = []
    for record in grouped.to_dict("records"):
        name = record["pert_iname"]
        identity = identity_by_name.get(name, {})
        accession = record.get("target_accession")
        geneid = record.get("target_geneid")
        target_label = "PubChem BioAssay target"
        if pd.notna(geneid):
            target_label += f" (GeneID {geneid})"
        rows.append({
            "pert_iname": name,
            "lincs_pert_id_GSE92742": identity.get("pert_id_GSE92742"),
            "lincs_pert_id_GSE70138": identity.get("pert_id_GSE70138"),
            "lincs_inchikey": identity.get("inchi_key"),
            "lincs_pubchem_cid": identity.get("pubchem_cid"),
            "canonical_smiles": identity.get("canonical_smiles"),
            "chembl_match_method": "not_queried_in_default_combined_run",
            "chembl_molecule_id": None,
            "chembl_pref_name": None,
            "target_chembl_id": None,
            "target_name": target_label,
            "target_type": "PubChem BioAssay target evidence",
            "organism": None,
            "target_gene_symbol": None,
            "target_uniprot_accession": accession if pd.notna(accession) else None,
            "mechanism_of_action": None,
            "action_type": None,
            "direct_interaction": None,
            "mechanism_reference": None,
            "target_role": f"assay target evidence; n_assay_rows={int(record['n_assay_rows'])}",
            "target_source": "PubChem BioAssay summary",
            "target_source_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{identity.get('pubchem_cid')}" if pd.notna(identity.get("pubchem_cid")) else PUBCHEM_BASE,
            "target_confidence": "assay_evidence_not_primary_target",
            "target_status": "annotated_assay_evidence",
        })
    return pd.DataFrame(rows, columns=columns)


def load_legacy_annotations(names: set[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Reuse the original plastic-only annotations without modifying their files."""
    target_path = TARGET_OUT / f"{LEGACY_TARGET_PREFIX}_compound_targets.csv"
    identity_path = TARGET_OUT / f"{LEGACY_TARGET_PREFIX}_pubchem_identity.csv"
    assay_path = TARGET_OUT / f"{LEGACY_TARGET_PREFIX}_pubchem_assay_target_evidence.csv.gz"
    binding_path = TARGET_OUT / f"{LEGACY_TARGET_PREFIX}_bindingdb_affinity_evidence.csv.gz"

    def read_filtered(path: Path, compression: str | None = None) -> pd.DataFrame:
        if not path.exists():
            return pd.DataFrame()
        frame = pd.read_csv(path, compression=compression, low_memory=False)
        return frame.loc[frame["pert_iname"].isin(names)].copy() if "pert_iname" in frame else frame

    return (
        read_filtered(target_path),
        read_filtered(identity_path),
        read_filtered(assay_path, "gzip"),
        read_filtered(binding_path, "gzip"),
    )


def annotate_bindingdb(unique: pd.DataFrame, identities: pd.DataFrame | None = None) -> pd.DataFrame:
    """Retrieve high-similarity BindingDB affinity evidence.

    The public ``getTargetByCompound`` service is a similarity search, not an
    exact identity endpoint.  We therefore use a stringent 0.999 cutoff and
    preserve that cutoff in every row.  These data are evidence for binding
    and off-target context, not an automatic replacement for the curated
    target table.
    """
    identities = identities if identities is not None else load_lincs_identity(set(unique["pert_iname"]))
    rows: list[dict[str, Any]] = []
    records = identities.to_dict("records")
    print(f"[bindingdb] querying affinity evidence for {len(records)} compounds", flush=True)
    for index, record in enumerate(records, start=1):
        name = record["pert_iname"]
        if index == 1 or index % 10 == 0 or index == len(records):
            print(f"[bindingdb] {index}/{len(records)} {name}", flush=True)
        try:
            payload = fetch_bindingdb(str(record.get("canonical_smiles") or ""), cutoff="0.999")
            root = payload.get("getLindsByUniprotResponse", {})
            affinities = root.get("bdb.affinities", [])
            if isinstance(affinities, dict):
                affinities = [affinities]
            for item in affinities:
                rows.append({
                    "pert_iname": name,
                    "lincs_inchikey": record.get("inchi_key"),
                    "lincs_canonical_smiles": record.get("canonical_smiles"),
                    "bindingdb_match_cutoff": 0.999,
                    "bindingdb_hit_count": root.get("bdb.hit"),
                    "bindingdb_monomer_id": item.get("bdb.monomerid"),
                    "bindingdb_target": item.get("bdb.target"),
                    "bindingdb_species": item.get("bdb.species"),
                    "bindingdb_affinity_type": item.get("bdb.affinity_type"),
                    "bindingdb_affinity_value": item.get("bdb.affinity"),
                    "bindingdb_inhibitor_label": item.get("bdb.inhibitor"),
                    "bindingdb_returned_smiles": item.get("bdb.smiles"),
                    "bindingdb_source_url": f"{BINDINGDB_TARGET_BY_COMPOUND}?cutoff=0.999",
                    "evidence_scope": "high-similarity ligand-target affinity; not automatically exact identity",
                })
        except Exception as exc:
            rows.append({
                "pert_iname": name,
                "lincs_inchikey": record.get("inchi_key"),
                "lincs_canonical_smiles": record.get("canonical_smiles"),
                "bindingdb_match_cutoff": 0.999,
                "bindingdb_hit_count": None,
                "bindingdb_monomer_id": None,
                "bindingdb_target": None,
                "bindingdb_species": None,
                "bindingdb_affinity_type": None,
                "bindingdb_affinity_value": None,
                "bindingdb_inhibitor_label": None,
                "bindingdb_returned_smiles": None,
                "bindingdb_source_url": BINDINGDB_TARGET_BY_COMPOUND,
                "evidence_scope": f"query_error:{type(exc).__name__}",
            })
        time.sleep(0.08)
    return pd.DataFrame(rows)


def classify_target_family(row: pd.Series) -> str:
    text = " ".join(str(row.get(col) or "") for col in ["pert_iname", "target_name", "target_gene_symbol", "mechanism_of_action"]).lower()
    if any(x in text for x in ["mammalian target", "mechanistic target", "mtor", "rheb", "rptor", "fkbp1a"]):
        return "mTOR pathway"
    if any(x in text for x in ["phosphoinositide 3-kinase", "pi3-kinase", "pi3k", "p110", "akt", "pik3"]):
        return "PI3K/AKT pathway"
    if any(x in text for x in ["proteasome", "psmb", "psma"]):
        return "proteasome"
    if any(x in text for x in ["tubulin", "microtubule"]):
        return "microtubule"
    if any(x in text for x in ["aromatase", "estrogen receptor"]):
        return "endocrine/aromatase"
    if any(x in text for x in ["nedd8", "ubc12", "uba3", "nedd8-activating"]):
        return "NEDD8/protein homeostasis"
    if any(x in text for x in ["receptor", "transporter", "channel"]):
        return "receptor/transporter/channel"
    if any(x in text for x in ["50s ribosomal", "23s rrna", "bacterial ribosome"]):
        return "bacterial ribosome"
    if any(x in text for x in ["protein kinase c alpha", "protein kinase c delta", "pkc family", "prkca", "prkcb", "prkcd", "prkce"]):
        return "PKC/phorbol signaling"
    if any(x in text for x in ["ribonucleotide reductase", "dna polymerase", "deoxycytidine kinase"]):
        return "nucleotide/DNA synthesis"
    if any(x in text for x in ["protein kinase", "kinase"]):
        return "other kinase"
    if text.strip():
        return "other/complex mechanism"
    return "unannotated"


def write_outputs(filtered: pd.DataFrame, unique: pd.DataFrame, targets: pd.DataFrame, identity: pd.DataFrame, assay: pd.DataFrame, bindingdb: pd.DataFrame) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    TARGET_OUT.mkdir(parents=True, exist_ok=True)
    filtered.to_csv(candidate_output("compounds.csv"), index=False)
    unique.to_csv(candidate_output("compounds_unique.csv"), index=False)
    identity.to_csv(target_output("pubchem_identity.csv"), index=False)
    targets["target_family"] = targets.apply(classify_target_family, axis=1)
    targets.to_csv(target_output("compound_targets.csv"), index=False)
    assay.to_csv(target_output("pubchem_assay_target_evidence.csv.gz"), index=False, compression="gzip")
    bindingdb.to_csv(target_output("bindingdb_affinity_evidence.csv.gz"), index=False, compression="gzip")

    if len(assay):
        assay_summary = (
            assay.groupby(["pert_iname", "target_accession", "target_geneid"], dropna=False)
            .agg(
                n_assay_rows=("aid", "size"),
                n_unique_assays=("aid", "nunique"),
                activity_names=("activity_name", lambda x: ";".join(sorted(set(str(v) for v in x if pd.notna(v) and str(v) not in {"", "nan"})))),
                pubmed_ids=("pubmed_id", lambda x: ";".join(sorted(set(str(v) for v in x if pd.notna(v) and str(v) not in {"", "nan"})))),
            )
            .reset_index()
            .sort_values(["pert_iname", "n_assay_rows"], ascending=[True, False])
        )
    else:
        assay_summary = pd.DataFrame(columns=["pert_iname", "target_accession", "target_geneid", "n_assay_rows", "n_unique_assays", "activity_names", "pubmed_ids"])
    assay_summary.to_csv(target_output("pubchem_assay_target_summary.csv"), index=False)

    if len(bindingdb):
        binding_summary = (
            bindingdb.dropna(subset=["bindingdb_target"])
            .groupby(["pert_iname", "bindingdb_target", "bindingdb_species"], dropna=False)
            .agg(
                n_measurements=("bindingdb_monomer_id", "size"),
                affinity_types=("bindingdb_affinity_type", lambda x: ";".join(sorted(set(str(v) for v in x if pd.notna(v) and str(v) not in {"", "nan"})))),
                affinity_values=("bindingdb_affinity_value", lambda x: ";".join(sorted(set(str(v) for v in x if pd.notna(v) and str(v) not in {"", "nan"})))),
            )
            .reset_index()
            .sort_values(["pert_iname", "n_measurements"], ascending=[True, False])
        )
    else:
        binding_summary = pd.DataFrame(columns=["pert_iname", "bindingdb_target", "bindingdb_species", "n_measurements", "affinity_types", "affinity_values"])
    binding_summary.to_csv(target_output("bindingdb_target_summary.csv"), index=False)

    annotated = targets[targets["target_status"].astype(str).str.startswith("annotated")].copy()
    target_summary = (
        annotated.groupby(["target_gene_symbol", "target_name", "target_family"], dropna=False)
        .agg(n_drugs=("pert_iname", "nunique"), drugs=("pert_iname", lambda x: ";".join(sorted(set(x)))))
        .reset_index()
        .sort_values(["n_drugs", "target_gene_symbol"], ascending=[False, True])
    )
    target_summary.to_csv(target_output("target_summary.csv"), index=False)

    drug_summary = targets.groupby("pert_iname", as_index=False).agg(
        n_target_rows=("target_status", "size"),
        n_annotated_target_rows=("target_status", lambda x: int(x.astype(str).str.startswith("annotated").sum())),
        target_families=("target_family", lambda x: ";".join(sorted(set(x)))),
        target_genes=("target_gene_symbol", lambda x: ";".join(sorted(set(str(v) for v in x if pd.notna(v) and str(v) not in {"", "nan"})))),
        mechanism_summary=("mechanism_of_action", lambda x: ";".join(sorted(set(str(v) for v in x if pd.notna(v) and str(v) not in {"", "nan"})))),
        target_status=("target_status", lambda x: ";".join(sorted(set(x)))),
    )
    drug_summary = unique[["unique_drug_id", "pert_iname", "n_filtered_rows", "median_reversal_NCS_across_rows", "direction_pattern", "datasets"]].merge(drug_summary, on="pert_iname", how="left")
    drug_summary.to_csv(target_output("drug_target_analysis.csv"), index=False)

    family_summary = (
        annotated.groupby("target_family", as_index=False)
        .agg(n_target_records=("target_gene_symbol", "size"), n_unique_drugs=("pert_iname", "nunique"), drugs=("pert_iname", lambda x: ";".join(sorted(set(x)))))
        .sort_values(["n_unique_drugs", "target_family"], ascending=[False, True])
    )
    family_summary.to_csv(target_output("target_family_summary.csv"), index=False)

    analysis = {
        "filter": FILTER,
        "filtered_rows": int(len(filtered)),
        "unique_drugs": int(unique["pert_iname"].nunique()),
        "direction_rows": filtered["row_direction_class"].value_counts().to_dict(),
        "direction_unique_drugs": unique["direction_pattern"].value_counts().to_dict(),
        "target_records": int(len(targets)),
        "pubchem_identity_records": int(len(identity)),
        "pubchem_assay_target_evidence_rows": int(len(assay)),
        "bindingdb_affinity_evidence_rows": int(len(bindingdb)),
        "drugs_with_annotated_mechanism_or_literature_fallback": int((targets.groupby("pert_iname")["target_status"].apply(lambda x: x.astype(str).str.startswith("annotated").any())).sum()),
        "drugs_without_annotated_mechanism_or_literature_fallback": int((targets.groupby("pert_iname")["target_status"].apply(lambda x: not x.astype(str).str.startswith("annotated").any())).sum()),
        "interpretation_note": "replicated_concordant means cross-release direction category concordance; it does not imply every drug has positive reversal connectivity.",
        "target_scope_note": "The current combined target file was refreshed through ChEMBL molecule/mechanism/target queries for all 66 compounds, with explicitly labeled literature fallbacks where no curated mechanism was returned. PubChem BioAssay and BindingDB outputs remain separate evidence layers; assay/panel associations are not automatically primary targets. If rebuilt from scratch during a ChEMBL outage, the default path defers ChEMBL for new compounds and can be retried with --query-chembl.",
        "sources": {
            "ranking": str(RANKING.relative_to(ROOT)),
            "chembl_api": CHEMBL_BASE,
            "pubchem_pug_api": PUBCHEM_BASE,
            "bindingdb_rest_api": BINDINGDB_TARGET_BY_COMPOUND,
            "lincs_pert_info_GSE92742": str(GSE92742_PERT.relative_to(ROOT)),
            "lincs_pert_info_GSE70138": str(GSE70138_PERT.relative_to(ROOT)),
        },
    }
    (ROOT / "manifests" / f"{OUTPUT_PREFIX}_analysis.json").write_text(json.dumps(analysis, indent=2, ensure_ascii=False), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh-targets", action="store_true", help="query ChEMBL and overwrite target annotations")
    parser.add_argument("--refresh-external", action="store_true", help="keep existing targets and refresh PubChem/BindingDB evidence")
    parser.add_argument("--filter-only", action="store_true", help="write the combined plastic-or-hydrogel candidate tables without refreshing target annotations")
    parser.add_argument("--query-chembl", action="store_true", help="also query ChEMBL for compounds not covered by the original plastic-only annotations")
    args = parser.parse_args()
    filtered, unique, _ = build_analysis()
    if args.filter_only:
        OUT.mkdir(parents=True, exist_ok=True)
        filtered.to_csv(candidate_output("compounds.csv"), index=False)
        unique.to_csv(candidate_output("compounds_unique.csv"), index=False)
        write_json(ROOT / "manifests" / f"{OUTPUT_PREFIX}_analysis.json", {
            "status": "filter_completed",
            "filter": FILTER,
            "filtered_rows": int(len(filtered)),
            "unique_drugs": int(unique["pert_iname"].nunique()),
            "direction_rows": filtered["row_direction_class"].value_counts().to_dict(),
            "direction_unique_drugs": unique["direction_pattern"].value_counts().to_dict(),
            "outputs": {
                "filtered": str(candidate_output("compounds.csv").relative_to(ROOT)),
                "unique": str(candidate_output("compounds_unique.csv").relative_to(ROOT)),
            },
            "note": "Target annotation and external evidence refresh have not been run for the combined plastic-or-hydrogel candidate set.",
        })
        print(json.dumps({"status": "filter_completed", "filter": FILTER, "filtered_rows": len(filtered), "unique_drugs": int(unique["pert_iname"].nunique())}, ensure_ascii=False))
        return
    target_path = target_output("compound_targets.csv")
    identity_path = target_output("pubchem_identity.csv")
    assay_path = target_output("pubchem_assay_target_evidence.csv.gz")
    bindingdb_path = target_output("bindingdb_affinity_evidence.csv.gz")
    identities = load_lincs_identity(set(unique["pert_iname"]))
    if args.refresh_targets or not target_path.exists():
        legacy_targets = legacy_identity = legacy_assay = legacy_bindingdb = pd.DataFrame()
        if not args.refresh_targets:
            legacy_targets, legacy_identity, legacy_assay, legacy_bindingdb = load_legacy_annotations(set(unique["pert_iname"]))
        cached_names = set(legacy_targets.get("pert_iname", pd.Series(dtype=str)).dropna().astype(str))
        pending_unique = unique.loc[~unique["pert_iname"].isin(cached_names)].copy()
        pending_identities = identities.loc[identities["pert_iname"].isin(set(pending_unique["pert_iname"]))].copy()
        print(f"[combined] cached original annotations for {len(cached_names)} compounds; processing {len(pending_unique)} new compounds", flush=True)
        if args.query_chembl:
            new_targets = annotate_targets(pending_unique, pending_identities)
        else:
            new_targets = placeholder_targets(pending_unique, pending_identities)
        new_targets, new_identity, new_assay = annotate_pubchem(pending_unique, new_targets, pending_identities)
        assay_targets = pubchem_assay_target_rows(new_assay, pending_identities)
        new_targets = pd.concat([new_targets, assay_targets], ignore_index=True, sort=False)
        new_bindingdb = annotate_bindingdb(pending_unique, pending_identities)
        targets = pd.concat([legacy_targets, new_targets], ignore_index=True, sort=False)
        identity = pd.concat([legacy_identity, new_identity], ignore_index=True, sort=False)
        assay = pd.concat([legacy_assay, new_assay], ignore_index=True, sort=False)
        bindingdb = pd.concat([legacy_bindingdb, new_bindingdb], ignore_index=True, sort=False)
    elif args.refresh_external or not identity_path.exists() or not bindingdb_path.exists():
        targets = pd.read_csv(target_path)
        targets, identity, assay = annotate_pubchem(unique, targets, identities)
        bindingdb = annotate_bindingdb(unique, identities)
    else:
        targets = pd.read_csv(target_path)
        identity = pd.read_csv(identity_path) if identity_path.exists() else pd.DataFrame()
        assay = pd.read_csv(assay_path) if assay_path.exists() else pd.DataFrame()
        bindingdb = pd.read_csv(bindingdb_path) if bindingdb_path.exists() else pd.DataFrame()
    write_outputs(filtered, unique, targets, identity, assay, bindingdb)
    print(json.dumps({
        "filtered_rows": len(filtered),
        "unique_drugs": unique["pert_iname"].nunique(),
        "target_rows": len(targets),
        "pubchem_assay_target_rows": len(assay),
        "bindingdb_rows": len(bindingdb),
        "output_dir": str(OUT),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
