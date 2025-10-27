#!/usr/bin/env python3
"""
Protein fetcher (rcsb-api only)

- No fallbacks to requests/HTTP.
- Uses rcsb-api Search API to find an entry (PDB ID), preferring UniProt-anchored hits.
- Downloads the full structure as mmCIF via rcsb-api's Model Server client.
- Keeps the public entrypoint: fetch(identifier: str, cfg: Dict[str, Any]) -> Dict[str, Path]
  so test_fetchers.py can validate the output file (checks for file existence & size).
"""

import re
from pathlib import Path
from typing import Dict, Any, Optional, List

# rcsb-api (install: pip install rcsb-api)
from rcsbapi.search import TextQuery, Sort, AttributeQuery, NestedAttributeQuery
from rcsbapi.search import search_attributes as attrs
from rcsbapi.model import ModelQuery

# -----------------------------
# Gene -> UniProt map (extend via cfg)
# -----------------------------
GENE_TO_UNIPROT = {
    "CYP3A4": "P08684",
}

# -----------------------------
# UniProt pattern (canonical or long accessions; optional isoform suffix)
# Examples: P08684, Q9Y6K9, A0A024RBG1, Q1HVC0-2
# -----------------------------
UNIPROT_REGEX = re.compile(
    r"^((?:[OPQ][0-9][A-Z0-9]{3}[0-9])|(?:[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}))(?:-[0-9]+)?$",
    flags=re.IGNORECASE,
)


def _looks_like_uniprot(s: str) -> bool:
    return bool(UNIPROT_REGEX.match(s.strip()))


def _resolve_uniprot(identifier: str, mapping: Dict[str, str]) -> Optional[str]:
    """Resolve UniProt from a gene symbol or accept a UniProt-looking string."""
    if not isinstance(identifier, str):
        print("[warn] Identifier is not a string.")
        return None
    ident = identifier.strip().upper()
    if ident in mapping:
        up = mapping[ident]
        print(f"[info] Mapped gene '{identifier}' to UniProt '{up}'.")
        return up
    if _looks_like_uniprot(ident):
        print(f"[info] Treating '{identifier}' as a UniProt accession.")
        return ident
    print(f"[warn] No UniProt for '{identifier}'. Falling back to text query.")
    return None


def _build_query(
    term: str,
    species_id: int,
    method: str,
    max_resolution: float,
    required_ligands: List[str],
    rows: int = 100,
):
    """Build and execute an RCSB query using rcsb-api; returns an iterator of PDB IDs."""
    method_u = method.upper()

    # Identifier: UniProt mapping via NESTED attribute pair (name + accession), else full-text.
    if _looks_like_uniprot(term):
        q = NestedAttributeQuery(
            AttributeQuery(
                "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name",
                "exact_match",
                "UniProt",
            ),
            AttributeQuery(
                "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "exact_match",
                term.upper(),
            ),
        )
    else:
        q = TextQuery(value=term)

    # Constrain by experimental method always; by species only for non-UniProt text queries
    q &= (attrs.exptl.method == method_u)
    if not _looks_like_uniprot(term):
        q &= (attrs.rcsb_entity_source_organism.ncbi_taxonomy_id == int(species_id))

    # Resolution filter + sorting (for X-ray/EM)
    sort = None
    if ("X-RAY" in method_u) or ("ELECTRON MICROSCOPY" in method_u):
        q &= (attrs.rcsb_entry_info.resolution_combined <= float(max_resolution))
        sort = Sort(sort_by="rcsb_entry_info.resolution_combined", direction="asc")

    # Required ligands (AND each comp_id)
    for comp in (required_ligands or []):
        q &= (attrs.rcsb_nonpolymer_entity_instance_container_identifiers.comp_id == comp.upper())

    # Execute: default return type is entry (PDB IDs)
    iterator = q(rows=rows, sort=sort) if sort else q(rows=rows)
    return iterator


def _select_best_entry_id(
    search_term: str,
    species_id: int,
    method: str,
    max_resolution: float,
    required_ligands: List[str],
    rows: int,
) -> Optional[str]:
    it = _build_query(search_term, species_id, method, max_resolution, required_ligands, rows=rows)
    first = next(iter(it), None)
    if first:
        print(f"[info] Selected PDB ID '{first}'.")
    else:
        print("[warn] No PDB IDs matched the criteria.")
    return first


def fetch(identifier: str, cfg: Dict[str, Any]) -> Dict[str, Path]:
    """Public API used by test_fetchers.py.

    Config keys under fetch.protein:
      - target_species_id (int, default 9606)
      - exp_method (str, default 'X-RAY DIFFRACTION')
      - max_resolution (float, default 3.0)
      - required_ligands (List[str], default [])
      - gene_to_uniprot_map (Dict[str,str], optional)
      - results_limit (int, default 100)

    Returns: {"best_match": Path}
    """
    pcfg = cfg.get("fetch", {}).get("protein", {})

    species_id = int(pcfg.get("target_species_id", 9606))
    method = str(pcfg.get("exp_method", "X-RAY DIFFRACTION"))
    max_resolution = float(pcfg.get("max_resolution", 3.0))
    required_ligands = list(pcfg.get("required_ligands", []))
    gene_map = dict(pcfg.get("gene_to_uniprot_map", GENE_TO_UNIPROT))
    rows = int(pcfg.get("results_limit", 100))

    # Resolve UniProt if we can; otherwise search on the text
    uniprot = _resolve_uniprot(identifier, gene_map)
    search_term = uniprot or identifier

    # Output target
    base = Path("data/input/proteins")
    base.mkdir(parents=True, exist_ok=True)
    id_part = (uniprot or identifier.replace(" ", "_")).upper()
    lig_suffix = f"_req_{'_'.join([l.upper() for l in required_ligands])}" if required_ligands else "_any"
    filename_base = f"{id_part}_{species_id}{lig_suffix}"

    print("[info] Querying RCSB PDB via rcsb-api…")
    pdb_id = _select_best_entry_id(search_term, species_id, method, max_resolution, required_ligands, rows)
    if not pdb_id:
        raise FileNotFoundError(
            f"No suitable PDB entry for '{identifier}' (Species={species_id}, Method='{method}', MaxRes={max_resolution}, Ligands={required_ligands})."
        )

    outp = base / f"{filename_base}_{pdb_id}.cif"

    if outp.exists():
        print(f"[info] File already exists: {outp}")
        return {"best_match": outp}

    print(f"[info] Downloading full structure for {pdb_id} → {outp} …")
    try:
        mq = ModelQuery()
        saved = mq.get_full_structure(
            entry_id=pdb_id,
            encoding="cif",
            download=True,
            filename=outp.name,
            file_directory=str(outp.parent),
        )
        out_path = Path(saved)
        # Sanity check
        if not out_path.exists() or out_path.stat().st_size < 100:
            raise RuntimeError("Downloaded file looks empty or missing.")
        print(f"[info] Saved: {out_path}")
        return {"best_match": out_path}
    except Exception as e:
        outp.unlink(missing_ok=True)
        raise RuntimeError(f"Failed to download mmCIF for {pdb_id}: {e}") from e
