"""RCSB-driven protein resolver with configurable preferences."""

from __future__ import annotations

from typing import Any, Dict, List

from .protein_fetcher import GENE_TO_UNIPROT, _resolve_uniprot, _select_best_entry_id


DEFAULT_PREFERENCES: Dict[str, Any] = {
    "target_species_id": 9606,  # Homo sapiens
    "exp_method": "X-RAY DIFFRACTION",
    "max_resolution": 2.5,
    "required_ligands": [],
    "results_limit": 25,
    "gene_to_uniprot_map": {},
    "fallback_map": {},
}


def resolve_pdb(protein_input: str, preferences: Dict[str, Any] | None = None) -> str:
    """Resolve a protein identifier to the best-fit PDB ID using RCSB heuristics."""

    if not isinstance(protein_input, str):
        raise TypeError("protein_input must be a string")

    token = protein_input.strip()
    if len(token) == 4 and token.isalnum():
        return token.upper()

    prefs = {**DEFAULT_PREFERENCES, **(preferences or {})}
    species_id = int(prefs.get("target_species_id", 9606))
    method = str(prefs.get("exp_method", "X-RAY DIFFRACTION"))
    max_resolution = float(prefs.get("max_resolution", 2.5))
    required_ligands: List[str] = list(prefs.get("required_ligands", []))
    results_limit = int(prefs.get("results_limit", 25))

    gene_map = dict(GENE_TO_UNIPROT)
    gene_map.update({k.upper(): v for k, v in prefs.get("gene_to_uniprot_map", {}).items()})

    search_term = _resolve_uniprot(token, gene_map) or token

    pdb_id = _select_best_entry_id(
        search_term,
        species_id,
        method,
        max_resolution,
        required_ligands,
        results_limit,
    )
    if pdb_id:
        return pdb_id

    fallback = prefs.get("fallback_map", {}) or {}
    fallback_upper = {k.upper(): v for k, v in fallback.items()}
    return fallback_upper.get(token.upper(), token.upper() if len(token) == 4 else token)
