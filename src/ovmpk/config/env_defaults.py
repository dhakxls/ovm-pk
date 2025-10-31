"""Utilities for loading default protein/ligand inputs."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Mapping
import ast

DEFAULT_ENV_PATH = Path("protein_ligand_env.md")

DEFAULT_RESOLVER_PREFERENCES: Dict[str, object] = {
    "target_species_id": 9606,
    "exp_method": "X-RAY DIFFRACTION",
    "max_resolution": 2.5,
    "results_limit": 25,
}

# Local import to avoid heavy dependencies at module import time
from ..fetchers.pdb_resolver import resolve_pdb

def load_protein_ligand_inputs(env_path: str | Path = DEFAULT_ENV_PATH) -> Dict[str, str]:
    """Parse the protein/ligand defaults file.

    The expected format is a simple key/value mapping, e.g.::

        protein = "CYP3A4"
        ligand = "ketoconazole"

    Lines beginning with ``#`` are ignored. Values are parsed with
    :func:`ast.literal_eval` to support quoted strings. If parsing fails, the
    raw string (stripped of surrounding quotes) is returned.
    """

    path = Path(env_path)
    if not path.exists():
        return {}

    defaults: Dict[str, str] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue

        key, raw_value = line.split("=", 1)
        key = key.strip()
        value = raw_value.strip()

        # Allow trailing periods from markdown sentences
        if value.endswith("."):
            value = value[:-1]

        try:
            parsed = ast.literal_eval(value)
        except Exception:
            parsed = value.strip('"').strip("'")

        if isinstance(parsed, str):
            defaults[key] = parsed
    return defaults


def ensure_defaults(data: Mapping[str, str], fallback: Mapping[str, str] | None = None) -> Dict[str, str]:
    """Merge parsed values with optional fallbacks."""
    merged: Dict[str, str] = dict(fallback or {})
    merged.update({k: v for k, v in data.items() if isinstance(v, str)})
    return merged


def resolve_default_inputs(
    env_path: str | Path = DEFAULT_ENV_PATH,
    resolver_preferences: Dict[str, object] | None = None,
) -> Dict[str, str]:
    """Return protein token, resolved PDB ID, ligand, and resolver prefs.

    The resolver preferences are merged with :data:`DEFAULT_RESOLVER_PREFERENCES`
    so callers can override resolution thresholds or required ligands while
    maintaining sensible defaults.
    """

    entries = ensure_defaults(load_protein_ligand_inputs(env_path), {})

    protein_token = entries.get("protein")
    ligand_name = entries.get("ligand")

    prefs = dict(DEFAULT_RESOLVER_PREFERENCES)
    if resolver_preferences:
        prefs.update(resolver_preferences)

    pdb_id = resolve_pdb(protein_token, prefs) if protein_token else None

    return {
        "protein_token": protein_token,
        "pdb_id": pdb_id,
        "ligand": ligand_name,
        "resolver_preferences": prefs,
    }
