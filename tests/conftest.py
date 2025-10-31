"""Pytest configuration for custom arguments sourced from env defaults."""

from pathlib import Path

from src.ovmpk.config.env_defaults import ensure_defaults, load_protein_ligand_inputs


def pytest_addoption(parser):
    defaults = ensure_defaults(
        load_protein_ligand_inputs(Path("protein_ligand_env.md")),
        {"protein": "", "ligand": ""},
    )

    parser.addoption(
        "--protein",
        action="store",
        default=defaults.get("protein", ""),
        help="PDB ID or protein name (defaults to value in protein_ligand_env.md)",
    )
    parser.addoption(
        "--ligand",
        action="store",
        default=defaults.get("ligand", ""),
        help="Ligand name (defaults to value in protein_ligand_env.md)",
    )
