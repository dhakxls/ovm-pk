"""
Robust metal parameterization with validation
Combines best of metal.py and metal_params.py
"""
from openmm import app
from openmm import System  # Correct import source
from pathlib import Path
from typing import Dict

METAL_ELEMENTS = {'FE', 'ZN', 'MG', 'CA', 'MN', 'CO', 'NI', 'CU'}

def validate_pdb(pdb_file: Path):
    """Verify PDB contains metal atoms"""
    with open(pdb_file) as f:
        content = f.read()
    if not any(f" {el} " in content.upper() for el in METAL_ELEMENTS):
        raise ValueError(f"No supported metal elements found in {pdb_file}")

def parameterize_metal(pdb_file: str, xml_template: str) -> System:
    """Parameterize metal-containing system with validation"""
    validate_pdb(Path(pdb_file))
    pdb = app.PDBFile(pdb_file)
    forcefield = app.ForceField('amber14-all.xml', xml_template)
    return forcefield.createSystem(pdb.topology)
