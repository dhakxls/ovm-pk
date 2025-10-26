from pathlib import Path
import os

def prepare(paths: dict, cfg: dict):
    # For MVP, just passthrough the apo PDB path.
    # Later: use PDBFixer to add hydrogens, missing atoms/residues, etc.
    return paths["apo"]
