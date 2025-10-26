from pathlib import Path
import os

def prepare(sdf_path: Path, cfg: dict):
    # For MVP: passthrough SDF.
    # Later: RDKit + protonation/charges/tautomer enumeration.
    return sdf_path
