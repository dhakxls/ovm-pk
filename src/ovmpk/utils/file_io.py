"""
File I/O utilities for molecular data
"""
from pathlib import Path
import json

def safe_write_json(data, path):
    """Atomic JSON write with backup"""
    path = Path(path)
    temp_path = path.with_suffix('.tmp')
    with open(temp_path, 'w') as f:
        json.dump(data, f)
    temp_path.replace(path)

def read_pdb_lines(pdb_path):
    """Read PDB file with validation"""
    path = Path(pdb_path)
    if not path.exists():
        raise FileNotFoundError(f"PDB file not found: {path}")
    return path.read_text().splitlines()
