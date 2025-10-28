from __future__ import annotations

from pathlib import Path
from typing import Optional


def ensure_pdbqt(path: Path) -> Path:
    """
    Utility to assert a receptor path *looks* like PDBQT (no conversion performed).
    Upstream code should create this with OpenBabel or other tools.
    """
    p = Path(path)
    if p.suffix.lower() != ".pdbqt":
        raise ValueError(f"Expected a .pdbqt receptor file, got: {p}")
    if not p.exists():
        raise FileNotFoundError(p)
    return p
