# src/ovmpk/docking/smina_wrapper.py
"""
Lightweight wrapper around the `smina` CLI (Vina-compatible) used by ovmpk.

Responsibilities
- Convert receptor PDB and ligand SDF to PDBQT via OpenBabel (`obabel`)
- Build and run a smina command line (boxed or blind)
- Write outputs (poses SDF + smina log) into data/work/docking
- Make filenames stable and suffix with an optional alias (OVM_ALIAS env)
- Fail fast with human-friendly diagnostics when inputs/tools are missing

Return
- A list containing the path to the SDF file with poses (as strings)

Notes
- This module is CPU-only (smina/vina do not use GPUs). Use --cpu for threads.
- Requires external binaries: `smina`, `obabel`
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any


WORK_DIR = Path("data/work/docking")


# ---------- Utilities ----------

def _need(bin_name: str) -> None:
    """Ensure an external binary is available on PATH."""
    if shutil.which(bin_name) is None:
        raise RuntimeError(
            f"Required binary '{bin_name}' not found on PATH.\n"
            f"  - Install it and ensure it's discoverable (e.g., conda install -c conda-forge {bin_name})."
        )


def _run(cmd: List[str], cwd: Optional[Path] = None) -> None:
    """Run a command, raising with full stdout/stderr on non-zero exit."""
    p = subprocess.run(
        cmd,
        cwd=str(cwd) if cwd else None,
        text=True,
        capture_output=True,
    )
    if p.returncode != 0:
        raise RuntimeError(
            "Command failed ({rc}): {cmd}\n--- stdout ---\n{out}\n--- stderr ---\n{err}".format(
                rc=p.returncode,
                cmd=" ".join(cmd),
                out=p.stdout,
                err=p.stderr,
            )
        )


def _has_atoms(pdbqt_path: Path) -> bool:
    """Quick check that a PDBQT actually contains atoms."""
    if not pdbqt_path.exists() or pdbqt_path.stat().st_size == 0:
        return False
    try:
        with pdbqt_path.open("r", errors="ignore") as fh:
            for line in fh:
                # PDBQT uses ATOM/HETATM records (columns are similar to PDB)
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    return True
    except Exception:
        return False
    return False


def _ensure_workdir() -> Path:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    return WORK_DIR


def _alias_suffix() -> str:
    """Optional filename suffix to avoid collisions (e.g., seed name)."""
    # e.g., export OVM_ALIAS="s42" to produce *_s42_*.*
    alias = os.environ.get("OVM_ALIAS", "").strip()
    return f"_{alias}" if alias else ""


def _cfg_docking(cfg: Dict[str, Any]) -> Dict[str, Any]:
    d = cfg.get("docking", {}) or {}
    return d


def _to_pdbqt_receptor(inp: Path, out: Path) -> None:
    """
    Convert receptor PDB -> PDBQT with hydrogens using OpenBabel.
    """
    _need("obabel")
    # Add hydrogens (-h); keep as-is otherwise. You may add options here if needed (e.g. -xr for waters).
    cmd = ["obabel", "-ipdb", str(inp), "-opdbqt", "-O", str(out), "-h"]
    _run(cmd)
    if not _has_atoms(out):
        raise RuntimeError(
            f"Receptor PDBQT has no atoms: {out}\n"
            f"  - Check receptor PDB at {inp}\n"
            f"  - Ensure it contains heavy atoms (and not an empty/stripped file)"
        )


def _to_pdbqt_ligand(inp: Path, out: Path) -> None:
    """
    Convert ligand SDF -> PDBQT with hydrogens using OpenBabel.
    """
    _need("obabel")
    # -h adds hydrogens. OpenBabel will set atom types and PDBQT charges as needed.
    cmd = ["obabel", "-isdf", str(inp), "-opdbqt", "-O", str(out), "-h"]
    _run(cmd)
    if not _has_atoms(out):
        raise RuntimeError(
            f"Ligand PDBQT has no atoms: {out}\n"
            f"  - Check ligand SDF at {inp} (is it valid SDF with 3D coords?)\n"
            f"  - If you generated protonation variants, ensure the file is non-empty."
        )


# ---------- Core ----------

def run(protein_pdb: str | Path, ligand_sdf: str | Path, cfg: Dict[str, Any]) -> List[str]:
    """
    Execute a smina docking run and return the SDF file(s) with poses.
    Parameters
    ----------
    protein_pdb : path to receptor PDB
    ligand_sdf  : path to ligand SDF
    cfg         : configuration dict (expects 'docking', 'runtime', 'reporting' entries)
    Returns
    -------
    List[str] : [path_to_poses_sdf]
    """
    _need("smina")
    _ensure_workdir()

    protein_pdb = Path(protein_pdb)
    ligand_sdf  = Path(ligand_sdf)

    if not protein_pdb.exists() or protein_pdb.stat().st_size == 0:
        raise RuntimeError(f"Receptor PDB not found or empty: {protein_pdb}")
    if not ligand_sdf.exists() or ligand_sdf.stat().st_size == 0:
        raise RuntimeError(f"Ligand SDF not found or empty: {ligand_sdf}")

    docking_cfg = _cfg_docking(cfg)
    runtime_cfg = (cfg.get("runtime") or {})
    threads     = int(runtime_cfg.get("cpu_threads", max(1, (os.cpu_count() or 8) // 2)))

    alias = _alias_suffix()

    # File layout
    work = _ensure_workdir()
    receptor_pdbqt = work / f"{protein_pdb.stem}{alias}.pdbqt"
    ligand_pdbqt   = work / f"{ligand_sdf.stem}{alias}.pdbqt"
    out_sdf        = work / f"{protein_pdb.stem}_{ligand_sdf.stem}{alias}_poses.sdf"
    log_txt        = work / f"{protein_pdb.stem}_{ligand_sdf.stem}{alias}_smina.log"

    # Convert inputs
    _to_pdbqt_receptor(protein_pdb, receptor_pdbqt)
    _to_pdbqt_ligand(ligand_sdf, ligand_pdbqt)

    # Schedules / params
    poses           = int(docking_cfg.get("poses", 20))
    exhaustiveness  = int(docking_cfg.get("exhaustiveness", 8))
    seed            = int(docking_cfg.get("seed", 42))
    blind           = bool(docking_cfg.get("blind_mode", False))

    # Build smina command
    cmd: List[str] = [
        "smina",
        "-r", str(receptor_pdbqt),
        "-l", str(ligand_pdbqt),
        "--out", str(out_sdf),
        "--num_modes", str(poses),
        "--exhaustiveness", str(exhaustiveness),
        "--seed", str(seed),
        "--cpu", str(threads),
        "--log", str(log_txt),
    ]

    # Box handling
    if blind:
        # Let smina auto-box on the receptor; you can make this tighter with --autobox_add if desired
        cmd += ["--autobox", str(receptor_pdbqt)]
    else:
        box = (docking_cfg.get("box") or {})
        center = box.get("center")
        size   = box.get("size")
        if not (center and size and len(center) == 3 and len(size) == 3):
            raise RuntimeError(
                "Boxed docking requested but 'docking.box.center' and 'docking.box.size' "
                "are missing or malformed in config."
            )
        cx, cy, cz = [float(center[0]), float(center[1]), float(center[2])]
        sx, sy, sz = [float(size[0]), float(size[1]), float(size[2])]
        cmd += [
            "--center_x", f"{cx:.3f}",
            "--center_y", f"{cy:.3f}",
            "--center_z", f"{cz:.3f}",
            "--size_x", f"{sx:.3f}",
            "--size_y", f"{sy:.3f}",
            "--size_z", f"{sz:.3f}",
        ]

    # Run smina
    _run(cmd)

    # Sanity: ensure we produced something
    if not out_sdf.exists() or out_sdf.stat().st_size == 0:
        # If smina succeeded but wrote 0 poses (rare), still surface clearly
        raise RuntimeError(
            f"smina finished but produced no poses: {out_sdf}\n"
            f"Check log: {log_txt}"
        )

    return [out_sdf.as_posix()]
