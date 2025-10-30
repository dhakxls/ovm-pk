# src/ovmpk/docking/smina_wrapper.py
"""
Lightweight wrapper around the `smina` CLI (Vina-compatible) used by ovmpk.
Includes optional PDB cleaning, PDBFixer preprocessing, and multi-seed runs.

New:
- Auto-center docking box on Fe atom(s) of heme (HEM/HEC) if configured:
    docking:
      box:
        auto_center: heme_fe
        fallback_center: [x, y, z]
        size: [12.0, 12.0, 12.0]
"""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any, Union
import random  # For multi-seed runs

# Define a work directory for docking outputs
WORK_DIR = Path("data/work/docking")

# --- Try importing optional tools ---
try:
    import pdbfixer
    from openmm import app
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False
# ------------------------------------
__all__ = ["run", "autocenter_box_from_pdb", "SminaScorer"]

# ---------- Utilities ----------

def _need(bin_name: str) -> Optional[str]:
    """Check if an external binary is available on PATH, return path or None."""
    return shutil.which(bin_name)

def _run_subprocess(cmd: list[str], cwd: Optional[Path] = None, log_path: Optional[Path] = None) -> None:
    """Run a command, raising with full stdout/stderr on non-zero exit."""
    print(f"[cmd] Running: {' '.join(cmd)}")
    try:
        process = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            text=True,
            capture_output=True,
            check=True,  # Raise CalledProcessError on non-zero exit
        )
        if log_path:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_content = (
                f"Command: {' '.join(cmd)}\n"
                f"Return Code: {process.returncode}\n"
                f"--- stdout ---\n{process.stdout}\n"
                f"--- stderr ---\n{process.stderr}\n"
            )
            log_path.write_text(log_content)

    except FileNotFoundError:
        raise RuntimeError(f"Command not found: {cmd[0]}. Ensure it's installed and in PATH.")
    except subprocess.CalledProcessError as e:
        if log_path:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            log_content = (
                f"Command: {' '.join(e.cmd)}\n"
                f"Return Code: {e.returncode}\n"
                f"--- stdout ---\n{e.stdout}\n"
                f"--- stderr ---\n{e.stderr}\n"
            )
            log_path.write_text(log_content)
        raise RuntimeError(
            "Command failed ({rc}): {cmd}\n--- stdout ---\n{out}\n--- stderr ---\n{err}".format(
                rc=e.returncode, cmd=" ".join(e.cmd), out=e.stdout, err=e.stderr
            )
        )
    except Exception as e:
        raise RuntimeError(f"An unexpected error occurred while running command {' '.join(cmd)}: {e}")

def _has_atoms(pdbqt_path: Path) -> bool:
    """Quick check that a PDBQT actually contains atoms."""
    if not pdbqt_path.exists() or pdbqt_path.stat().st_size == 0:
        print(f"[debug] PDBQT check: File not found or empty: {pdbqt_path}")
        return False
    try:
        with pdbqt_path.open("r", errors="ignore") as fh:
            for line in fh:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    print(f"[debug] PDBQT check: Found ATOM/HETATM line in: {pdbqt_path}")
                    return True
    except Exception as e:
        print(f"[debug] PDBQT check: Error reading {pdbqt_path}: {e}")
        return False
    print(f"[debug] PDBQT check: No ATOM/HETATM lines found in: {pdbqt_path}")
    return False

def _ensure_workdir() -> Path:
    WORK_DIR.mkdir(parents=True, exist_ok=True)
    return WORK_DIR

def _alias_suffix(seed_val: Optional[int] = None) -> str:
    """Optional filename suffix (e.g., seed value)."""
    alias = os.environ.get("OVM_ALIAS", "").strip()
    parts = []
    if alias:
        parts.append(alias)
    if seed_val is not None:
        parts.append(f"s{seed_val}")
    return f"_{'_'.join(parts)}" if parts else ""

def _cfg_docking(cfg: Dict[str, Any]) -> Dict[str, Any]:
    return cfg.get("docking", {}) or {}

def _cfg_prep_protein(cfg: Dict[str, Any]) -> Dict[str, Any]:
    return cfg.get("prep", {}).get("protein", {}) or {}

# ---------- PDB helpers ----------

def autocenter_box_from_pdb(
    pdb_path: Union[str, Path],
    auto_center: Optional[str] = "heme_fe",
    fallback_center: Optional[Tuple[float, float, float]] = None,
) -> Optional[Tuple[float, float, float]]:
    """
    Compute a docking box center from a receptor PDB.

    Parameters
    ----------
    pdb_path : str | Path
        Path to PDB file.
    auto_center : {"heme_fe", None}
        Strategy for auto-centering. Currently supports "heme_fe" → centroid
        of Fe atoms that belong to heme-like residues (HEM/HEC/etc.).
        If None (or "off"), returns fallback_center.
    fallback_center : tuple[float, float, float] | None
        Fallback center to return if the chosen strategy finds nothing.

    Returns
    -------
    (x, y, z) or None
    """
    pdb_path = Path(pdb_path)

    if auto_center in (None, "off", False):
        return fallback_center

    if auto_center == "heme_fe":
        # Parse fixed-width PDB columns:
        # atom name cols 12-16; resname 17-20; x 30-38; y 38-46; z 46-54; element 76-78
        fe_xyz: list[Tuple[float, float, float]] = []
        heme_like = {"HEM", "HEC", "HEA", "HEO", "HEM3", "HMS", "HDD", "HBI"}  # permissive set

        try:
            with pdb_path.open("r", errors="ignore") as fh:
                for line in fh:
                    rec = line[0:6].strip()
                    if rec not in ("ATOM", "HETATM"):
                        continue
                    resname = line[17:20].strip().upper()
                    atom_name = line[12:16].strip().upper()
                    element = line[76:78].strip().upper()

                    # Must be Fe; prefer heme-like residues but don’t strictly require
                    if element != "FE" and atom_name != "FE":
                        continue
                    if resname and (resname not in heme_like):
                        # Allow generic Fe hits if no heme tag present
                        # (still useful for meta-hemes or variant 3-letter codes)
                        pass

                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        fe_xyz.append((x, y, z))
                    except Exception:
                        continue
        except FileNotFoundError:
            return fallback_center

        if fe_xyz:
            # centroid of all Fe atoms
            sx = sum(p[0] for p in fe_xyz)
            sy = sum(p[1] for p in fe_xyz)
            sz = sum(p[2] for p in fe_xyz)
            n = float(len(fe_xyz))
            return (sx / n, sy / n, sz / n)

        return fallback_center

    # Unknown strategy → just return fallback
    return fallback_center

# ---------- PDB Processing Steps ----------

def _clean_pdb(input_pdb: Path, output_pdb: Path, prep_cfg: Dict[str, Any]) -> Path:
    """Extracts specified chains and HETATMs from a PDB file."""
    keep_chains = set(prep_cfg.get("keep_chains", ["A"]))  # Default to chain A
    keep_hetatms = set(prep_cfg.get("keep_hetatms", ["HEM"]))  # Default to Heme
    print(f"[info] Cleaning PDB: Keeping chains {keep_chains} and HETATMs {keep_hetatms} from {input_pdb} -> {output_pdb}")

    wrote_atom = False
    in_model = False
    model_count = 0

    try:
        with input_pdb.open('r') as infile, output_pdb.open('w') as outfile:
            for line in infile:
                if line.startswith("MODEL"):
                    model_count += 1
                    if model_count == 1:
                        in_model = True
                        outfile.write(line)
                    continue
                if line.startswith("ENDMDL"):
                    if in_model:
                        outfile.write(line)
                    in_model = False
                    continue
                if model_count > 1 and not in_model:
                    continue

                record_type = line[0:6].strip()
                if record_type in ("ATOM", "ANISOU"):
                    chain_id = line[21:22].strip()
                    if chain_id in keep_chains:
                        outfile.write(line)
                        wrote_atom = True
                elif record_type == "TER":
                    if wrote_atom and line[21:22].strip() in keep_chains:
                        outfile.write(line)
                elif record_type == "HETATM":
                    res_name = line[17:20].strip()
                    if res_name in keep_hetatms:
                        outfile.write(line)
                        wrote_atom = True
                elif record_type == "END":
                    outfile.write(line)
                    break

        if not wrote_atom:
            print(f"[warn] PDB cleaning resulted in an empty file for {output_pdb}. Check keep_chains/keep_hetatms.")
        else:
            print(f"[info] Cleaned PDB saved to {output_pdb}")

        return output_pdb

    except FileNotFoundError:
        print(f"[error] Input PDB not found during cleaning: {input_pdb}")
        raise
    except Exception as e:
        print(f"[error] Failed during PDB cleaning: {e}")
        return input_pdb

def _preprocess_receptor_pdbfixer(input_pdb: Path, output_pdb: Path, prep_cfg: Dict[str, Any]) -> Path:
    """Uses PDBFixer to add hydrogens and missing atoms."""
    if not HAS_PDBFIXER:
        print("[warn] PDBFixer library not found. Skipping PDBFixer preprocessing.")
        return input_pdb

    target_ph = float(prep_cfg.get("ph", 7.4))
    print(f"[info] Running PDBFixer on {input_pdb} (pH={target_ph})...")
    try:
        fixer = pdbfixer.PDBFixer(filename=str(input_pdb))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(target_ph)

        with open(output_pdb, 'w') as f:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        print(f"[info] PDBFixer output saved to {output_pdb}")
        return output_pdb

    except Exception as e:
        print(f"[warn] PDBFixer failed: {e}. Using input PDB {input_pdb} for obabel conversion.")
        return input_pdb

def _to_pdbqt_receptor(inp_pdb: Path, out_pdbqt: Path) -> None:
    """Convert receptor PDB -> PDBQT with hydrogens using OpenBabel."""
    _need("obabel")
    print(f"[info] Converting receptor PDB to PDBQT: {inp_pdb} -> {out_pdbqt}")
    obabel_log = out_pdbqt.with_suffix(".obabel_receptor.log")
    cmd1 = ["obabel", "-ipdb", str(inp_pdb), "-opdbqt", "-O", str(out_pdbqt), "-h", "-p", "-a"]
    try:
        _run_subprocess(cmd1, log_path=obabel_log)
        if _has_atoms(out_pdbqt):
            print("[info] obabel receptor conversion successful.")
            return
    except Exception as e:
        print(f"[warn] obabel command failed with -h -p -a: {e}. Trying without -p.")

    cmd2 = ["obabel", "-ipdb", str(inp_pdb), "-opdbqt", "-O", str(out_pdbqt), "-h", "-a"]
    try:
        _run_subprocess(cmd2, log_path=obabel_log.with_suffix(".log2"))
        if _has_atoms(out_pdbqt):
            print("[info] obabel receptor conversion successful (fallback without -p).")
            return
    except Exception as e:
        print(f"[warn] obabel command failed without -p: {e}. Trying minimal flags.")

    cmd3 = ["obabel", "-ipdb", str(inp_pdb), "-opdbqt", "-O", str(out_pdbqt), "-h"]
    try:
        _run_subprocess(cmd3, log_path=obabel_log.with_suffix(".log3"))
        if _has_atoms(out_pdbqt):
            print("[info] obabel receptor conversion successful (fallback minimal -h).")
            return
    except Exception as e:
        print(f"[error] All obabel conversion attempts failed for receptor: {e}")

    if not _has_atoms(out_pdbqt):
        raise RuntimeError(
            f"Receptor PDBQT conversion failed or produced no atoms: {out_pdbqt}\n"
            f"  - Check input PDB used for conversion (cleaned/fixed): {inp_pdb}\n"
            f"  - Check obabel version and logs (e.g., {obabel_log}).\n"
            f"  - Ensure input PDB contains standard residues + HETATMs obabel can handle."
        )

def _to_pdbqt_ligand(inp_sdf: Path, out_pdbqt: Path) -> None:
    """Convert prepared ligand SDF -> PDBQT with hydrogens and charges using OpenBabel."""
    _need("obabel")
    print(f"[info] Converting ligand SDF to PDBQT: {inp_sdf} -> {out_pdbqt}")
    obabel_log = out_pdbqt.with_suffix(".obabel_ligand.log")
    cmd = ["obabel", "-isdf", str(inp_sdf), "-opdbqt", "-O", str(out_pdbqt), "-h", "--partialcharge", "gasteiger"]
    try:
        _run_subprocess(cmd, log_path=obabel_log)
    except Exception as e:
        print(f"[error] obabel ligand conversion failed: {e}")

    if not _has_atoms(out_pdbqt):
        raise RuntimeError(
            f"Ligand PDBQT conversion failed or produced no atoms: {out_pdbqt}\n"
            f"  - Check prepared ligand SDF at {inp_sdf}.\n"
            f"  - Check obabel logs (e.g., {obabel_log})."
        )
    print("[info] obabel ligand conversion successful.")

def _validate_sdf(sdf_path: Path) -> bool:
    """Thorough validation of SDF file contents."""
    if not sdf_path.exists():
        print(f"[error] SDF missing: {sdf_path}")
        return False
    if sdf_path.stat().st_size == 0:
        print(f"[error] Empty SDF: {sdf_path}")
        return False
    
    try:
        from rdkit import Chem
        suppl = Chem.SDMolSupplier(str(sdf_path))
        valid_mols = [m for m in suppl if m is not None]
        if not valid_mols:
            print(f"[error] SDF contains no valid molecules: {sdf_path}")
        return bool(valid_mols)
    except Exception as e:
        print(f"[error] SDF read failed ({sdf_path}): {str(e)}")
        return False

# ---------- Core Docking Function ----------

def run(protein_pdb_in: str | Path, ligand_sdf_in: str | Path, cfg: Dict[str, Any]) -> List[str]:
    """
    Execute smina docking, potentially multiple times with different seeds.
    Includes preprocessing steps for the receptor PDB (clean → PDBFixer → PDBQT)
    and ligand SDF → PDBQT conversion. Supports auto-centering the docking box
    from the receptor PDB (e.g., heme Fe centroid) with an explicit fallback.

    Args:
        protein_pdb_in: Path to the initial receptor PDB (already fetched/prepared upstream).
        ligand_sdf_in:  Path to the prepared ligand SDF (already protonated/charged upstream).
        cfg:             Full pipeline configuration dictionary.

    Returns:
        List[str]: Paths to all generated SDF files containing poses.
    """
    # --- Sanity / setup ---
    if _need("smina") is None:
        raise RuntimeError("Required binary 'smina' not found on PATH.")
    work = _ensure_workdir()

    protein_pdb_in = Path(protein_pdb_in)
    ligand_sdf_in  = Path(ligand_sdf_in)

    if not protein_pdb_in.exists() or protein_pdb_in.stat().st_size == 0:
        raise RuntimeError(f"Input receptor PDB not found or empty: {protein_pdb_in}")
    if not ligand_sdf_in.exists() or ligand_sdf_in.stat().st_size == 0:
        raise RuntimeError(f"Input ligand SDF not found or empty: {ligand_sdf_in}")

    prep_cfg     = _cfg_prep_protein(cfg)
    docking_cfg  = _cfg_docking(cfg)
    runtime_cfg  = cfg.get("runtime", {}) or {}
    threads      = int(runtime_cfg.get("cpu_threads", max(1, (os.cpu_count() or 2) // 2)))

    # --- Receptor preparation pipeline ---
    # 1) Clean PDB (retain chains/HETATMs per prep config)
    pdb_cleaned = work / f"{protein_pdb_in.stem}_cleaned.pdb"
    pdb_after_cleaning = _clean_pdb(protein_pdb_in, pdb_cleaned, prep_cfg)

    # 2) PDBFixer pass (optional if pdbfixer available)
    pdb_fixed = work / f"{pdb_after_cleaning.stem}_fixed.pdb"
    pdb_for_obabel = _preprocess_receptor_pdbfixer(pdb_after_cleaning, pdb_fixed, prep_cfg)

    # 3) Convert receptor PDB → PDBQT
    receptor_pdbqt = work / f"{pdb_for_obabel.stem}.pdbqt"
    _to_pdbqt_receptor(pdb_for_obabel, receptor_pdbqt)

    # --- Ligand conversion (prepared SDF → PDBQT) ---
    ligand_pdbqt = work / f"{ligand_sdf_in.stem}.pdbqt"
    _to_pdbqt_ligand(ligand_sdf_in, ligand_pdbqt)

    # --- Docking parameter resolution ---
    poses          = int(docking_cfg.get("poses", 20))
    exhaustiveness = int(docking_cfg.get("exhaustiveness", 8))
    default_seed   = int(docking_cfg.get("seed", 42))
    num_seeds      = int(docking_cfg.get("num_seeds", 1))
    blind          = bool(docking_cfg.get("blind_mode", False))
    box_cfg        = docking_cfg.get("box", {}) or {}

    # Normalize size to 3-vector if provided as scalar
    size = box_cfg.get("size")
    if isinstance(size, (int, float)):
        size = [size, size, size]

    # Resolve center for boxed docking (if not blind)
    center: Optional[Tuple[float, float, float]] = None
    if not blind:
        auto_key = box_cfg.get("auto_center")  # e.g., "heme_fe"
        fallback = box_cfg.get("fallback_center")
        if isinstance(fallback, (list, tuple)) and len(fallback) == 3:
            fallback = (float(fallback[0]), float(fallback[1]), float(fallback[2]))
        else:
            fallback = None

        # Try auto-centering first from the cleaned/fixed receptor PDB
        center = autocenter_box_from_pdb(pdb_for_obabel, auto_center=auto_key, fallback_center=fallback)

        # If auto-center failed, try explicit box.center
        if center is None:
            explicit = box_cfg.get("center")
            if isinstance(explicit, (list, tuple)) and len(explicit) == 3:
                center = (float(explicit[0]), float(explicit[1]), float(explicit[2]))

        # Validate we have both center and size
        if center is None or not (isinstance(size, (list, tuple)) and len(size) == 3):
            raise RuntimeError(
                "Boxed docking requested but a valid center/size could not be resolved.\n"
                "Provide one of:\n"
                "  - docking.box.auto_center (e.g., 'heme_fe') with optional docking.box.fallback_center and docking.box.size, or\n"
                "  - explicit docking.box.center and docking.box.size."
            )

        print(f"[info] auto-center: using center: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")

    autobox_add = float(box_cfg.get("autobox_add", 8.0))  # used only for blind mode

    # --- Execute one or many smina runs ---
    output_sdf_files: List[str] = []
    for i in range(num_seeds):
        # Choose seed: fixed if 1 run; random per run if multiple
        current_seed = default_seed if num_seeds == 1 else random.randint(1, 10_000_000)
        print(f"\n[info] Starting smina run {i+1}/{num_seeds} with seed {current_seed}...")

        run_alias = _alias_suffix(seed_val=current_seed if num_seeds > 1 else None)
        out_sdf   = work / f"{receptor_pdbqt.stem}_{ligand_pdbqt.stem}{run_alias}_poses.sdf"
        log_txt   = work / f"{receptor_pdbqt.stem}_{ligand_pdbqt.stem}{run_alias}_smina.log"

        cmd_base: List[str] = [
            "smina",
            "-r", str(receptor_pdbqt),
            "-l", str(ligand_pdbqt),
            "--out", str(out_sdf),
            "--num_modes", str(poses),
            "--exhaustiveness", str(exhaustiveness),
            "--seed", str(current_seed),
            "--cpu", str(threads),
        ]

        if blind:
            # autobox around ligand with padding
            cmd_base += ["--autobox_ligand", str(ligand_pdbqt), "--autobox_add", f"{autobox_add:.3f}"]
        else:
            # use resolved center/size
            cmd_base += [
                "--center_x", f"{center[0]:.3f}",
                "--center_y", f"{center[1]:.3f}",
                "--center_z", f"{center[2]:.3f}",
                "--size_x",   f"{float(size[0]):.3f}",
                "--size_y",   f"{float(size[1]):.3f}",
                "--size_z",   f"{float(size[2]):.3f}",
            ]

        try:
            _run_subprocess(cmd_base, log_path=log_txt)
            
            # Validate SDF output
            if not out_sdf.exists():
                print(f"[error] Output SDF missing: {out_sdf}")
                continue
            if out_sdf.stat().st_size == 0:
                print(f"[error] Empty SDF output: {out_sdf}")
                continue
                
            try:
                from rdkit import Chem
                if not any(mol is not None for mol in Chem.SDMolSupplier(str(out_sdf))):
                    print(f"[error] No valid molecules in SDF: {out_sdf}")
                    continue
            except Exception as e:
                print(f"[error] SDF validation failed: {e}")
                continue
            
            output_sdf_files.append(out_sdf.as_posix())
            print(f"[success] Valid docking output: {out_sdf}")
            
            if cfg.get("physics_module"):
                physics = physics_registry.get(cfg["physics_module"])
                if physics:
                    scores = physics.score_pose(
                        str(pdb_for_obabel), 
                        str(out_sdf)
                    )
                    # Store scores in output JSON
            
        except Exception as e:
            print(f"[error] smina run failed: {e}")
            continue

    if not output_sdf_files:
        raise RuntimeError("All smina docking runs failed to produce valid output SDF files.")

    return output_sdf_files

class SminaScorer:
    """Wrapper for smina scoring functions."""
    def __init__(self):
        self.name = "smina"
    
    def score_pose(self, protein_pdb: str, ligand_sdf: str) -> dict:
        """Default smina scoring implementation."""
        return {"total": -10.0}  # Example

__all__ = ["run", "autocenter_box_from_pdb", "SminaScorer"]
