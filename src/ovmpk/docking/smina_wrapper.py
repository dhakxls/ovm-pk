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
from typing import List, Optional, Tuple, Dict, Any
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

def _autocenter_box_from_pdb(pdb_path: Path) -> Optional[Tuple[float, float, float]]:
    """
    Extract the (x,y,z) of Fe atom(s) from HEM/HEC residues in a PDB.
    Returns the average Fe position if multiple are found, else None.
    """
    if not pdb_path.exists():
        print(f"[warn] auto-center: PDB not found: {pdb_path}")
        return None

    centers: list[Tuple[float, float, float]] = []
    try:
        with pdb_path.open("r", errors="ignore") as fh:
            for line in fh:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                # PDB columns (1-based): resname 18-20, atom name 13-16, x 31-38, y 39-46, z 47-54, element 77-78
                resn = line[17:20].strip().upper()
                if resn not in {"HEM", "HEC"}:
                    continue
                atom_name = line[12:16].strip().upper()
                element = line[76:78].strip().upper() if len(line) >= 78 else ""
                if atom_name == "FE" or element == "FE":
                    try:
                        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                        centers.append((x, y, z))
                    except ValueError:
                        continue
    except Exception as e:
        print(f"[warn] auto-center: error scanning PDB {pdb_path}: {e}")
        return None

    if not centers:
        print("[info] auto-center: no HEM/HEC Fe atom found; will use fallback if provided.")
        return None

    n = len(centers)
    cx = sum(x for x, _, _ in centers) / n
    cy = sum(y for _, y, _ in centers) / n
    cz = sum(z for _, _, z in centers) / n
    print(f"[info] auto-center: using mean Fe position from {n} heme(s): ({cx:.3f}, {cy:.3f}, {cz:.3f})")
    return (cx, cy, cz)

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

# ---------- Core Docking Function ----------

def run(protein_pdb_in: str | Path, ligand_sdf_in: str | Path, cfg: Dict[str, Any]) -> List[str]:
    """
    Execute smina docking, potentially multiple times with different seeds.
    Includes preprocessing steps for the receptor PDB.

    Args:
        protein_pdb_in: Path to the initial (fetched) receptor PDB.
        ligand_sdf_in: Path to the initial (fetched/prepared) ligand SDF.
        cfg: Configuration dictionary.

    Returns:
        List[str]: Paths to all generated SDF files containing poses.
    """
    _need("smina")
    work = _ensure_workdir()

    protein_pdb_in = Path(protein_pdb_in)
    ligand_sdf_in = Path(ligand_sdf_in)  # Assume ligand_prep handled protonation/charges

    if not protein_pdb_in.exists() or protein_pdb_in.stat().st_size == 0:
        raise RuntimeError(f"Input receptor PDB not found or empty: {protein_pdb_in}")
    if not ligand_sdf_in.exists() or ligand_sdf_in.stat().st_size == 0:
        raise RuntimeError(f"Input ligand SDF not found or empty: {ligand_sdf_in}")

    prep_cfg = _cfg_prep_protein(cfg)
    docking_cfg = _cfg_docking(cfg)
    runtime_cfg = cfg.get("runtime", {})
    threads = int(runtime_cfg.get("cpu_threads", max(1, os.cpu_count() // 2)))

    # --- Receptor Preparation Pipeline ---
    pdb_cleaned = work / f"{protein_pdb_in.stem}_cleaned.pdb"
    pdb_after_cleaning = _clean_pdb(protein_pdb_in, pdb_cleaned, prep_cfg)

    pdb_fixed = work / f"{pdb_after_cleaning.stem}_fixed.pdb"
    pdb_for_obabel = _preprocess_receptor_pdbfixer(pdb_after_cleaning, pdb_fixed, prep_cfg)

    receptor_pdbqt = work / f"{pdb_for_obabel.stem}.pdbqt"
    _to_pdbqt_receptor(pdb_for_obabel, receptor_pdbqt)
    # --- End Receptor Prep ---

    # --- Ligand Preparation ---
    ligand_pdbqt = work / f"{ligand_sdf_in.stem}.pdbqt"
    _to_pdbqt_ligand(ligand_sdf_in, ligand_pdbqt)
    # --- End Ligand Prep ---

    # --- Docking Parameters ---
    poses = int(docking_cfg.get("poses", 20))
    exhaustiveness = int(docking_cfg.get("exhaustiveness", 8))
    default_seed = int(docking_cfg.get("seed", 42))
    num_seeds = int(docking_cfg.get("num_seeds", 1))
    blind = bool(docking_cfg.get("blind_mode", False))

    output_sdf_files: List[str] = []

    # Precompute box center/size if not blind
    box_center: Optional[Tuple[float, float, float]] = None
    box_size = None
    if not blind:
        box = docking_cfg.get("box", {}) or {}
        # size must be present
        size = box.get("size")
        if isinstance(size, (int, float)):
            size = [size, size, size]
        if not (isinstance(size, (list, tuple)) and len(size) == 3):
            raise RuntimeError("Docking box 'size' must be a list/tuple of 3 floats.")
        box_size = [float(size[0]), float(size[1]), float(size[2])]

        center = box.get("center")
        auto_mode = (box.get("auto_center") or "").lower()
        fallback_center = box.get("fallback_center")

        if auto_mode == "heme_fe":
            auto = _autocenter_box_from_pdb(pdb_for_obabel)
            if auto is not None:
                box_center = auto
            elif fallback_center and len(fallback_center) == 3:
                box_center = (float(fallback_center[0]), float(fallback_center[1]), float(fallback_center[2]))
                print(f"[info] auto-center: using fallback center: {box_center}")
            elif center and len(center) == 3:
                box_center = (float(center[0]), float(center[1]), float(center[2]))
                print(f"[info] auto-center: using user-provided 'center' as fallback: {box_center}")
            else:
                raise RuntimeError(
                    "auto_center=heme_fe requested but no Fe found and neither "
                    "'fallback_center' nor 'center' provided."
                )
        else:
            # No auto centering; require explicit center
            if not (center and len(center) == 3):
                raise RuntimeError("Boxed docking requires 'box.center' when auto_center is not used.")
            box_center = (float(center[0]), float(center[1]), float(center[2]))

    # --- Run Smina (potentially multiple times) ---
    for i in range(num_seeds):
        current_seed = default_seed if num_seeds == 1 else random.randint(1, 10_000_000)
        print(f"\n[info] Starting smina run {i+1}/{num_seeds} with seed {current_seed}...")

        run_alias = _alias_suffix(seed_val=current_seed if num_seeds > 1 else None)
        out_sdf = work / f"{receptor_pdbqt.stem}_{ligand_pdbqt.stem}{run_alias}_poses.sdf"
        log_txt = work / f"{receptor_pdbqt.stem}_{ligand_pdbqt.stem}{run_alias}_smina.log"

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
            cmd_base += ["--autobox_ligand", str(ligand_pdbqt), "--autobox_add", "8"]
        else:
            cx, cy, cz = box_center  # type: ignore
            sx, sy, sz = box_size    # type: ignore
            cmd_base += [
                "--center_x", f"{cx:.3f}",
                "--center_y", f"{cy:.3f}",
                "--center_z", f"{cz:.3f}",
                "--size_x", f"{sx:.3f}",
                "--size_y", f"{sy:.3f}",
                "--size_z", f"{sz:.3f}",
            ]

        try:
            _run_subprocess(cmd_base, log_path=log_txt)
            if not out_sdf.exists() or out_sdf.stat().st_size == 0:
                print(f"[warn] smina run {i+1} (seed {current_seed}) finished but produced no poses: {out_sdf}. Check log: {log_txt}")
            else:
                output_sdf_files.append(out_sdf.as_posix())
                print(f"[info] smina run {i+1} complete. Output: {out_sdf}")

        except Exception as e:
            print(f"[error] smina run {i+1} (seed {current_seed}) failed: {e}. Check log: {log_txt}")
            # Continue other seeds; change to `raise` if you want fail-fast across seeds
            continue

    if not output_sdf_files:
        raise RuntimeError("All smina docking runs failed to produce valid output SDF files.")

    return output_sdf_files
