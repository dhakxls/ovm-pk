# src/ovmpk/docking/smina_wrapper.py
"""
Lightweight wrapper around the `smina` CLI (Vina-compatible) used by ovmpk.
Includes optional PDB cleaning, PDBFixer preprocessing, and multi-seed runs.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional, Tuple, Dict, Any
import random # For multi-seed runs

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
        # Capture output to variables
        process = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            text=True,
            capture_output=True,
            check=True, # Raise CalledProcessError on non-zero exit
        )
        # Write stdout/stderr to log file if requested
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
        # Attempt to write logs even on failure before raising
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
                rc=e.returncode,
                cmd=" ".join(e.cmd),
                out=e.stdout,
                err=e.stderr,
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
        parts.append(f"s{seed_val}") # Use decimal seed for clarity here
    return f"_{'_'.join(parts)}" if parts else ""

def _cfg_docking(cfg: Dict[str, Any]) -> Dict[str, Any]:
    return cfg.get("docking", {}) or {}

def _cfg_prep_protein(cfg: Dict[str, Any]) -> Dict[str, Any]:
    return cfg.get("prep", {}).get("protein", {}) or {}

# ---------- PDB Processing Steps ----------

def _clean_pdb(input_pdb: Path, output_pdb: Path, prep_cfg: Dict[str, Any]) -> Path:
    """Extracts specified chains and HETATMs from a PDB file."""
    keep_chains = set(prep_cfg.get("keep_chains", ["A"])) # Default to chain A
    keep_hetatms = set(prep_cfg.get("keep_hetatms", ["HEM"])) # Default to Heme
    print(f"[info] Cleaning PDB: Keeping chains {keep_chains} and HETATMs {keep_hetatms} from {input_pdb} -> {output_pdb}")

    wrote_atom = False
    in_model = False
    model_count = 0

    try:
        with input_pdb.open('r') as infile, output_pdb.open('w') as outfile:
            for line in infile:
                # Handle MODEL records for NMR structures (keep first model only)
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
                if model_count > 1 and not in_model: # Skip lines between models > 1
                    continue

                # Process ATOM, HETATM, TER, END based on chain/residue
                record_type = line[0:6].strip()
                if record_type in ("ATOM", "ANISOU"):
                    chain_id = line[21:22].strip()
                    if chain_id in keep_chains:
                        outfile.write(line)
                        wrote_atom = True
                elif record_type == "TER":
                     # Write TER card only if the *previous* ATOM was kept (approximate)
                     # PDB format can be complex here, this is a basic heuristic
                     if wrote_atom and line[21:22].strip() in keep_chains:
                          outfile.write(line)
                     # Reset flag after TER
                     # wrote_atom = False # Keep commented - TER applies to last written atom of the chain
                elif record_type == "HETATM":
                    res_name = line[17:20].strip()
                    if res_name in keep_hetatms:
                        outfile.write(line)
                        wrote_atom = True # HETATM counts as an atom written
                elif record_type == "END":
                    outfile.write(line)
                    break # Stop processing after END
                # else: pass # Keep other records like REMARK, HEADER etc.? Maybe not needed.

        if not wrote_atom:
             print(f"[warn] PDB cleaning resulted in an empty file for {output_pdb}. Check keep_chains/keep_hetatms.")
             # Optionally delete empty file or return original path
             # output_pdb.unlink(missing_ok=True)
             # return input_pdb # Fallback if cleaning fails
        else:
             print(f"[info] Cleaned PDB saved to {output_pdb}")

        return output_pdb # Return path to cleaned file

    except FileNotFoundError:
        print(f"[error] Input PDB not found during cleaning: {input_pdb}")
        raise
    except Exception as e:
        print(f"[error] Failed during PDB cleaning: {e}")
        return input_pdb # Fallback to original on error

def _preprocess_receptor_pdbfixer(input_pdb: Path, output_pdb: Path, prep_cfg: Dict[str, Any]) -> Path:
    """Uses PDBFixer to add hydrogens and missing atoms."""
    if not HAS_PDBFIXER:
        print("[warn] PDBFixer library not found. Skipping PDBFixer preprocessing.")
        return input_pdb # Return original path if fixer not available

    target_ph = float(prep_cfg.get("ph", 7.4))
    print(f"[info] Running PDBFixer on {input_pdb} (pH={target_ph})...")
    try:
        fixer = pdbfixer.PDBFixer(filename=str(input_pdb))
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms() # Keep heterogens by default
        fixer.addMissingHydrogens(target_ph)

        with open(output_pdb, 'w') as f:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
        print(f"[info] PDBFixer output saved to {output_pdb}")
        return output_pdb # Return path to fixed file

    except Exception as e:
        print(f"[warn] PDBFixer failed: {e}. Using input PDB {input_pdb} for obabel conversion.")
        return input_pdb # Fallback to input path if fixer fails

def _to_pdbqt_receptor(inp_pdb: Path, out_pdbqt: Path) -> None:
    """Convert receptor PDB -> PDBQT with hydrogens using OpenBabel."""
    _need("obabel")
    print(f"[info] Converting receptor PDB to PDBQT: {inp_pdb} -> {out_pdbqt}")
    # Try common flags: -h (add hydrogens), -p (pH model), -a (preserve atom types for known residues)
    # OpenBabel can sometimes struggle with Heme, hence trying different flag combinations.
    obabel_log = out_pdbqt.with_suffix(".obabel_receptor.log")
    cmd1 = ["obabel", "-ipdb", str(inp_pdb), "-opdbqt", "-O", str(out_pdbqt), "-h", "-p", "-a"]
    try:
        _run_subprocess(cmd1, log_path=obabel_log)
        if _has_atoms(out_pdbqt):
            print("[info] obabel receptor conversion successful.")
            return # Success
    except Exception as e:
        print(f"[warn] obabel command failed with -h -p -a: {e}. Trying without -p.")

    # Fallback 1: Try without pH model (-p)
    cmd2 = ["obabel", "-ipdb", str(inp_pdb), "-opdbqt", "-O", str(out_pdbqt), "-h", "-a"]
    try:
        _run_subprocess(cmd2, log_path=obabel_log.with_suffix(".log2")) # Use different log
        if _has_atoms(out_pdbqt):
            print("[info] obabel receptor conversion successful (fallback without -p).")
            return # Success
    except Exception as e:
        print(f"[warn] obabel command failed without -p: {e}. Trying minimal flags.")

    # Fallback 2: Minimal flags (just add hydrogens)
    cmd3 = ["obabel", "-ipdb", str(inp_pdb), "-opdbqt", "-O", str(out_pdbqt), "-h"]
    try:
        _run_subprocess(cmd3, log_path=obabel_log.with_suffix(".log3")) # Use different log
        if _has_atoms(out_pdbqt):
            print("[info] obabel receptor conversion successful (fallback minimal -h).")
            return # Success
    except Exception as e:
        print(f"[error] All obabel conversion attempts failed for receptor: {e}")
        # Error is raised below based on _has_atoms check

    # Final check
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
    # -h adds hydrogens. --partialcharge gasteiger assigns charges (obabel should handle this automatically too)
    obabel_log = out_pdbqt.with_suffix(".obabel_ligand.log")
    cmd = ["obabel", "-isdf", str(inp_sdf), "-opdbqt", "-O", str(out_pdbqt), "-h", "--partialcharge", "gasteiger"]
    try:
         _run_subprocess(cmd, log_path=obabel_log)
    except Exception as e:
        print(f"[error] obabel ligand conversion failed: {e}")
         # Error raised below

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
    ligand_sdf_in  = Path(ligand_sdf_in) # Assume ligand_prep handled protonation/charges

    if not protein_pdb_in.exists() or protein_pdb_in.stat().st_size == 0:
        raise RuntimeError(f"Input receptor PDB not found or empty: {protein_pdb_in}")
    if not ligand_sdf_in.exists() or ligand_sdf_in.stat().st_size == 0:
        raise RuntimeError(f"Input ligand SDF not found or empty: {ligand_sdf_in}")

    prep_cfg = _cfg_prep_protein(cfg)
    docking_cfg = _cfg_docking(cfg)
    runtime_cfg = cfg.get("runtime", {})
    threads = int(runtime_cfg.get("cpu_threads", max(1, os.cpu_count() // 2)))

    # --- Receptor Preparation Pipeline ---
    # 1. Clean PDB (optional based on config or defaults)
    pdb_cleaned = work / f"{protein_pdb_in.stem}_cleaned.pdb"
    pdb_after_cleaning = _clean_pdb(protein_pdb_in, pdb_cleaned, prep_cfg)

    # 2. Preprocess with PDBFixer (optional)
    pdb_fixed = work / f"{pdb_after_cleaning.stem}_fixed.pdb"
    pdb_for_obabel = _preprocess_receptor_pdbfixer(pdb_after_cleaning, pdb_fixed, prep_cfg)

    # 3. Convert final processed PDB to PDBQT
    receptor_pdbqt = work / f"{pdb_for_obabel.stem}.pdbqt" # Base name for PDBQT
    _to_pdbqt_receptor(pdb_for_obabel, receptor_pdbqt)
    # --- End Receptor Prep ---

    # --- Ligand Preparation ---
    # Convert input ligand SDF to PDBQT (assuming ligand_prep.py handled internal prep)
    ligand_pdbqt = work / f"{ligand_sdf_in.stem}.pdbqt" # Base name
    _to_pdbqt_ligand(ligand_sdf_in, ligand_pdbqt)
    # --- End Ligand Prep ---

    # --- Docking Parameters ---
    poses = int(docking_cfg.get("poses", 20))
    exhaustiveness = int(docking_cfg.get("exhaustiveness", 8))
    default_seed = int(docking_cfg.get("seed", 42)) # Seed for single run
    num_seeds = int(docking_cfg.get("num_seeds", 1)) # Number of runs
    blind = bool(docking_cfg.get("blind_mode", False))

    output_sdf_files: List[str] = []

    # --- Run Smina (potentially multiple times) ---
    for i in range(num_seeds):
        # Determine seed for this run
        current_seed = default_seed if num_seeds == 1 else random.randint(1, 10_000_000)
        print(f"\n[info] Starting smina run {i+1}/{num_seeds} with seed {current_seed}...")

        # Generate unique filenames for this run using alias suffix
        run_alias = _alias_suffix(seed_val=current_seed if num_seeds > 1 else None)
        out_sdf = work / f"{receptor_pdbqt.stem}_{ligand_pdbqt.stem}{run_alias}_poses.sdf"
        log_txt = work / f"{receptor_pdbqt.stem}_{ligand_pdbqt.stem}{run_alias}_smina.log"

        # Build base smina command
        cmd_base: List[str] = [
            "smina",
            "-r", str(receptor_pdbqt),
            "-l", str(ligand_pdbqt),
            "--out", str(out_sdf),
            "--num_modes", str(poses),
            "--exhaustiveness", str(exhaustiveness),
            "--seed", str(current_seed),
            "--cpu", str(threads),
            # "--log", str(log_txt), # Log is now handled by _run_subprocess
        ]

        # Box handling
        if blind:
            cmd_base += ["--autobox_ligand", str(ligand_pdbqt), "--autobox_add", "8"] # Auto box around ligand + padding
        else:
            box = docking_cfg.get("box", {})
            center = box.get("center")
            size = box.get("size")
            # Ensure size is list/tuple of 3 numbers
            if isinstance(size, (int, float)): size = [size, size, size]

            if not (center and size and len(center) == 3 and len(size) == 3):
                raise RuntimeError(
                    "Boxed docking requested but 'docking.box.center' (list of 3 floats) "
                    "and/or 'docking.box.size' (list of 3 floats or single float) "
                    "are missing or malformed in config."
                )
            cmd_base += [
                "--center_x", f"{float(center[0]):.3f}",
                "--center_y", f"{float(center[1]):.3f}",
                "--center_z", f"{float(center[2]):.3f}",
                "--size_x", f"{float(size[0]):.3f}",
                "--size_y", f"{float(size[1]):.3f}",
                "--size_z", f"{float(size[2]):.3f}",
            ]

        # Execute smina command for this seed
        try:
            _run_subprocess(cmd_base, log_path=log_txt)
            # Sanity check output
            if not out_sdf.exists() or out_sdf.stat().st_size == 0:
                print(f"[warn] smina run {i+1} (seed {current_seed}) finished but produced no poses: {out_sdf}. Check log: {log_txt}")
            else:
                output_sdf_files.append(out_sdf.as_posix())
                print(f"[info] smina run {i+1} complete. Output: {out_sdf}")

        except Exception as e:
             print(f"[error] smina run {i+1} (seed {current_seed}) failed: {e}. Check log: {log_txt}")
             # Decide whether to continue other seeds or raise immediately
             # For now, let's continue to try other seeds
             continue

    if not output_sdf_files:
        raise RuntimeError("All smina docking runs failed to produce valid output SDF files.")

    return output_sdf_files

