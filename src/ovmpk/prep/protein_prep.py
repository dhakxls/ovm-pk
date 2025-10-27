from pathlib import Path
import os
from typing import Dict, Any

# Define a work directory for protein prep outputs
WORK_DIR = Path("data/work/protein_prep")

# Try importing PDBFixer and OpenMM
try:
    import pdbfixer
    from openmm import app
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False

def prepare(paths: Dict[str, Path], cfg: Dict[str, Any]) -> Path:
    """
    Prepares a protein PDB file for docking or simulation using PDBFixer.

    Steps (if PDBFixer is available):
    1. Reads the input PDB file (typically the 'apo' structure).
    2. Finds missing residues and atoms.
    3. Adds missing heavy atoms.
    4. Adds missing hydrogens based on a specified pH.
    5. Importantly, keeps heterogens (like Heme) in the structure.
    6. Writes the processed structure to a new PDB file in the work directory.

    Args:
        paths: Dictionary containing input paths, expects "apo" key.
        cfg: Configuration dictionary, expects a 'prep.protein' section.

    Returns:
        Path to the prepared PDB file (either the fixed one or the original).
    """
    input_pdb_path = paths.get("apo")
    if input_pdb_path is None or not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB file not found in paths dictionary or path invalid: {input_pdb_path}")

    # Get config parameters
    prep_cfg = cfg.get("prep", {}).get("protein", {})
    target_ph = float(prep_cfg.get("ph", 7.4))
    output_suffix = prep_cfg.get("output_suffix", f"_fixed_ph{target_ph}")
    run_fixer = prep_cfg.get("run_pdbfixer", True) # Option to disable fixer

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    outp_pdb = WORK_DIR / f"{input_pdb_path.stem}{output_suffix}.pdb"

    if HAS_PDBFIXER and run_fixer:
        print(f"[info] Running PDBFixer on {input_pdb_path} (target pH: {target_ph})...")
        try:
            fixer = pdbfixer.PDBFixer(filename=str(input_pdb_path))

            # Find missing elements but keep heterogens like Heme
            fixer.findMissingResidues()
            fixer.findNonstandardResidues() # Identify non-standard ones
            fixer.findMissingAtoms()

            # Add missing heavy atoms, DO NOT remove heterogens
            fixer.addMissingAtoms()

            # Add missing hydrogens at the target pH
            fixer.addMissingHydrogens(target_ph)

            # Write the fixed PDB file
            with open(outp_pdb, 'w') as f:
                app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True) # keepIds helps maintain residue/atom numbering

            print(f"[info] PDBFixer complete. Output: {outp_pdb}")
            return outp_pdb

        except Exception as e:
            print(f"[warn] PDBFixer failed: {e}. Returning original PDB path: {input_pdb_path}")
            # Fallback to original PDB if fixer fails
            return input_pdb_path
    else:
        if not run_fixer:
            print("[info] PDBFixer step explicitly disabled in config.")
        else:
            print("[warn] PDBFixer library not found. Skipping protein preparation/fixing step.")
        # Return the original PDB path if fixer isn't run
        return input_pdb_path
