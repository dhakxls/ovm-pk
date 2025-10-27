from pathlib import Path
import os
import subprocess
import shutil
from typing import Dict, Any, Optional

# RDKit imports
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# Define a work directory for ligand prep outputs
WORK_DIR = Path("data/work/ligand_prep")

def _need(bin_name: str) -> Optional[str]:
    """Check if an external binary is available on PATH, return path or None."""
    return shutil.which(bin_name)

def _run_subprocess(cmd: list[str], cwd: Optional[Path] = None) -> None:
    """Run a command, raising with full stdout/stderr on non-zero exit."""
    print(f"[cmd] Running: {' '.join(cmd)}")
    try:
        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            text=True,
            capture_output=True,
            check=True, # Raise CalledProcessError on non-zero exit
        )
        # Optional: print stdout/stderr for debugging
        # print(f"--- stdout ---\n{p.stdout}\n--- stderr ---\n{p.stderr}")
    except FileNotFoundError:
         raise RuntimeError(f"Command not found: {cmd[0]}. Ensure it's installed and in PATH.")
    except subprocess.CalledProcessError as e:
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


def prepare(sdf_path: Path, cfg: Dict[str, Any]) -> Path:
    """
    Prepares a ligand SDF file for docking.

    Steps:
    1. Protonates the ligand at a specified pH using dimorphite-dl (if available).
    2. Calculates partial charges (Gasteiger default) using RDKit.
    3. Generates 3D coordinates if needed (basic embed).
    4. Writes the processed molecule to a new SDF file in the work directory.

    Args:
        sdf_path: Path to the input ligand SDF file.
        cfg: Configuration dictionary, expects a 'prep.ligand' section.

    Returns:
        Path to the prepared SDF file.
    """
    if not HAS_RDKIT:
        raise RuntimeError("RDKit is required for ligand preparation but not found. Please install it.")

    # Get config parameters
    prep_cfg = cfg.get("prep", {}).get("ligand", {})
    target_ph = float(prep_cfg.get("ph", 7.4))
    charge_method = prep_cfg.get("charge_method", "gasteiger") # Currently only gasteiger via RDKit is implemented here
    output_suffix = prep_cfg.get("output_suffix", f"_prepared_ph{target_ph}")
    force_3d = prep_cfg.get("force_3d", True) # Attempt to generate 3D if input lacks it

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    outp_sdf = WORK_DIR / f"{sdf_path.stem}{output_suffix}.sdf"

    print(f"[info] Starting ligand preparation for {sdf_path}...")

    # --- Step 1: Protonation using dimorphite-dl ---
    dimorphite_bin = _need("dimorphite_dl")
    temp_protonated_sdf = WORK_DIR / f"{sdf_path.stem}_temp_protonated.sdf"

    if dimorphite_bin:
        print(f"[info] Running dimorphite-dl for pH {target_ph}...")
        dimorphite_cmd = [
            dimorphite_bin,
            "--input_file", str(sdf_path),
            "--output_file", str(temp_protonated_sdf),
            "--min_ph", str(target_ph),
            "--max_ph", str(target_ph),
            "--pka_precision", "1", # Use fewer states for simplicity
            "--label_states", "False", # Don't add labels
            "--enumerate", "False" # Take the most likely state at target pH
        ]
        try:
            _run_subprocess(dimorphite_cmd)
            current_sdf_path = temp_protonated_sdf # Use protonated file for next steps
            print(f"[info] Protonation complete. Using {current_sdf_path}")
        except Exception as e:
            print(f"[warn] dimorphite-dl failed: {e}. Proceeding with original input SDF: {sdf_path}")
            current_sdf_path = sdf_path # Fallback to original if protonation fails
    else:
        print("[warn] dimorphite-dl not found in PATH. Skipping pH-specific protonation.")
        current_sdf_path = sdf_path

    # --- Step 2: Load molecule, generate 3D if needed, calculate charges ---
    print(f"[info] Reading molecule(s) from {current_sdf_path}...")
    # Use ForwardSDMolSupplier to handle potential errors gracefully
    suppl = Chem.ForwardSDMolSupplier(str(current_sdf_path), removeHs=False, sanitize=True) # Keep Hs, try sanitizing

    processed_mols = []
    has_3d = False # Track if at least one molecule has 3D coords

    for mol in suppl:
        if mol is None:
            print("[warn] Failed to read a molecule from SDF, skipping.")
            continue

        # Check for 3D coordinates
        try:
            conf = mol.GetConformer()
            if conf.Is3D():
                has_3d = True
            else:
                 print("[info] Molecule lacks 3D coordinates.")
        except ValueError: # No conformer exists
            print("[info] Molecule lacks 3D coordinates (no conformer).")
            has_3d = False # Ensure this is reset if any molecule lacks coords

        # Generate 3D if needed and requested
        if not has_3d and force_3d:
            print("[info] Attempting to generate 3D coordinates...")
            try:
                mol_h = Chem.AddHs(mol, addCoords=True) # Add Hs explicitly before embedding
                AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
                AllChem.MMFFOptimizeMolecule(mol_h) # Basic optimization
                mol = mol_h # Use the molecule with generated 3D coords
                print("[info] 3D coordinates generated.")
                has_3d = True # Mark as having 3D now
            except Exception as e:
                print(f"[warn] Failed to generate 3D coordinates: {e}. Proceeding without 3D.")


        # Calculate partial charges
        print(f"[info] Calculating '{charge_method}' partial charges...")
        try:
            if charge_method.lower() == 'gasteiger':
                AllChem.ComputeGasteigerCharges(mol)
                print("[info] Gasteiger charges computed.")
                # Charges are stored as properties on atoms, SDF writer should pick them up.
            else:
                print(f"[warn] Charge method '{charge_method}' not implemented in this script. Skipping charge calculation.")
        except Exception as e:
             print(f"[warn] Failed to compute charges: {e}. Proceeding without charges.")

        processed_mols.append(mol)

    if not processed_mols:
        raise RuntimeError(f"No valid molecules could be processed from the input SDF: {current_sdf_path}")

    # --- Step 3: Write processed molecules to output SDF ---
    print(f"[info] Writing {len(processed_mols)} processed molecule(s) to {outp_sdf}...")
    with Chem.SDWriter(str(outp_sdf)) as writer:
        writer.SetForceV3000(True) # Use V3000 format for better compatibility
        for mol in processed_mols:
            # Set molecule name if possible
            if mol.HasProp("_Name"):
                pass # Already has name
            elif sdf_path.stem:
                 mol.SetProp("_Name", sdf_path.stem)

            # RDKit's SDWriter automatically includes computed charges and coordinates.
            writer.write(mol)

    # --- Cleanup temp file ---
    if temp_protonated_sdf.exists() and temp_protonated_sdf != current_sdf_path:
        try:
            temp_protonated_sdf.unlink()
        except OSError:
            print(f"[warn] Could not remove temporary file: {temp_protonated_sdf}")

    print(f"[info] Ligand preparation complete. Output: {outp_sdf}")
    return outp_sdf
