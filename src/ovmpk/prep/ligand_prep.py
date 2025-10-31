# ovmpk/src/ovmpk/prep/ligand_prep.py (Corrected dimorphite_dl call)

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
    1. Reads the first molecule from input SDF to get SMILES.
    2. Protonates the ligand SMILES at a specified pH using dimorphite-dl (if available).
    3. Reads the protonated SMILES output.
    4. Converts SMILES back to RDKit Mol object, adds Hydrogens.
    5. Calculates partial charges (Gasteiger default) using RDKit.
    6. Generates 3D coordinates if needed/requested.
    7. Writes the processed molecule to a new SDF file in the work directory.

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
    charge_method = prep_cfg.get("charge_method", "gasteiger")
    output_suffix = prep_cfg.get("output_suffix", f"_prepared_ph{target_ph}")
    force_3d = prep_cfg.get("force_3d", True)

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    outp_sdf = WORK_DIR / f"{sdf_path.stem}{output_suffix}.sdf"

    print(f"[info] Starting ligand preparation for {sdf_path}...")

    # --- Read Input SMILES ---
    initial_mol = None
    initial_smiles = None
    mol_name = sdf_path.stem # Default name
    try:
        # Read the first molecule to get its SMILES string
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
        if suppl and len(suppl) > 0:
            initial_mol = suppl[0]
            if initial_mol:
                if initial_mol.HasProp("_Name"):
                     mol_name = initial_mol.GetProp("_Name")
                initial_smiles = Chem.MolToSmiles(initial_mol)
                print(f"[info] Extracted initial SMILES: {initial_smiles}")
            else:
                raise ValueError("First molecule in SDF could not be read.")
        else:
            raise ValueError("SDF file is empty or could not be read.")
    except Exception as e:
        raise RuntimeError(f"Error reading initial SMILES from {sdf_path}: {e}")

    # --- Step 1: Protonation using dimorphite-dl ---
    dimorphite_bin = _need("dimorphite_dl")
    # Output file for dimorphite SMILES result
    temp_protonated_smiles_file = WORK_DIR / f"{sdf_path.stem}_temp_protonated.smi"
    current_smiles = initial_smiles # Start with the original SMILES

    if dimorphite_bin and initial_smiles:
        print(f"[info] Running dimorphite-dl for pH {target_ph}...")
        # --- CORRECTED dimorphite_cmd using SMILES input ---
        dimorphite_cmd = [
            dimorphite_bin,
            "--ph_min", str(target_ph),
            "--ph_max", str(target_ph),
            "--output_file", str(temp_protonated_smiles_file),
            "--precision", "1", # Match argument from previous attempt
            # The SMILES string should be the last positional argument
            initial_smiles
        ]
        # Arguments like --label_states and --enumerate seem unsupported based on error
        # --- END CORRECTION ---
        try:
            _run_subprocess(dimorphite_cmd)
            # Read the protonated SMILES back from the output file
            if temp_protonated_smiles_file.exists():
                with open(temp_protonated_smiles_file, 'r') as f:
                    # Read first line, strip whitespace, get SMILES (usually first part)
                    line = f.readline().strip()
                    if line:
                        protonated_smiles = line.split()[0] # Assume SMILES is the first part
                        if protonated_smiles:
                            current_smiles = protonated_smiles
                            print(f"[info] Protonation complete. Using protonated SMILES: {current_smiles}")
                        else:
                            print("[warn] Dimorphite output file was empty. Using original SMILES.")
                    else:
                        print("[warn] Dimorphite output file was empty. Using original SMILES.")
            else:
                 print("[warn] Dimorphite output file not found. Using original SMILES.")

        except Exception as e:
            print(f"[warn] dimorphite-dl failed: {e}. Proceeding with original input SMILES: {initial_smiles}")
            current_smiles = initial_smiles # Fallback to original SMILES
    else:
        if not initial_smiles:
             print("[warn] Could not extract initial SMILES. Skipping protonation.")
        elif not dimorphite_bin:
             print("[warn] dimorphite-dl not found in PATH. Skipping pH-specific protonation.")
        current_smiles = initial_smiles

    # --- Step 2: Reconstruct Mol, Add H, Generate 3D if needed, Calculate charges ---
    print(f"[info] Processing SMILES: {current_smiles}")
    mol = Chem.MolFromSmiles(current_smiles)
    if mol is None:
        raise RuntimeError(f"Could not generate RDKit molecule from SMILES: {current_smiles}")

    # Add Hydrogens (important after SMILES conversion and before 3D/charges)
    mol = Chem.AddHs(mol, addCoords=force_3d) # Add coords only if forcing 3D
    print("[info] Added hydrogens to molecule.")

    # Check for 3D coordinates
    has_3d = False
    try:
        conf = mol.GetConformer()
        # Check if coordinates are meaningful (not all zero) - basic check
        if conf.Is3D() and any(abs(p.x) > 1e-6 or abs(p.y) > 1e-6 or abs(p.z) > 1e-6 for p in conf.GetPositions()):
             has_3d = True
        else:
             print("[info] Molecule lacks 3D coordinates after AddHs.")
    except ValueError: # No conformer exists
        print("[info] Molecule lacks 3D coordinates (no conformer).")

    # Generate 3D if needed and requested
    if not has_3d and force_3d:
        print("[info] Attempting to generate 3D coordinates...")
        try:
            # AddHs should have already happened, but ensure explicit H for embedding
            # mol_h = Chem.AddHs(mol, addCoords=True) # Already done
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
            AllChem.MMFFOptimizeMolecule(mol) # Basic optimization
            print("[info] 3D coordinates generated.")
            has_3d = True
        except Exception as e:
            print(f"[warn] Failed to generate 3D coordinates: {e}. Proceeding without guaranteed 3D.")
    elif not has_3d and not force_3d:
        print("[warn] force_3d is False and molecule lacks 3D coordinates. Output may not be suitable for docking.")


    # Calculate partial charges
    print(f"[info] Calculating '{charge_method}' partial charges...")
    try:
        if charge_method.lower() == 'gasteiger':
            # ComputeGasteigerCharges requires sanitized mol with Hydrogens
            # Sanitization happens implicitly in MolFromSmiles/AddHs if possible
            AllChem.ComputeGasteigerCharges(mol)
            print("[info] Gasteiger charges computed.")
        else:
            print(f"[warn] Charge method '{charge_method}' not implemented. Skipping charge calculation.")
    except Exception as e:
         print(f"[warn] Failed to compute charges: {e}. Proceeding without charges.")

    # Set molecule name
    mol.SetProp("_Name", mol_name)

    # --- Step 3: Write processed molecule to output SDF ---
    print(f"[info] Writing processed molecule to {outp_sdf}...")
    with Chem.SDWriter(str(outp_sdf)) as writer:
        writer.SetForceV3000(True)
        writer.write(mol)

    # --- Cleanup temp file ---
    if temp_protonated_smiles_file.exists():
        try:
            temp_protonated_smiles_file.unlink()
        except OSError:
            print(f"[warn] Could not remove temporary file: {temp_protonated_smiles_file}")

    print(f"[info] Ligand preparation complete. Output: {outp_sdf}")
    return outp_sdf


"""Ligand prep with microstate support."""
from pathlib import Path
from typing import Dict
from .microstates import generate_microstates, dominant_microstate

def prepare_ligand(sdf_file: Path, config: Dict) -> Dict:
    """Process ligand with pH awareness."""
    mol = Chem.MolFromMolFile(str(sdf_file))
    
    # Generate microstates
    microstates = generate_microstates(
        mol, 
        config['conditions']['pH']
    )
    
    # Select dominant state for docking
    dominant = dominant_microstate(microstates)
    
    return {
        'microstates': microstates,
        'dominant': dominant,
        'original_file': sdf_file
    }
