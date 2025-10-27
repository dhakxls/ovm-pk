# tests/test_docking.py
import argparse
import os
import logging
import yaml
from pathlib import Path
import sys
import shutil # Import shutil for copying results

# Assuming the main docking logic is in a function/class within smina_wrapper
from ovmpk.docking.smina_wrapper import run as run_smina_docking # Import 'run' and alias it
from ovmpk.utils.logging import get_logger # Use the logger setup

def find_latest_file_with_suffix(directory, suffix, extensions):
    """ Helper to find the most recently modified file matching suffix and extensions. """
    latest_file = None
    latest_mtime = 0
    search_path = Path(directory)

    if not search_path.is_dir():
         print(f"[warn] Directory not found for searching: {directory}")
         return None

    for ext in extensions:
        pattern = f"*{suffix}{ext}"
        matches = list(search_path.glob(pattern))
        for match in matches:
             try:
                 mtime = match.stat().st_mtime
                 if mtime > latest_mtime:
                     latest_mtime = mtime
                     latest_file = match
             except FileNotFoundError:
                 continue # Ignore if file disappears between glob and stat
             except Exception as e:
                 # Using print as logger might not be fully configured yet
                 print(f"[warn] Error stating file {match}: {e}")

    if latest_file:
        return str(latest_file)
    else:
        # Fallback: find latest file with just the extension if suffix match fails
        print(f"[warn] Could not find file with specific suffix '{suffix}'. Falling back to latest match with extensions {extensions}.")
        for ext in extensions:
             matches = list(search_path.glob(f"*{ext}"))
             # Sort by modification time, newest first
             matches.sort(key=lambda x: x.stat().st_mtime, reverse=True)
             if matches:
                  print(f"[warn] Using fallback file: {matches[0]}")
                  return str(matches[0]) # Return the latest modified file with the right extension
        return None # No matching file found

def main():
    """
    Main function to run smina docking test using a config file.
    """
    parser = argparse.ArgumentParser(description="Test smina docking using a config file.")
    parser.add_argument("--config", default="configs/prod_test.yaml", help="Path to the YAML configuration file.")
    parser.add_argument("--log_level", default="INFO", help="Logging level (e.g., DEBUG, INFO, WARNING)")

    args = parser.parse_args()

    # --- Load Configuration ---
    try:
        with open(args.config, 'r') as f:
            config = yaml.safe_load(f)
        print(f"[info] Using config file: {args.config}")
    except Exception as e:
        print(f"[error] Error loading config file {args.config}: {e}")
        return

    # --- Determine Output Directory ---
    base_output_dir = Path(config.get('output_base_dir', 'data/output/docking_test'))
    base_output_dir.mkdir(parents=True, exist_ok=True)

    # --- Configure Logging ---
    log_file = base_output_dir / "test_docking.log"
    log_level = args.log_level.upper()
    numeric_level = getattr(logging, log_level, None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')

    logger = get_logger("ovmpk_test_docking")
    logger.setLevel(numeric_level)

    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        ch = logging.StreamHandler(sys.stdout)
        ch.setLevel(numeric_level)
        formatter = logging.Formatter("[%(levelname)s] %(message)s")
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        print(f"[debug] Added console handler to logger.")

    try:
        if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(log_file) for h in logger.handlers):
            fh = logging.FileHandler(str(log_file), mode='w')
            fh.setLevel(numeric_level)
            file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            fh.setFormatter(file_formatter)
            logger.addHandler(fh)
            print(f"[debug] Added file handler ({log_file}) to logger.")
        else:
             print(f"[debug] File handler for {log_file} already exists.")
    except Exception as e:
        print(f"[error] Failed to configure file logging: {e}")
        logger.error("File logging configuration failed.")

    logger.info("Starting docking test...")
    logger.info(f"Using config: {args.config}")

    # --- Locate Input Files ---
    prep_config = config.get('prep', {})
    protein_prep_config = prep_config.get('protein', {})
    ligand_prep_config = prep_config.get('ligand', {})
    docking_config = config.get('docking', {})

    # Directory where protein_prep.py *saved* its output (needed by smina_wrapper internal prep)
    protein_prep_work_dir = Path("data/work/protein_prep")
    protein_suffix = protein_prep_config.get('output_suffix', '_fixed_ph7.4')
    # The smina wrapper's run() function expects the *initial* PDB/CIF and does its own PDBQT conversion.
    # So, find the PDB file generated by protein_prep.py (which includes PDBFixer fixes).
    # We rely on the smina_wrapper to handle the conversion of this PDB to PDBQT.
    prepared_protein_pdb = find_latest_file_with_suffix(
        protein_prep_work_dir,
        protein_suffix,
        ['.pdb']
    )
    if not prepared_protein_pdb or not os.path.exists(prepared_protein_pdb):
        logger.error(f"Could not find prepared protein PDB file in '{protein_prep_work_dir}' with suffix '{protein_suffix}'.")
        # Try finding the original CIF as a last resort? Maybe not needed if prep must succeed first.
        prepared_protein_pdb = None


    # Directory where ligand_prep.py saved its output
    ligand_prep_work_dir = Path("data/work/ligand_prep")
    ligand_suffix = ligand_prep_config.get('output_suffix', '_prepared_ph7.4')
    # The smina wrapper's run() function expects the prepared ligand SDF.
    prepared_ligand_sdf = find_latest_file_with_suffix(
        ligand_prep_work_dir,
        ligand_suffix,
        ['.sdf']
    )
    if not prepared_ligand_sdf or not os.path.exists(prepared_ligand_sdf):
        logger.error(f"Could not find prepared ligand SDF file in '{ligand_prep_work_dir}' with suffix '{ligand_suffix}'.")
        prepared_ligand_sdf = None


    # Define where final docking results should be collected
    docking_out_collect_dir = base_output_dir / "final_docking_poses"
    docking_out_collect_dir.mkdir(parents=True, exist_ok=True)

    # --- Run Docking ---
    # Pass the PDB from protein_prep and the SDF from ligand_prep to the run function.
    # The smina_wrapper's run() function will handle PDB->PDBQT and SDF->PDBQT conversions.
    if prepared_protein_pdb and prepared_ligand_sdf and docking_config:
        logger.info(f"Using Prepared Protein: {prepared_protein_pdb}")
        logger.info(f"Using Prepared Ligand: {prepared_ligand_sdf}")
        # Note: smina_wrapper will create outputs in its own WORK_DIR ('data/work/docking')

        try:
            list_of_output_sdfs = run_smina_docking(
                protein_pdb_in=Path(prepared_protein_pdb), # Pass the PDB processed by protein_prep
                ligand_sdf_in=Path(prepared_ligand_sdf),   # Pass the SDF processed by ligand_prep
                cfg=config
            )

            if list_of_output_sdfs and isinstance(list_of_output_sdfs, list):
                all_exist = all(os.path.exists(f) for f in list_of_output_sdfs)
                if all_exist:
                    logger.info(f"Docking finished successfully. Output files located in: {Path(list_of_output_sdfs[0]).parent}") # Log the directory
                    for f_path in list_of_output_sdfs:
                        logger.info(f"  - {Path(f_path).name}")
                        try:
                            # Copy results to the test output directory
                            shutil.copy(f_path, docking_out_collect_dir)
                            logger.info(f"    Copied to {docking_out_collect_dir}")
                        except Exception as copy_e:
                            logger.warning(f"    Failed to copy output file {f_path}: {copy_e}")
                else:
                    logger.error("Docking process completed but some output files are missing.")
                    logger.error(f"Expected files: {list_of_output_sdfs}")
            else:
                 logger.error("Docking process completed but did not return the expected list of output paths.")

        except Exception as e:
            logger.error(f"Error during docking: {e}", exc_info=True)
    else:
        missing = []
        if not prepared_protein_pdb: missing.append("prepared protein PDB file from protein_prep step")
        if not prepared_ligand_sdf: missing.append("prepared ligand SDF file from ligand_prep step")
        if not docking_config: missing.append("'docking' section in config")
        logger.warning(f"Skipping docking - missing prerequisites: {', '.join(missing)}.")

    logger.info("Docking test finished.")

if __name__ == "__main__":
    main()

    