# tests/test_prep.py
import argparse
import os
import logging
import yaml
import glob
from pathlib import Path
import sys # Import sys

# --- CORRECTED IMPORTS ---
from ovmpk.prep.protein_prep import prepare as prepare_protein_function
from ovmpk.prep.ligand_prep import prepare as prepare_ligand_function
# --- CORRECTED LOGGING IMPORT ---
from ovmpk.utils.logging import get_logger
# ---

def find_file_by_identifier(directory, identifier_part):
    """ Helper to find a file in a directory based on part of its name. """
    matches = list(Path(directory).glob(f"*{identifier_part}*"))
    if matches:
        preferred_order = ['.cif', '.pdb', '.sdf']
        for ext in preferred_order:
            for match in matches:
                if match.suffix.lower() == ext:
                    return str(match)
        return str(matches[0])
    return None

def main():
    """
    Main function to run protein and ligand preparation tests using a config file.
    """
    parser = argparse.ArgumentParser(description="Test protein and ligand preparation using a config file.")
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
    base_output_dir = Path(config.get('output_base_dir', 'data/output/prep_test'))
    base_output_dir.mkdir(parents=True, exist_ok=True)

    # --- CORRECTED LOGGING SETUP ---
    log_file = base_output_dir / "test_prep.log"
    log_level = args.log_level.upper()
    numeric_level = getattr(logging, log_level, None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level}')

    # Get the logger using the provided function
    logger = get_logger("ovmpk_test_prep") # Use a specific name if desired
    logger.setLevel(numeric_level)

    # Remove existing handlers (optional, prevents duplicate console logs if run multiple times)
    # for handler in logger.handlers[:]:
    #     logger.removeHandler(handler)

    # Add console handler (the get_logger function might already do this)
    # Check if a StreamHandler already exists
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        ch = logging.StreamHandler(sys.stdout) # Explicitly use stdout
        ch.setLevel(numeric_level)
        formatter = logging.Formatter("[%(levelname)s] %(message)s")
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        print(f"[debug] Added console handler to logger.") # Simple print for immediate feedback

    # Add file handler
    try:
        # Check if a FileHandler for the same file already exists
        if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(log_file) for h in logger.handlers):
            fh = logging.FileHandler(str(log_file), mode='w') # Use 'w' to overwrite log each run
            fh.setLevel(numeric_level)
            file_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
            fh.setFormatter(file_formatter)
            logger.addHandler(fh)
            print(f"[debug] Added file handler ({log_file}) to logger.") # Simple print
        else:
             print(f"[debug] File handler for {log_file} already exists.")

    except Exception as e:
        print(f"[error] Failed to configure file logging: {e}")
        logger.error("File logging configuration failed.")

    # Now use logger instead of logging directly
    logger.info("Starting preparation test...")
    logger.info(f"Using config: {args.config}")
    # --- END LOGGING SETUP CORRECTION ---


    # --- Locate Input Files based on Config ---
    # Find protein file
    protein_fetcher_config = config.get('fetch', {}).get('protein', {})
    protein_prep_config = config.get('prep', {}).get('protein', {})
    protein_input_dir = protein_prep_config.get('input_dir', 'data/input/proteins')
    protein_id = protein_fetcher_config.get('pdb_id', None)
    if not protein_id:
         uniprot_id = "P08684" # Hardcoded based on 5VCC example
         protein_identifier_part = uniprot_id
         logger.warning(f"pdb_id not explicitly in fetch config, searching for file containing '{protein_identifier_part}' in {protein_input_dir}")
    else:
         protein_identifier_part = protein_id
    protein_file = find_file_by_identifier(protein_input_dir, protein_identifier_part)
    if not protein_file or not os.path.exists(protein_file):
        logger.error(f"Could not find input protein file matching '{protein_identifier_part}' in '{protein_input_dir}'. Searched pattern: *{protein_identifier_part}*")
        protein_file = next(iter(Path(protein_input_dir).glob('*.cif')), None) or \
                       next(iter(Path(protein_input_dir).glob('*.pdb')), None)
        if protein_file:
             protein_file = str(protein_file)
             logger.warning(f"Falling back to using first found protein file: {protein_file}")
        else:
             logger.error("No protein file found in input directory.")
             protein_file = None

    # Find ligand file
    ligand_fetcher_config = config.get('fetch', {}).get('ligand', {})
    ligand_prep_config = config.get('prep', {}).get('ligand', {})
    ligand_input_dir = ligand_prep_config.get('input_dir', 'data/input/ligands')
    ligand_identifier_part = ligand_fetcher_config.get('identifier', 'ligand').replace(" ", "_")
    ligand_file = find_file_by_identifier(ligand_input_dir, f"{ligand_identifier_part}_name")
    if not ligand_file or not os.path.exists(ligand_file):
        logger.error(f"Could not find input ligand file matching '{ligand_identifier_part}_name' in '{ligand_input_dir}'. Searched pattern: *{ligand_identifier_part}_name*")
        ligand_file = next(iter(Path(ligand_input_dir).glob('*.sdf')), None)
        if ligand_file:
             ligand_file = str(ligand_file)
             logger.warning(f"Falling back to using first found ligand file: {ligand_file}")
        else:
             logger.error("No ligand file found in input directory.")
             ligand_file = None

    # --- Protein Preparation Test ---
    if protein_file and 'protein' in config.get('prep', {}):
        logger.info(f"Preparing protein: {protein_file}")
        try:
            protein_out_dir = Path(config['prep']['protein'].get('output_dir', base_output_dir / "prepared_protein"))
            protein_out_dir.mkdir(parents=True, exist_ok=True)
            protein_input_paths = {'apo': Path(protein_file)}
            prepared_protein_file_path = prepare_protein_function(
                paths=protein_input_paths,
                cfg=config
            )
            prepared_protein_file = str(prepared_protein_file_path) if prepared_protein_file_path else None
            if prepared_protein_file and os.path.exists(prepared_protein_file):
                logger.info(f"Protein preparation successful. Output: {prepared_protein_file}")
            elif prepared_protein_file == str(Path(protein_file)):
                 logger.info(f"Protein preparation skipped or failed, using original file: {prepared_protein_file}")
            else:
                logger.error("Protein preparation failed or produced no output file.")
        except Exception as e:
            logger.error(f"Error during protein preparation: {e}", exc_info=True)
    else:
        logger.warning("Skipping protein preparation - input file not found or 'prep.protein' section missing in config.")

    # --- Ligand Preparation Test ---
    if ligand_file and 'ligand' in config.get('prep', {}):
        logger.info(f"Preparing ligand: {ligand_file}")
        try:
            ligand_out_dir = Path(config['prep']['ligand'].get('output_dir', base_output_dir / "prepared_ligand"))
            ligand_out_dir.mkdir(parents=True, exist_ok=True)
            prepared_ligand_file_path = prepare_ligand_function(
                sdf_path=Path(ligand_file),
                cfg=config
            )
            prepared_ligand_file = str(prepared_ligand_file_path) if prepared_ligand_file_path else None
            if prepared_ligand_file and os.path.exists(prepared_ligand_file):
                logger.info(f"Ligand preparation successful. Output: {prepared_ligand_file}")
            elif prepared_ligand_file == str(Path(ligand_file)):
                 logger.info(f"Ligand preparation skipped or failed, using original file: {prepared_ligand_file}")
            else:
                logger.error("Ligand preparation failed or produced no output file.")
        except Exception as e:
            logger.error(f"Error during ligand preparation: {e}", exc_info=True)
    else:
        logger.warning("Skipping ligand preparation - input file not found or 'prep.ligand' section missing in config.")

    logger.info("Preparation test finished.")

if __name__ == "__main__":
    main()
    