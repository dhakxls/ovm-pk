#!/usr/bin/env python3
"""
Simple script to test protein and ligand fetcher functions directly.
Loads config, calls fetchers, checks basic output validity.
"""

import sys
import os
from pathlib import Path
import yaml

# --- Diagnostics ---
print("--- Environment Diagnostics ---")
print(f"Python Executable: {sys.executable}")
print("sys.path:")
for p in sys.path:
    print(f"  {p}")
# Try importing rcsbapi directly here
try:
    import rcsbapi
    from rcsbapi.search import TextQuery # Try a specific import too
    print("Direct import of rcsbapi in test_fetchers.py: SUCCESS")
    HAS_RCSB_API_TEST = True
except ImportError as e:
    print(f"Direct import of rcsbapi in test_fetchers.py: FAILED - {e}")
    HAS_RCSB_API_TEST = False
print("--- End Diagnostics ---\n")
# --- End Diagnostics ---


# Assuming the script runs from the main project directory (ovm-pk)
# Adjust if running from a different location
try:
    from src.ovmpk.fetchers import protein_fetcher, ligand_fetcher
except ImportError as e:
    print(f"Error importing fetcher modules: {e}")
    print("Make sure you are running this script from the main 'ovm-pk' directory.")
    sys.exit(1)

def main():
    print("--- Testing Fetchers ---")
    config_path = Path("configs/prod_test.yaml")
    print(f"Using config: {config_path}")
    if not config_path.exists():
        print(f"Error: Config file not found at {config_path}")
        return 1

    try:
        cfg = yaml.safe_load(config_path.read_text())
        print("Config loaded successfully.")
    except Exception as e:
        print(f"Error loading config file {config_path}: {e}")
        return 1

    # --- Test Protein Fetcher ---
    protein_success = False
    protein_identifier = "CYP3A4" # Example identifier from config/usage
    print(f"\nTesting protein fetcher for identifier: {protein_identifier}")
    try:
        # Create directories manually first to isolate fetcher issues
        Path("data/input/proteins").mkdir(parents=True, exist_ok=True)

        protein_result_dict = protein_fetcher.fetch(protein_identifier, cfg)
        print(f"Protein fetcher result: {protein_result_dict}")

        # Check the result - expects a dict like {'best_match': Path(...)}
        result_key = next(iter(protein_result_dict)) # Get the key ('best_match' or 'dry_run_placeholder')
        protein_file_path = protein_result_dict.get(result_key)

        if protein_file_path and isinstance(protein_file_path, Path) and protein_file_path.exists() and protein_file_path.stat().st_size > 0:
            print(f"Success: Protein file exists and is not empty: {protein_file_path}")
            protein_success = True
            # Optional: Add check for placeholder content if needed
            # if protein_fetcher._is_placeholder(protein_file_path):
            #     print("[info] Fetched file is a placeholder (as expected in dry run).")
            # else:
            #     print("[info] Fetched file is likely real data.")
        elif protein_file_path:
            print(f"Failure: Protein file path returned ({protein_file_path}), but file is missing or empty.")
        else:
            print("Failure: Protein fetcher did not return a valid file path.")

    except ImportError as e:
        # Catch the specific error mentioned in the log output
        print(f"Error during protein fetching: {e}")
        if "rcsb-api" in str(e):
             print(" >> This suggests the 'rcsb-api' package is still not found by the fetcher module.")
             if HAS_RCSB_API_TEST:
                  print(" >> However, direct import in test_fetchers.py SUCCEEDED. Check for environment inconsistencies or module caching.")
             else:
                  print(" >> Direct import in test_fetchers.py also FAILED. Please ensure rcsb-api is installed correctly in the active environment.")
        else:
             print(f" >> An unexpected ImportError occurred: {e}")

    except Exception as e:
        print(f"Error during protein fetching: {e}")
        # Optional: Print traceback for more detail
        # import traceback
        # traceback.print_exc()

    # --- Test Ligand Fetcher ---
    ligand_success = False
    ligand_identifier = "ketoconazole" # Example identifier
    print(f"\nTesting ligand fetcher for identifier: {ligand_identifier}")
    try:
        # Create directories manually first
        Path("data/input/ligands").mkdir(parents=True, exist_ok=True)

        ligand_file_path = ligand_fetcher.fetch(ligand_identifier, cfg)
        print(f"Ligand fetcher result: {ligand_file_path}")

        if ligand_file_path and isinstance(ligand_file_path, Path) and ligand_file_path.exists() and ligand_file_path.stat().st_size > 0:
            print(f"Success: Ligand file exists and is not empty: {ligand_file_path}")
            ligand_success = True
            # Optional: Add check for placeholder content
            # if ligand_fetcher._is_ligand_placeholder(ligand_file_path):
            #      print("[info] Fetched file is a placeholder.")
            # else:
            #      print("[info] Fetched file is likely real data.")
        elif ligand_file_path:
            print(f"Failure: Ligand file path returned ({ligand_file_path}), but file is missing or empty.")
        else:
            print("Failure: Ligand fetcher did not return a valid file path.")

    except Exception as e:
        print(f"Error during ligand fetching: {e}")

    print("\n--- Fetcher Test Complete ---")
    return 0 if protein_success and ligand_success else 1

if __name__ == "__main__":
    # Ensure script runs relative to project root if possible
    project_root = Path(__file__).parent.resolve()
    if not (project_root / "src").is_dir():
         print(f"Warning: Running from {Path.cwd()}. Ensure this is the project root ('ovm-pk').")

    # Change working directory to project root if necessary? (Risky if paths aren't relative)
    # os.chdir(project_root)

    sys.exit(main())

