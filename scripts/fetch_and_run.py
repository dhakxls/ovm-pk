"""Fetch protein and run pipeline."""
import os
import sys
import shutil
import requests
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from scripts.run_with_logging import CONFIG, Pipeline

def fetch_pdb(pdb_id: str):
    """Download PDB file directly from RCSB."""
    os.makedirs('data/input/proteins', exist_ok=True)
    input_path = Path(f"data/input/proteins/{pdb_id}.pdb")
    
    if not input_path.exists():
        print(f"Downloading {pdb_id}.pdb from RCSB...")
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url)
        response.raise_for_status()
        input_path.write_text(response.text)
        print(f"Saved to {input_path}")
    else:
        print(f"{pdb_id}.pdb already exists")
    
    # Copy to work directory
    work_path = Path(f"data/work/protein/{pdb_id}.pdb")
    work_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(input_path, work_path)
    
    # Verify copy
    if not work_path.exists():
        raise FileNotFoundError(f"Failed to copy to {work_path}")
    print(f"Verified copy at {work_path}")

def run():
    """Run the pipeline."""
    # Clean previous runs
    os.system('rm -rf data/work/ data/output/latest/ runlogs/pipeline.log')
    os.makedirs('data/work/protein', exist_ok=True)
    
    # Execute pipeline
    print("Starting pipeline...")
    Pipeline(CONFIG, work_dir="data/work").run()

if __name__ == "__main__":
    fetch_pdb("5VCC")  # From default.yaml
    run()
