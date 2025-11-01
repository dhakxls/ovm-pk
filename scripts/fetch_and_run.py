"""Fetch protein and run pipeline using run-scoped directories."""

import sys
from datetime import datetime
from pathlib import Path

import requests

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from scripts.run_with_logging import CONFIG, Pipeline
from ovmpk.utils.run_dirs import set_run_context, stage_dir


def _ensure_run_context() -> None:
    run_root = Path("test_run").resolve()
    run_id = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    set_run_context(run_root, run_id)


def fetch_pdb(pdb_id: str) -> Path:
    """Download PDB file directly from RCSB into the run-scoped Stage 1 directory."""

    input_dir = Path("cache/input/proteins")
    input_dir.mkdir(parents=True, exist_ok=True)
    input_path = input_dir / f"{pdb_id}.pdb"

    if not input_path.exists():
        print(f"Downloading {pdb_id}.pdb from RCSB...")
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        input_path.write_text(response.text)
        print(f"Saved to {input_path}")
    else:
        print(f"{pdb_id}.pdb already cached at {input_path}")

    stage_path = stage_dir("stage1/protein") / f"{pdb_id}.pdb"
    stage_path.parent.mkdir(parents=True, exist_ok=True)
    stage_path.write_text(input_path.read_text())
    print(f"Staged PDB at {stage_path}")
    return stage_path


def run() -> None:
    """Run the pipeline in the current run directory."""

    print("Starting pipeline...")
    Pipeline(CONFIG).run()


if __name__ == "__main__":
    _ensure_run_context()
    fetch_pdb("5VCC")  # From default.yaml
    run()
