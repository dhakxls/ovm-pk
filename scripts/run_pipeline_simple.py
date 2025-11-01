"""Run pipeline with complete configuration using run directories."""

import sys
from datetime import datetime
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from ovmpk import Pipeline
from ovmpk.utils.run_dirs import set_run_context

# Complete configuration
CONFIG = {
    "physics": {
        "metal_model": {
            "mode": "MCPB",
            "oxidation_state": "Fe3+",
            "coordination": 6,
            "qmmm_method": "B3LYP/def2-TZVP",
        }
    },
    "system": {
        "pdb_id": "5VCC",
        "chain_id": "A",
    },
    "reporting": {
        "outdir": "results/reports",
        "run_root": "results/runs",
    },
}


def _ensure_run_context() -> None:
    run_root = Path("test_run").resolve()
    run_id = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    set_run_context(run_root, run_id)


if __name__ == "__main__":
    _ensure_run_context()
    Pipeline(CONFIG).run()
