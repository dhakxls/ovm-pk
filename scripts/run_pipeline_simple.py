"""Run pipeline with complete configuration."""
import os
import sys
from pathlib import Path

# Add project root to path
sys.path.append(str(Path(__file__).parent.parent))

from ovmpk import Pipeline

# Complete configuration
CONFIG = {
    "physics": {
        "metal_model": {
            "mode": "MCPB",
            "oxidation_state": "Fe3+",
            "coordination": 6,
            "qmmm_method": "B3LYP/def2-TZVP"
        }
    },
    "system": {
        "pdb_id": "5VCC",
        "chain_id": "A"
    },
    "reporting": {
        "outdir": "results/reports",
        "run_root": "results/runs"
    }
}

if __name__ == "__main__":
    # Clean previous runs
    os.system('rm -rf data/work/ data/output/latest/ runlogs/pipeline.log')
    os.makedirs('data/work/protein', exist_ok=True)
    
    # Execute pipeline
    Pipeline(CONFIG, work_dir="data/work").run()
