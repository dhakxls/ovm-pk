"""Script with proper Pipeline initialization."""
import logging
from pathlib import Path
from ovmpk import Pipeline

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    filename='runlogs/pipeline.log'
)

# Hardcoded config matching requirements
CONFIG = {
    "physics": {
        "metal_model": {
            "mode": "MCPB",
            "oxidation_state": "Fe3+",
            "coordination": 6,
            "qmmm_method": "B3LYP/def2-TZVP"
        }
    },
    "docking": {
        "box_size_A": 18.0,
        "exhaustiveness": 32
    }
}

if __name__ == "__main__":
    try:
        logging.info("Starting pipeline")
        Pipeline(CONFIG, work_dir="data/work").run()
    except Exception as e:
        logging.error(f"Pipeline failed: {str(e)}")
        raise
