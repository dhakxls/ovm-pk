"""Scripted pipeline runner with logging and run-scoped directories."""

import logging
from datetime import datetime
from pathlib import Path

from ovmpk import Pipeline
from ovmpk.utils.run_dirs import set_run_context

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    filename="runlogs/pipeline.log",
)

# Hardcoded config matching requirements
CONFIG = {
    "physics": {
        "metal_model": {
            "mode": "MCPB",
            "oxidation_state": "Fe3+",
            "coordination": 6,
            "qmmm_method": "B3LYP/def2-TZVP",
        }
    },
    "docking": {
        "box_size_A": 18.0,
        "exhaustiveness": 32,
    },
}


def _ensure_run_context() -> None:
    run_root = Path("test_run").resolve()
    run_id = datetime.now().strftime("run_%Y%m%d_%H%M%S")
    set_run_context(run_root, run_id)


if __name__ == "__main__":
    _ensure_run_context()
    try:
        logging.info("Starting pipeline")
        Pipeline(CONFIG).run()
    except Exception as exc:  # pragma: no cover - CLI usage
        logging.error("Pipeline failed: %s", exc)
        raise
