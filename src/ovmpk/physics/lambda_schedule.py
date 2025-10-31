"""Lambda schedule generation for ABFE."""
import numpy as np

def generate_lambda_schedule(config: dict) -> dict:
    """Create λ windows for decoupling."""
    return {
        "coulomb": np.linspace(0, 1, config["lambda_elec"]).tolist(),
        "vdw": np.linspace(0, 1, config["lambda_lj"]).tolist()
    }
