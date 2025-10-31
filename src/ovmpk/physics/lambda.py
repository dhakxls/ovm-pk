"""Lambda schedule generation for ABFE."""
from typing import List, Dict
import numpy as np

def generate_lambda_schedule(config: Dict) -> Dict[str, List[float]]:
    """Create λ windows for decoupling.
    
    Args:
        config: Must contain 'lambda_elec' and 'lambda_lj' window counts
    
    Returns:
        Dictionary with 'coulomb', 'vdw', 'restraint' schedules
    """
    return {
        "coulomb": np.linspace(0, 1, config["lambda_elec"]).tolist(),
        "vdw": np.linspace(0, 1, config["lambda_lj"]).tolist(),
        "restraint": [1.0, 0.8, 0.6, 0.4, 0.2, 0.0]  # Boresch release
    }
