"""Adaptive sampling for ABFE."""
import numpy as np
from typing import Dict, List
from pymbar import MBAR

class AdaptiveSampler:
    """Dynamically allocate sampling based on uncertainty."""
    
    def __init__(self, config: Dict):
        self.config = config
        self.min_samples = 100  # Minimum steps per window
        
    def analyze(self, trajectories: List, energies: List) -> Dict[int, float]:
        """Identify under-sampled λ windows."""
        mbar = MBAR(energies)
        uncertainties = mbar.computeUncertainties()
        
        # Return dict of {window_idx: additional_time_ns}
        return {
            i: max(0.5, 2 * unc)
            for i, unc in enumerate(uncertainties)
            if unc > 0.2  # Threshold
        }
