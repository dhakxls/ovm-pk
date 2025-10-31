"""Simulation analysis utilities."""
from pymbar import MBAR
import numpy as np
from typing import List, Dict

def calculate_overlaps(energies: List[np.ndarray]) -> Dict[str, np.ndarray]:
    """Compute state overlap matrix.
    
    Returns:
        Dict with 'matrix' and 'neff' (effective samples)
    """
    mbar = MBAR(energies)
    matrix = mbar.computeOverlap()['matrix']
    
    return {
        'matrix': matrix,
        'neff': np.diag(mbar.computeEffectiveSamplesNumber())
    }

def check_convergence(overlaps: Dict, threshold: float = 0.3) -> bool:
    """Determine if sampling is converged."""
    min_overlap = np.min(overlaps['matrix'][np.nonzero(overlaps['matrix'])])
    return min_overlap > threshold
