"""Affinity calculations with temperature and multi-state support."""
import math
import numpy as np
from typing import Optional, List, Dict

R_KCAL_PER_MOLK = 1.98720425864083e-3  # kcal/(mol·K)

def dg_to_ki(dg_kcal: float, temperature_k: float = 310.0) -> float:
    """Convert ΔG (kcal/mol) to Ki (M) at given temperature."""
    return math.exp(dg_kcal / (R_KCAL_PER_MOLK * temperature_k))

def ki_to_dg(ki_m: float, temperature_k: float = 310.0) -> float:
    """Convert Ki (M) to ΔG (kcal/mol) at given temperature."""
    return R_KCAL_PER_MOLK * temperature_k * math.log(ki_m)

def ki_nm_to_dg(ki_nm: float, temperature_k: float = 310.0) -> float:
    """Convert Ki (nM) to ΔG (kcal/mol)."""
    return ki_to_dg(ki_nm * 1e-9, temperature_k)

def boltzmann_weighted_dg(microstate_results: List[Dict]) -> float:
    """Calculate pH-weighted ΔG from microstate simulations.
    
    Args:
        microstate_results: List of {'deltaG': float, 'fraction': float}
    """
    dgs = np.array([x['deltaG'] for x in microstate_results])
    weights = np.array([x['fraction'] for x in microstate_results])
    
    # Boltzmann-weighted average
    return -np.log(np.sum(weights * np.exp(-dgs)))

def ki_with_microstates(microstate_results: List[Dict], temp_k: float) -> float:
    """Calculate net Ki considering all microstates."""
    net_dg = boltzmann_weighted_dg(microstate_results)
    return dg_to_ki(net_dg, temp_k)
