"""Affinity conversion utilities."""
import math

def ki_to_dg(ki_nm: float, temp_k: float = 298.15) -> float:
    """Convert Ki in nM to ΔG in kcal/mol."""
    R = 0.0019872041  # kcal/(mol·K)
    ki_m = ki_nm * 1e-9
    return R * temp_k * math.log(ki_m)
