"""
Core calculation utilities for molecular workflows
"""
from openmm import unit

def require_openmm():
    """Validate OpenMM installation"""
    try:
        import openmm
        return True
    except ImportError:
        raise RuntimeError("OpenMM is required but not installed")

def calculate_binding_affinity(deltaG: float) -> float:
    """Convert ΔG to predicted K_i (nM)"""
    R = 8.314 * unit.joule/(unit.mole*unit.kelvin)  # Gas constant
    T = 300 * unit.kelvin  # Standard temperature
    return (1e9 * unit.mole/unit.liter * unit.exp(deltaG/(R*T))).value_in_unit(unit.nanomolar)
