from __future__ import annotations
import pytest
pytestmark = pytest.mark.order(12)

# tests/test_ki_to_dg.py
import math

R_KCAL = 0.0019872041  # kcal/mol/K

def ki_to_dg(k_i_molar: float, T: float = 298.15) -> float:
    """Standard ΔG° (kcal/mol) from Ki (M). ΔG° = RT ln(Kd/1M)."""
    return R_KCAL * T * math.log(k_i_molar)

def pki_to_dg(pKi: float, T: float = 298.15) -> float:
    """ΔG° from pKi where Ki = 10^-pKi (M)."""
    Ki = 10.0 ** (-pKi)
    return ki_to_dg(Ki, T)

def test_pki_roundtrip_room_temp():
    # pKi = 6 → Ki = 1e-6 M → ΔG° ~ -8.2 kcal/mol @ 298K
    dg = pki_to_dg(6.0, T=298.15)
    assert abs(dg - (-8.2)) < 0.2  # loose tolerance for constant precision

def test_ki_limits_monotonic():
    # Stronger binders (smaller Ki) should have more negative ΔG°
    dg_nM = ki_to_dg(1e-9)   # 1 nM
    dg_uM = ki_to_dg(1e-6)   # 1 µM
    dg_mM = ki_to_dg(1e-3)   # 1 mM
    assert dg_nM < dg_uM < dg_mM  # more negative < less negative

def test_temperature_effect_direction():
    # At higher T, magnitude of ΔG° grows for same Ki (since RT ln K)
    ki = 1e-6
    dg_cold = ki_to_dg(ki, T=278.15)
    dg_warm = ki_to_dg(ki, T=310.15)
    # both are negative; warm should be a bit more negative
    assert dg_warm < dg_cold