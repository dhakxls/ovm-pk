# tests/test_md_equil.py
import os
from pathlib import Path
import json
import math
import numpy as np

import openmm as mm
from openmm import app, XmlSerializer, unit

REPO = Path(__file__).resolve().parents[1]
PREP_DIR = REPO / "data" / "output" / "md_prepare_step_test"   # from the previous test
RUN_DIR  = REPO / "data" / "output" / "md_equil_test"

def _pick_platform(preferred: str | None = "CUDA") -> mm.Platform:
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    order = ([preferred] if preferred else []) + ["CUDA", "OpenCL", "CPU"]
    for name in order:
        if name in names:
            return mm.Platform.getPlatformByName(name)
    return mm.Platform.getPlatformByName("CPU")

def _positions_are_finite(positions):
    for v in positions:
        x, y, z = v.value_in_unit(unit.nanometer)
        if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
            return False
    return True

def _quantity_is_finite(q):
    try:
        val = q.value_in_unit(unit.kilojoule_per_mole)
        return math.isfinite(val)
    except Exception:
        return False

def test_short_equilibration():
    # --- Inputs from the previous stage ---
    assert PREP_DIR.exists(), f"Missing prep dir: {PREP_DIR}"
    sysxml = PREP_DIR / "system.xml"
    pdbmin = PREP_DIR / "complex_solvated.pdb"  # post-minimize
    if not pdbmin.exists():
        # fallback: pre-minimize positions (rarely needed)
        pdbmin = PREP_DIR / "complex_solvated_premin.pdb"
    prep_meta = PREP_DIR / "prepare.json"

    for p in [sysxml, pdbmin, prep_meta]:
        assert p.exists(), f"Missing required file: {p}"

    RUN_DIR.mkdir(parents=True, exist_ok=True)

    # --- Load prepared system & topology/positions ---
    system = XmlSerializer.deserialize(sysxml.read_text())
    pdb = app.PDBFile(str(pdbmin))
    topology = pdb.topology
    positions = pdb.positions

    # --- MD configuration (tiny, fast) ---
    temperature = 310 * unit.kelvin
    friction   = 1.0 / unit.picosecond
    dt         = 0.002 * unit.picoseconds

    # NVT warmup then NPT settle (very short)
    steps_nvt = int(os.environ.get("OVMPK_EQUIL_STEPS_NVT", "2000"))  # ~4 ps
    steps_npt = int(os.environ.get("OVMPK_EQUIL_STEPS_NPT", "3000"))  # ~6 ps

    # --- Platform (GPU if available) ---
    platform = _pick_platform("CUDA")

    # --- NVT: Langevin, set velocities, quick shakeout ---
    integrator = mm.LangevinMiddleIntegrator(temperature, friction, dt)
    sim = app.Simulation(topology, system, integrator, platform)
    sim.context.setPositions(positions)
    sim.context.setVelocitiesToTemperature(temperature)

    # Reporters (lightweight)
    dcd_path = RUN_DIR / "equil.dcd"
    log_path = RUN_DIR / "equil.log"
    sim.reporters.append(app.DCDReporter(str(dcd_path), 100))
    sim.reporters.append(app.StateDataReporter(
        str(log_path), 100, step=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, density=True, speed=True, separator=","
    ))

    # Run NVT
    sim.step(steps_nvt)

    # Sanity after NVT
    state = sim.context.getState(getPositions=True, getEnergy=True)
    assert _positions_are_finite(state.getPositions()), "NVT produced NaN positions"
    assert _quantity_is_finite(state.getPotentialEnergy()), "NVT potential energy not finite"
    assert _quantity_is_finite(state.getKineticEnergy()), "NVT kinetic energy not finite"

    # --- NPT: add barostat, continue briefly ---
    baro = mm.MonteCarloBarostat(1.0 * unit.atmosphere, temperature, 25)
    system.addForce(baro)
    # Recreate Simulation to ensure integrator/system updated cleanly
    integrator2 = mm.LangevinMiddleIntegrator(temperature, friction, dt)
    sim = app.Simulation(topology, system, integrator2, platform)
    sim.context.setPositions(state.getPositions())
    sim.context.setVelocitiesToTemperature(temperature)

    sim.reporters.append(app.DCDReporter(str(dcd_path), 100))  # continue same file
    sim.reporters.append(app.StateDataReporter(
        str(log_path), 100, step=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, density=True, speed=True, separator=",", append=True
    ))

    sim.step(steps_npt)

    # Final checks
    final = sim.context.getState(getPositions=True, getEnergy=True, getVelocities=True, enforcePeriodicBox=True)
    assert _positions_are_finite(final.getPositions()), "NPT produced NaN positions"
    assert _quantity_is_finite(final.getPotentialEnergy()), "NPT potential energy not finite"
    assert _quantity_is_finite(final.getKineticEnergy()), "NPT kinetic energy not finite"

    # Files were written and non-empty
    assert dcd_path.exists() and dcd_path.stat().st_size > 0, "No trajectory written"
    assert log_path.exists() and log_path.stat().st_size > 0, "No log written"

    # Metadata echo (optional): confirm platform selection recorded previously
    meta = json.loads(prep_meta.read_text())
    assert meta.get("platform") in ("CUDA", "OpenCL", "CPU")
