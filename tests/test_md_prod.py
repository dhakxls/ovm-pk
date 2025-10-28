from __future__ import annotations
import pytest
pytestmark = pytest.mark.order(6)

# tests/test_md_prod.py
import os
from pathlib import Path
import math
import json

import openmm as mm
from openmm import app, XmlSerializer, unit

REPO = Path(__file__).resolve().parents[1]
PREP_DIR = REPO / "data" / "output" / "md_prepare_step_test"
RUN_DIR  = REPO / "data" / "output" / "md_prod_test"

def _pick_platform(preferred: str | None = "CUDA") -> tuple[mm.Platform, dict]:
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    order = ([preferred] if preferred else []) + ["CUDA", "OpenCL", "CPU"]
    for name in order:
        if name in names:
            # Use mixed precision on CUDA if possible
            props = {"Precision": "mixed"} if name == "CUDA" else {}
            return mm.Platform.getPlatformByName(name), props
    return mm.Platform.getPlatformByName("CPU"), {}

def _positions_are_finite(positions):
    for v in positions:
        x, y, z = v.value_in_unit(unit.nanometer)
        if not (math.isfinite(x) and math.isfinite(y) and math.isfinite(z)):
            return False
    return True

def _quantity_is_finite(q, u=unit.kilojoule_per_mole):
    try:
        val = q.value_in_unit(u)
        return math.isfinite(val)
    except Exception:
        return False

def test_short_production():
    # Inputs from prep
    sysxml = PREP_DIR / "system.xml"
    pdbmin = PREP_DIR / "complex_solvated.pdb"
    if not pdbmin.exists():
        pdbmin = PREP_DIR / "complex_solvated_premin.pdb"
    prep_meta = PREP_DIR / "prepare.json"
    for p in [sysxml, pdbmin, prep_meta]:
        assert p.exists(), f"Missing required file: {p}"

    RUN_DIR.mkdir(parents=True, exist_ok=True)

    # Load system & positions
    system = XmlSerializer.deserialize(sysxml.read_text())
    pdb = app.PDBFile(str(pdbmin))
    topology = pdb.topology
    positions = pdb.positions

    # Attach barostat for NPT production
    temperature = 310 * unit.kelvin
    baro = mm.MonteCarloBarostat(1.0 * unit.atmosphere, temperature, 25)
    system.addForce(baro)

    friction = 1.0 / unit.picosecond
    dt = 0.002 * unit.picoseconds
    steps = int(os.environ.get("OVMPK_PROD_STEPS", "10000"))  # ~20 ps

    platform, props = _pick_platform("CUDA")
    integrator = mm.LangevinMiddleIntegrator(temperature, friction, dt)
    sim = app.Simulation(topology, system, integrator, platform, props)
    sim.context.setPositions(positions)
    sim.context.setVelocitiesToTemperature(temperature)

    # Reporters
    dcd = RUN_DIR / "prod.dcd"
    log = RUN_DIR / "prod.log"
    chk = RUN_DIR / "prod.chk"

    sim.reporters.append(app.DCDReporter(str(dcd), 250))
    sim.reporters.append(app.StateDataReporter(
        str(log), 250, step=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, density=True, speed=True, separator=","
    ))

    # Run & checkpoint
    sim.step(steps)
    sim.saveCheckpoint(str(chk))

    # Final sanity
    st = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True, enforcePeriodicBox=True)
    assert _positions_are_finite(st.getPositions()), "Production produced NaN positions"
    assert _quantity_is_finite(st.getPotentialEnergy()), "Potential energy not finite"
    assert _quantity_is_finite(st.getKineticEnergy()), "Kinetic energy not finite"

    # Files exist & non-empty
    assert dcd.exists() and dcd.stat().st_size > 0, "No trajectory written"
    assert log.exists() and log.stat().st_size > 0, "No log written"
    assert chk.exists() and chk.stat().st_size > 0, "No checkpoint written"

    # Temperature roughly around target (read from reporter log)
    with open(log) as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]
    header = lines[0]
    headers = [h.strip() for h in header.split(",")]
    temp_col = next((i for i, h in enumerate(headers) if "Temperature" in h), None)
    assert temp_col is not None, f"Temperature column not found in log header: {headers}"
    last_vals = [x.strip() for x in lines[-1].split(",")]
    T = float(last_vals[temp_col])
    assert 240.0 <= T <= 380.0, f"Temperature out of reasonable range: {T:.1f} K"