# tests/test_md_resume.py
import os
from pathlib import Path

import openmm as mm
from openmm import app, XmlSerializer, unit

REPO = Path(__file__).resolve().parents[1]
PREP_DIR = REPO / "data" / "output" / "md_prepare_step_test"
PROD_DIR = REPO / "data" / "output" / "md_prod_test"
RUN_DIR  = REPO / "data" / "output" / "md_resume_test"

def _pick_platform(preferred: str | None = "CUDA") -> tuple[mm.Platform, dict]:
    names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
    order = ([preferred] if preferred else []) + ["CUDA", "OpenCL", "CPU"]
    for name in order:
        if name in names:
            return mm.Platform.getPlatformByName(name), ({"Precision":"mixed"} if name=="CUDA" else {})
    return mm.Platform.getPlatformByName("CPU"), {}

def test_resume_from_checkpoint():
    sysxml = PREP_DIR / "system.xml"
    pdbref = PREP_DIR / "complex_solvated.pdb"
    chk    = PROD_DIR / "prod.chk"
    for p in [sysxml, pdbref, chk]:
        assert p.exists(), f"Missing required file: {p}"

    RUN_DIR.mkdir(parents=True, exist_ok=True)

    system = XmlSerializer.deserialize(sysxml.read_text())
    pdb = app.PDBFile(str(pdbref))
    topology = pdb.topology

    platform, props = _pick_platform("CUDA")

    integrator = mm.LangevinMiddleIntegrator(310*unit.kelvin, 1.0/unit.picosecond, 0.002*unit.picoseconds)
    sim = app.Simulation(topology, system, integrator, platform, props)
    sim.loadCheckpoint(str(chk))

    # Short continue
    steps2 = int(os.environ.get("OVMPK_RESUME_STEPS", "5000"))
    dcd2 = RUN_DIR / "resume.dcd"
    log2 = RUN_DIR / "resume.log"

    sim.reporters.append(app.DCDReporter(str(dcd2), 250))
    sim.reporters.append(app.StateDataReporter(
        str(log2), 250, step=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, density=True, speed=True, separator=","
    ))

    sim.step(steps2)

    assert dcd2.exists() and dcd2.stat().st_size > 0, "No resumed trajectory written"
    assert log2.exists() and log2.stat().st_size > 0, "No resumed log written"
