#!/usr/bin/env python3
import argparse, json
from pathlib import Path
from simtk import unit
from openmm import app, Platform, XmlSerializer
import openmm as mm

def pick_platform(prefer="CUDA"):
    names = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    for cand in (prefer, "OpenCL", "CPU"):
        if cand in names:
            return Platform.getPlatformByName(cand)
    return Platform.getPlatformByName("CPU")

def main():
    ap = argparse.ArgumentParser(description="Run GPU MD (equil + production) from prepared XMLs")
    ap.add_argument("--in", dest="indir", required=True, help="Folder with system.xml, state_min.xml")
    ap.add_argument("--steps", type=int, default=2_000_000, help="Production steps (2 fs -> 4 ns)")
    ap.add_argument("--t", type=float, default=310.0, help="Temperature (K)")
    ap.add_argument("--dt", type=float, default=0.002, help="Timestep (ps)")
    ap.add_argument("--friction", type=float, default=1.0, help="Friction (1/ps)")
    ap.add_argument("--eq-steps", type=int, default=100_000, help="Equilibration steps before production")
    ap.add_argument("--report", type=int, default=2_500, help="Report every N steps")
    ap.add_argument("--precision", default="mixed", choices=["single", "mixed", "double"], help="GPU precision")
    ap.add_argument("--dcd", default="traj.dcd", help="DCD filename")
    ap.add_argument("--log", default="md.log", help="StateData log")
    ap.add_argument("--chk", default="md.chk", help="Checkpoint file")
    args = ap.parse_args()

    indir = Path(args.indir)
    system = XmlSerializer.deserialize((indir / "system.xml").read_text())
    state  = XmlSerializer.deserialize((indir / "state_min.xml").read_text())
    pdb    = app.PDBFile(str(indir / "complex_solvated.pdb"))

    integrator = mm.LangevinMiddleIntegrator(
        args.t * unit.kelvin,
        args.friction / unit.picosecond,
        args.dt * unit.picoseconds
    )

    platform = pick_platform()
    props = {}
    if platform.getName() == "CUDA":
        props["Precision"] = args.precision
    print(f"[info] platform={platform.getName()} props={props}")

    sim = app.Simulation(pdb.topology, system, integrator, platform, props)
    sim.context.setPositions(state.getPositions())
    if state.getVelocities() is not None:
        sim.context.setVelocities(state.getVelocities())

    # Reporters
    sim.reporters.append(app.DCDReporter(str(indir / args.dcd), args.report))
    sim.reporters.append(app.StateDataReporter(
        str(indir / args.log), args.report,
        step=True, potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
        temperature=True, density=True, progress=True, speed=True, elapsedTime=True, separator="\t"
    ))
    sim.reporters.append(app.CheckpointReporter(str(indir / args.chk), args.report * 10))

    # Short equilibration
    print("[info] equilibration…")
    sim.step(args.eq_steps)

    # Production
    print("[info] production…")
    sim.step(args.steps)

    # Save final state
    final = sim.context.getState(getPositions=True, getVelocities=True)
    (indir / "state_final.xml").write_text(mm.XmlSerializer.serialize(final))
    print("[done] MD finished.")

if __name__ == "__main__":
    main()
