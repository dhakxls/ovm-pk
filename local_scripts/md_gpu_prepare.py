#!/usr/bin/env python3
import argparse, sys, os, json
from pathlib import Path

from simtk import unit
from openmm import app, Platform, XmlSerializer
import openmm as mm

# Toolkit
from rdkit import Chem
from openff.toolkit.topology import Molecule
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

try:
    import pdbfixer
    HAS_FIXER = True
except Exception:
    HAS_FIXER = False

def pick_platform():
    for name in ("CUDA", "OpenCL", "CPU"):
        if name in [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]:
            return name
    return "CPU"

def load_ligand_from_sdf(sdf_path: Path, pose_index: int):
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = None
    # pose_index is 1-based (as we’ve been reporting)
    idx0 = max(0, pose_index - 1)
    if idx0 < len(suppl):
        mol = suppl[idx0]
    if mol is None:
        raise RuntimeError(f"Could not read pose {pose_index} from {sdf_path}")
    # RDKit -> OpenFF Molecule
    offmol = Molecule.from_rdkit(mol, allow_undefined_stereo=True)
    return mol, offmol

def write_rdkit_pdb(mol, out_pdb: Path):
    w = Chem.rdmolfiles.PDBWriter(str(out_pdb))
    try:
        w.write(mol)
    finally:
        w.close()

def try_build_forcefield(off_lig: Molecule):
    """
    Amber14SB for protein, TIP3P water, SMIRNOFF for ligand.
    Try to add heme parameters if available.
    """
    smirnoff = SMIRNOFFTemplateGenerator(molecules=[off_lig], forcefield='openff-2.0.0.offxml')
    heme_ffxml_candidates = [
        # common locations shipped by openmmforcefields
        "amber/hemhemea.xml",
        "amber/heme.xml",
        "openmmforcefields/ffxml/amber/hemhemea.xml",
        "openmmforcefields/ffxml/amber/heme.xml",
        "openmmforcefields/ffxml/heme.xml",
    ]
    base = ["amber14/protein.ff14SB.xml", "amber14/tip3pfb.xml"]
    ff = None
    heme_loaded = None
    for h in [None] + heme_ffxml_candidates:
        try:
            files = base if h is None else base + [h]
            ff = app.ForceField(*(lambda _ff=[*files]: (print(f"[info] base FF: {_ff}"), _ff.extend(_find_heme_ffxmls()) or print(f"[info] +heme FFs: {_find_heme_ffxmls()}"), tuple(_ff))[-1])())
            ff.registerTemplateGenerator(smirnoff.generator)
            heme_loaded = h
            break
        except Exception:
            continue
    if heme_loaded:
        print(f"[info] Heme parameters loaded from: {heme_loaded}")
    else:
        print("[warn] Could not load a heme ffxml. If your PDB has HEM/HEC, "
              "system creation may fail. Install openmmforcefields fully or provide a heme ffxml.")
    return ff

def main():
    ap = argparse.ArgumentParser(description="Prepare solvated MD system (GPU-ready) from receptor + docked SDF pose")
    ap.add_argument("--receptor", required=True, help="Receptor PDB (keep HEM in place, e.g. results/.../receptor.pdb)")
    ap.add_argument("--sdf", required=True, help="Docked poses SDF")
    ap.add_argument("--pose-index", type=int, required=True, help="1-based index of best pose in SDF")
    ap.add_argument("--out", required=True, help="Output directory, e.g. results/md/ketoconazole_pH74")
    ap.add_argument("--ph", type=float, default=7.4, help="Protonation pH for adding hydrogens (PDBFixer)")
    ap.add_argument("--padding", type=float, default=10.0, help="Solvent padding (Å)")
    ap.add_argument("--ionic", type=float, default=0.15, help="Ionic strength (M)")
    args = ap.parse_args()

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    rec_pdb = Path(args.receptor)
    sdf = Path(args.sdf)
    if not rec_pdb.exists():
        sys.exit(f"[FATAL] receptor not found: {rec_pdb}")
    if not sdf.exists():
        sys.exit(f"[FATAL] sdf not found: {sdf}")

    # --- ligand ---
    rdkit_lig, off_lig = load_ligand_from_sdf(sdf, args.pose_index)
    lig_pdb = outdir / "ligand_pose.pdb"
    write_rdkit_pdb(rdkit_lig, lig_pdb)

    # --- protein (with heme) ---
    print("[info] loading receptor PDB")
    prot = app.PDBFile(str(rec_pdb))

    # Optional: fix, add hydrogens (keeps HEM if present)
    if HAS_FIXER:
        print("[info] running PDBFixer (addH, missing atoms)")
        fixer = pdbfixer.PDBFixer(filename=str(rec_pdb))
        # DO NOT remove heterogens so we keep HEM
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=args.ph)
        with open(outdir / "receptor_fixed.pdb", "w") as fh:
            app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)
        prot = app.PDBFile(str(outdir / "receptor_fixed.pdb"))

    # --- merge receptor + ligand into a single Modeller ---
    print("[info] merging receptor + ligand")
    lig = app.PDBFile(str(lig_pdb))
    modeller = app.Modeller(prot.topology, prot.positions)
    modeller.add(lig.topology, lig.positions)

    # --- parameterize (Amber14SB + TIP3P + SMIRNOFF ligand + heme if available) ---
    ff = try_build_forcefield(off_lig)

    # --- solvate ---
    print("[info] adding solvent and ions")
    modeller.addSolvent(
        ff,
        model="tip3p",
        padding=args.padding * unit.angstroms,
        ionicStrength=args.ionic * unit.molar,
    )

    # --- system ---
    print("[info] creating System")
    system = ff.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0 * unit.nanometer,
        constraints=app.HBonds,
        hydrogenMass=3.0 * unit.amu,  # HMR allows 4 fs later if you want
    )
    system.addForce(mm.MonteCarloBarostat(1.0 * unit.bar, 310 * unit.kelvin))

    # --- quick minimization ---
    print("[info] minimizing")
    integrator = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
    platform = Platform.getPlatformByName(pick_platform())
    sim = app.Simulation(modeller.topology, system, integrator, platform)
    sim.context.setPositions(modeller.positions)
    sim.minimizeEnergy(maxIterations=2000)

    # --- serialize everything we need for production ---
    with open(outdir / "system.xml", "w") as f: f.write(XmlSerializer.serialize(system))
    with open(outdir / "integrator_min.xml", "w") as f: f.write(XmlSerializer.serialize(integrator))
    state = sim.context.getState(getPositions=True, getVelocities=True)
    with open(outdir / "state_min.xml", "w") as f: f.write(XmlSerializer.serialize(state))

    with open(outdir / "complex_solvated.pdb", "w") as fh:
        app.PDBFile.writeFile(modeller.topology, sim.context.getState(getPositions=True).getPositions(), fh)

    meta = {
        "receptor": str(rec_pdb),
        "sdf": str(sdf),
        "pose_index": args.pose_index,
        "ph": args.ph,
        "padding_A": args.padding,
        "ionic_M": args.ionic,
        "platform": platform.getName(),
    }
    (outdir / "prepare.json").write_text(json.dumps(meta, indent=2))
    print(f"[done] Prepared MD folder: {outdir}")

if __name__ == "__main__":
    main()
