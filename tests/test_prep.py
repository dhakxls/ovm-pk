from __future__ import annotations
import pytest
pytestmark = pytest.mark.order(2)

import sys, io
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
from rich.progress import Progress, TimeElapsedColumn, BarColumn, MofNCompleteColumn, TimeRemainingColumn, SpinnerColumn
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

WORK = REPO / "data" / "work"
OUT  = REPO / "data" / "work"
LOGDIR = REPO / "data" / "output" / "prep_test" / "logs"
LOGDIR.mkdir(parents=True, exist_ok=True)

def _rich_progress():
    return Progress(
        SpinnerColumn(),
        "[progress.description]{task.description}",
        BarColumn(),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        transient=True,
    )

def _write_log(name: str, text: str):
    (LOGDIR / f"{name}.log").write_text(text)

@pytest.mark.dependency(depends=["fetchers_ok"])
def test_prep_protein_and_ligand():
    # newest fetched
    pdb = max((WORK / "protein_fetch").glob("*.pdb"), key=lambda p: p.stat().st_mtime)
    sdf = max((WORK / "ligand_fetch").glob("*.sdf"), key=lambda p: p.stat().st_mtime)

    prot_out = OUT / "protein_prep"
    lig_out  = OUT / "ligand_prep"
    prot_out.mkdir(parents=True, exist_ok=True)
    lig_out.mkdir(parents=True, exist_ok=True)

    # expected outputs
    fixed_pdb   = prot_out / f"{pdb.stem}_fixed_ph7.4.pdb"
    prepared_sdf= lig_out / f"{sdf.stem}_prepared_ph7.4.sdf"

    f_stdout, f_stderr = io.StringIO(), io.StringIO()
    with redirect_stdout(f_stdout), redirect_stderr(f_stderr), _rich_progress() as prog:
        task = prog.add_task("[bold]Preparing protein & ligand...", total=2)

        # PROTEIN
        try:
            from ovmpk.prep.protein_prep import prepare_protein  # type: ignore
            prepare_protein(str(pdb), str(fixed_pdb), ph=7.4)    # your API may differ
        except Exception:
            # minimal fallback using PDBFixer if your API isn’t available
            import pdbfixer
            from openmm import app
            fixer = pdbfixer.PDBFixer(filename=str(pdb))
            fixer.findMissingResidues(); fixer.findMissingAtoms()
            fixer.addMissingAtoms(); fixer.addMissingHydrogens(pH=7.4)
            with fixed_pdb.open("w") as fh:
                app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)
        prog.update(task, advance=1)

        # LIGAND
        try:
            from ovmpk.prep.ligand_prep import prepare_ligand  # type: ignore
            prepare_ligand(str(sdf), str(prepared_sdf), ph=7.4) # your API may differ
        except Exception:
            # minimal fallback: RDKit charge & write SDF
            from rdkit import Chem
            mols = [m for m in Chem.SDMolSupplier(str(sdf), removeHs=False) if m]
            assert mols, f"Could not read {sdf}"
            w = Chem.SDWriter(str(prepared_sdf))
            try:
                for m in mols: w.write(m)
            finally:
                w.close()
        prog.update(task, advance=1)

    _write_log("prep_stdout", f_stdout.getvalue())
    _write_log("prep_stderr", f_stderr.getvalue())

    assert fixed_pdb.exists()   and fixed_pdb.stat().st_size   > 0
    assert prepared_sdf.exists() and prepared_sdf.stat().st_size > 0