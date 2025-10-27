import sys, io, json, math
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
from rich.progress import Progress, TimeElapsedColumn, BarColumn, MofNCompleteColumn, TimeRemainingColumn, SpinnerColumn
from rdkit import Chem

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

WORK = REPO / "data" / "work"
OUTDIR = REPO / "data" / "output" / "pose_selection_test"
LOGDIR = OUTDIR / "logs"
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

def _fe_n_distance(mol, fe_xyz):
    # finds the closest N to the given Fe coords (tuple)
    conf = mol.GetConformer()
    import numpy as np
    f = np.array(fe_xyz)
    best = math.inf
    for a in mol.GetAtoms():
        if a.GetSymbol() != "N": continue
        p = conf.GetAtomPosition(a.GetIdx())
        d = np.linalg.norm(f - np.array([p.x, p.y, p.z]))
        best = min(best, d)
    return best

def test_pose_selection_min_fe_n_distance():
    # Inputs: newest docking SDFs
    dock_dir = WORK / "docking"
    cand = sorted(dock_dir.glob("*.sdf"))
    assert cand, f"No docking outputs in {dock_dir}"
    # Fe center: infer from filename or use last used default
    # You can also read from configs/prod_test.yaml if you store it there.
    fe_xyz = (-15.846, -23.032, -11.293)

    sel_dir = OUTDIR / "selected_pose"
    sel_dir.mkdir(parents=True, exist_ok=True)

    stdout, stderr = io.StringIO(), io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr), _rich_progress() as prog:
        t = prog.add_task("[bold]Scanning poses...", total=len(cand))
        best = (None, None, None)  # (path, idx, dist)
        for sdf in cand:
            suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
            for i, m in enumerate(suppl, start=1):
                if m is None: continue
                try:
                    d = _fe_n_distance(m, fe_xyz)
                except Exception:
                    continue
                if best[2] is None or d < best[2]:
                    best = (sdf, i, d)
            prog.advance(t)

        assert best[0] is not None, "No valid pose found"
        src_sdf, pose_idx, dist = best

        # extract that pose to SDF
        sup2 = Chem.SDMolSupplier(str(src_sdf), removeHs=False)
        m = sup2[pose_idx - 1]
        out_sdf = sel_dir / f"ligand_best_pose_{pose_idx}_from_{src_sdf.stem}.sdf"
        w = Chem.SDWriter(str(out_sdf))
        try:
            w.write(m)
        finally:
            w.close()

        # optional: SMILES snapshot & selection.json (helps md_prepare discover ref SMILES)
        try:
            s = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(m)), isomericSmiles=True)
        except Exception:
            s = None
        (sel_dir / "selection.json").write_text(json.dumps({
            "source_sdf": str(src_sdf),
            "pose_index": pose_idx,
            "metric": "min_FeN_distance",
            "score": dist,
            "smiles": s,
        }, indent=2))

    _write_log("pose_select_stdout", stdout.getvalue())
    _write_log("pose_select_stderr", stderr.getvalue())

    # Assertions
    assert (sel_dir / "selection.json").exists()
    assert any(sel_dir.glob("ligand_best_pose_*_from_*.sdf"))
