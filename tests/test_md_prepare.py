# tests/test_md_prepare.py

import logging
import os
import sys
import json
from pathlib import Path
import subprocess
import yaml
import hashlib

# Optional: RDKit for pulling a reference SMILES from prepared ligand SDF
try:
    from rdkit import Chem
    HAS_RDKIT = True
except Exception:
    HAS_RDKIT = False

from pathlib import Path
import sys
REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from ovmpk.utils.logging import get_logger
from ovmpk.utils.smiles_fetch import get_trusted_smiles

REPO = Path(__file__).resolve().parents[1]
MD_SCRIPT = REPO / "scripts" / "md_gpu_prepare.py"

STAGES = ["ligand", "receptor", "merge", "forcefield", "solvate", "system", "minimize", "all"]


def sha1sum(path: Path) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def _latest(path_glob):
    """Return the most recently modified Path matching the glob, or None if none exist."""
    paths = list(path_glob)
    if not paths:
        return None
    return max(paths, key=lambda p: p.stat().st_mtime)


def _get_ref_smiles_from_prepared_sdf() -> str | None:
    """
    Find the newest *_prepared_*.sdf under data/work/ligand_prep/ and return a non-isomeric SMILES.

    We prefer non-isomeric canonical SMILES as a robust graph template for fixing bond orders
    (protonation/tautomer handling is done inside md_gpu_prepare.py).
    """
    if not HAS_RDKIT:
        return None
    lig_dir = REPO / "data" / "work" / "ligand_prep"
    if not lig_dir.exists():
        return None
    sdf = _latest(lig_dir.glob("*_prepared_*.sdf"))
    if not sdf:
        return None
    suppl = Chem.SDMolSupplier(str(sdf), removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        return None
    try:
        # Use heavy-atom graph only for template; avoid stereo/protons here
        smiles = Chem.MolToSmiles(Chem.RemoveHs(Chem.Mol(mol)), isomericSmiles=False)
        return smiles or None
    except Exception:
        return None


def get_reference_smiles(cfg: dict) -> str | None:
    # 1) env/config (keeps prior behavior)
    s = os.environ.get("OVMPK_LIGAND_SMILES")
    if s:
        return s
    s = cfg.get("ligand", {}).get("reference_smiles")
    if s:
        return s

    # 2) selection.json if present (pose selection might have written it)
    sel_json = REPO / "data" / "output" / "pose_selection_test" / "selected_pose" / "selection.json"
    if sel_json.exists():
        try:
            data = json.loads(sel_json.read_text())
            if data.get("smiles"):
                return data["smiles"]
        except Exception:
            pass

    # 3) Auto-resolve from the same SDF we’ll pass to md_gpu_prepare (hands-off path)
    sdf_dir = REPO / "data" / "output" / "pose_selection_test" / "selected_pose"
    candidates = [p for p in sdf_dir.glob("*.sdf") if not p.name.endswith(".deradicalized.sdf")]
    if candidates:
        sdf = max(candidates, key=lambda p: p.stat().st_mtime)
        try:
            smiles, source = get_trusted_smiles(sdf, pose_index=1)  # returns (smiles, "PubChem"/"ChEMBL"/etc.)
            if smiles:
                print(f"[info] Auto-fetched SMILES from {source}: {smiles}")
                return smiles
        except Exception:
            # ok to fall through; md_gpu_prepare can also --auto-ligand-smiles
            pass

    # 4) As a last resort, fall back to a prepared SDF (if any)
    return _get_ref_smiles_from_prepared_sdf()


def run_stage(
    stage: str,
    receptor: Path,
    sdf: Path,
    outdir: Path,
    cfg: dict,
    pose_index: int = 1,
    ref_smiles: str | None = None,
) -> None:
    """
    Invoke md_gpu_prepare.py for a single stage and assert it completes.
    Writes a detailed {stage}.log with provenance, cmd, RC, and stdout/stderr.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, str(MD_SCRIPT),
        "--receptor", str(receptor),
        "--sdf", str(sdf),
        "--pose-index", str(pose_index),
        "--out", str(outdir),
        "--ph", str(cfg.get("prep", {}).get("protein", {}).get("ph", 7.4)),
        "--padding", str(cfg.get("md", {}).get("padding", 10.0)),
        "--ionic", str(cfg.get("md", {}).get("ionic", 0.15)),
        "--water-model", "tip3p",
        "--strip-heme",
        "--strip-resnames", "GOL,SO4,EDO,PEG,TRS,ACT,PG4,PO4,CL,NA",
        "--platform", "CUDA",
        "--stop-after", stage,
    ]
    if ref_smiles:
        cmd += ["--ligand-smiles", ref_smiles]
    else:
        # exercise the script’s hands-off resolver path
        cmd += ["--auto-ligand-smiles"]

    # Include script SHA to guarantee we invoked the expected file
    script_sha = sha1sum(MD_SCRIPT)[:12]
    log = outdir / f"{stage}.log"

    # Run and capture output
    res = subprocess.run(cmd, text=True, capture_output=True)
    log.write_text(
        f"SCRIPT: {MD_SCRIPT} sha1={script_sha}\n"
        f"CMD: {' '.join(cmd)}\nRC: {res.returncode}\n"
        f"--- STDOUT ---\n{res.stdout}\n--- STDERR ---\n{res.stderr}\n"
    )

    if res.returncode != 0:
        tail_out = "\n".join(res.stdout.strip().splitlines()[-40:])
        tail_err = "\n".join(res.stderr.strip().splitlines()[-40:])
        raise AssertionError(
            f"{stage} failed (rc={res.returncode})\n"
            f"script sha1={script_sha}\n"
            f"stdout tail:\n{tail_out}\n\nstderr tail:\n{tail_err}\n\nSee log: {log}"
        )

    ok_marker = outdir / f"{stage}.ok"
    assert ok_marker.exists(), f"Missing stage marker: {ok_marker}"


def test_md_prepare_step_by_step(tmp_path: Path):
    """
    End-to-end stepwise test of scripts/md_gpu_prepare.py:
      ligand → receptor → merge → forcefield → solvate → system → minimize → all
    """
    logger = get_logger("ovmpk_test_md_step")
    logger.setLevel(logging.INFO)

    cfg_path = REPO / "configs" / "prod_test.yaml"
    cfg = yaml.safe_load(cfg_path.read_text())

    # Inputs (most recent artifacts from earlier tests)
    receptor_dir = REPO / "data" / "work" / "protein_prep"

    sdf_dir = REPO / "data" / "output" / "pose_selection_test" / "selected_pose"
    candidates = [p for p in sdf_dir.glob("*.sdf") if not p.name.endswith(".deradicalized.sdf")]
    assert candidates, f"No selected-pose SDFs found (excluding *.deradicalized.sdf) in {sdf_dir}"
    sdf = max(candidates, key=lambda p: p.stat().st_mtime)

    receptor = _latest(receptor_dir.glob("*.pdb"))
    assert receptor and receptor.exists(), f"No receptor PDB found in {receptor_dir}"

    assert sdf and sdf.exists(), f"No selected pose SDF found in {sdf_dir}"

    outdir = REPO / "data" / "output" / "md_prepare_step_test"

    # Clean output directory (except keep directory itself)
    if outdir.exists():
        for p in outdir.glob("*"):
            try:
                p.unlink()
            except IsADirectoryError:
                pass
    outdir.mkdir(parents=True, exist_ok=True)

    # Determine a trusted ligand SMILES (env/config/selection.json → resolver)
    ref_smiles = get_reference_smiles(cfg)
    if ref_smiles:
        logger.info(f"Using reference SMILES for MD prep: {ref_smiles}")
    else:
        logger.warning("No reference SMILES available; continuing without --ligand-smiles (will use auto-resolver).")

    # Run through all stages
    run_stage("ligand", receptor, sdf, outdir, cfg, ref_smiles=ref_smiles)
    assert (outdir / "ligand_pose.pdb").exists()
    assert (outdir / "ligand_smiles.txt").exists()

    run_stage("receptor", receptor, sdf, outdir, cfg, ref_smiles=ref_smiles)
    assert (outdir / "receptor_fixed.pdb").exists()

    run_stage("merge", receptor, sdf, outdir, cfg, ref_smiles=ref_smiles)
    assert (outdir / "complex_merged.pdb").exists()

    run_stage("forcefield", receptor, sdf, outdir, cfg, ref_smiles=ref_smiles)
    assert (outdir / "ff.info.json").exists()

    run_stage("solvate", receptor, sdf, outdir, cfg, ref_smiles=ref_smiles)
    assert (outdir / "complex_solvated_premin.pdb").exists() or (outdir / "complex_solvated.pdb").exists()


    run_stage("system", receptor, sdf, outdir, cfg, ref_smiles=ref_smiles)
    assert (outdir / "system.xml").exists()

    run_stage("minimize", receptor, sdf, outdir, cfg, ref_smiles=ref_smiles)
    assert (outdir / "integrator_min.xml").exists()
    assert (outdir / "state_min.xml").exists()
    assert (outdir / "complex_solvated.pdb").exists()

    run_stage("all", receptor, sdf, outdir, cfg, ref_smiles=ref_smiles)
    assert (outdir / "prepare.json").exists()
