# tests/test_docking.py
from __future__ import annotations
import os, sys, logging, shutil, yaml
from pathlib import Path
import pytest

from ovmpk.docking.smina_wrapper import run as run_smina_docking
from ovmpk.utils.logging import get_logger

REPO = Path(__file__).resolve().parents[1]

def _find_latest(directory: Path, suffix: str, exts: list[str]) -> Path | None:
    directory = Path(directory)
    if not directory.is_dir():
        return None
    best = None
    best_mtime = -1
    for ext in exts:
        for p in directory.glob(f"*{suffix}{ext}"):
            try:
                mtime = p.stat().st_mtime
                if mtime > best_mtime:
                    best, best_mtime = p, mtime
            except FileNotFoundError:
                continue
    if best:
        return best
    # fallback: newest by extension only
    for ext in exts:
        cands = sorted(directory.glob(f"*{ext}"), key=lambda x: x.stat().st_mtime, reverse=True)
        if cands:
            return cands[0]
    return None

@pytest.mark.order(3)
def test_docking_end_to_end(run_paths):
    """
    Runs smina docking using prepared protein/ligand artifacts,
    verifies outputs, and copies final SDFs into this run's output dir.
    """
    # --- config & dirs --------------------------------------------------------
    cfg_path = REPO / "configs" / "prod_test.yaml"
    assert cfg_path.exists(), f"Missing config: {cfg_path}"
    cfg = yaml.safe_load(cfg_path.read_text())

    out_root = Path(run_paths["out_root"]) / "docking"
    out_root.mkdir(parents=True, exist_ok=True)

    log_file = Path(run_paths["log_root"]) / "test_docking.log"
    logger = get_logger("ovmpk_test_docking")
    logger.setLevel(logging.INFO)

    # Console handler (once)
    if not any(isinstance(h, logging.StreamHandler) for h in logger.handlers):
        sh = logging.StreamHandler(sys.stdout)
        sh.setLevel(logging.INFO)
        sh.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        logger.addHandler(sh)
    # File handler (fresh file each run id)
    if not any(isinstance(h, logging.FileHandler) and h.baseFilename == str(log_file) for h in logger.handlers):
        fh = logging.FileHandler(str(log_file), mode="w", encoding="utf-8")
        fh.setLevel(logging.INFO)
        fh.setFormatter(logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s"))
        logger.addHandler(fh)

    # --- locate prepared inputs ----------------------------------------------
    protein_work = REPO / "data" / "work" / "protein_prep"
    ligand_work  = REPO / "data" / "work" / "ligand_prep"

    prot_suffix = cfg.get("prep", {}).get("protein", {}).get("output_suffix", "_fixed_ph7.4")
    lig_suffix  = cfg.get("prep", {}).get("ligand", {}).get("output_suffix",  "_prepared_ph7.4")

    protein_pdb = _find_latest(protein_work, prot_suffix, [".pdb"])
    ligand_sdf  = _find_latest(ligand_work,  lig_suffix,  [".sdf"])

    assert protein_pdb and protein_pdb.exists(), (
        f"No prepared protein PDB found under {protein_work} (suffix '{prot_suffix}')."
    )
    assert ligand_sdf and ligand_sdf.exists(), (
        f"No prepared ligand SDF found under {ligand_work} (suffix '{lig_suffix}')."
    )

    # --- progress UI ----------------------------------------------------------
    use_rich = False
    try:
        from rich.progress import Progress, SpinnerColumn, BarColumn, TimeElapsedColumn, TimeRemainingColumn, TextColumn
        use_rich = True
    except Exception:
        pass

    def _progress_ctx():
        if not use_rich:
            class Dummy:
                def __enter__(self): return None
                def __exit__(self, *a): return False
                def add_task(self, *a, **k): return 0
                def update(self, *a, **k): pass
            return Dummy()
        return Progress(
            SpinnerColumn(),
            TextColumn("[bold]Docking[/bold] — {task.description}"),
            BarColumn(),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            transient=False,
        )

    logger.info("Starting docking test…")
    logger.info(f"Protein: {protein_pdb}")
    logger.info(f"Ligand:  {ligand_sdf}")

    with _progress_ctx() as prog:
        task = getattr(prog, "add_task", lambda *a, **k: 0)(
            "initializing", total=4
        )
        # Step 1: init
        getattr(prog, "update", lambda *a, **k: None)(task, description="validate inputs", advance=1)

        # Step 2: run docking
        getattr(prog, "update", lambda *a, **k: None)(task, description="running smina …")
        try:
            outputs = run_smina_docking(
                protein_pdb_in=protein_pdb,
                ligand_sdf_in=ligand_sdf,
                cfg=cfg
            )
        finally:
            getattr(prog, "update", lambda *a, **k: None)(task, advance=1)

        # Step 3: verify outputs
        getattr(prog, "update", lambda *a, **k: None)(task, description="verifying outputs")
        assert outputs and isinstance(outputs, list), "Docking returned no output paths."
        missing = [p for p in outputs if not Path(p).exists() or Path(p).stat().st_size == 0]
        if missing:
            logger.error("Docking outputs missing/empty: %s", missing)
        assert not missing, "Some docking outputs are missing or empty."
        getattr(prog, "update", lambda *a, **k: None)(task, advance=1)

        # Step 4: collect copies
        getattr(prog, "update", lambda *a, **k: None)(task, description=f"collecting to {out_root}")
        for p in outputs:
            dst = out_root / Path(p).name
            try:
                shutil.copy2(p, dst)
            except Exception as e:
                logger.warning("Failed to copy %s → %s: %s", p, dst, e)
        getattr(prog, "update", lambda *a, **k: None)(task, advance=1)

    # Final assertions / breadcrumbs
    copied = list(out_root.glob("*.sdf"))
    assert copied, f"No copied docking SDFs in {out_root}"
    logger.info("Docking test finished. Collected %d file(s) in %s", len(copied), out_root)

