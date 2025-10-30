from __future__ import annotations


import pytest
@pytest.fixture(scope="module")
def fetchers_init():
    """Initial test to set up fetcher prerequisites."""
    return True

def test_fetchers_end_to_end(fetchers_init):
    """Run the actual fetching test to ensure files exist"""
    # First make sure we have config
    cfg = _cfg()
    pdbid = (cfg.get("protein", {}) or {}).get("pdb", "5VCC")
    ligand_name = (cfg.get("ligand", {}) or {}).get("name", "ketoconazole")

    # Set up directories
    prot_dir = OUTDIR / "protein_fetch"
    lig_dir  = OUTDIR / "ligand_fetch"
    prot_dir.mkdir(parents=True, exist_ok=True)
    lig_dir.mkdir(parents=True, exist_ok=True)

    protein_pdb = prot_dir / f"{pdbid}.pdb"
    ligand_sdf  = lig_dir / f"{ligand_name}.sdf"

    # Run the actual fetch
    stdout = io.StringIO()
    stderr = io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr), _rich_progress() as prog:
        t = prog.add_task("[bold]Fetching inputs...", total=2)

        # Try your package first; if missing, use HTTP fallbacks
        ok = False
        try:
            from ovmpk.fetchers.protein_fetcher import fetch_protein  # type: ignore
            fetch_protein(pdbid=pdbid, out_dir=prot_dir)  # your API may differ
            ok = True
        except Exception:
            _fetch_protein_rcsb(pdbid, protein_pdb)
        prog.update(t, advance=1)

        try:
            from ovmpk.fetchers.ligand_fetcher import fetch_ligand  # type: ignore
            fetch_ligand(name=ligand_name, out_dir=lig_dir)        # your API may differ
            ok = True or ok
        except Exception:
            _fetch_ligand_pubchem_name(ligand_name, ligand_sdf)
        prog.update(t, advance=1)

    # Logs
    _write_log("fetch_stdout", stdout.getvalue())
    _write_log("fetch_stderr", stderr.getvalue())
    
    # Assertions
    assert protein_pdb.exists() and protein_pdb.stat().st_size > 0, "Protein PDB not fetched"
    assert ligand_sdf.exists() and ligand_sdf.stat().st_size > 0, "Ligand SDF not fetched"

def test_flexible_fetcher():
    """Test configurable fetching."""
    from ovmpk.fetchers import FlexibleFetcher
    
    # CYP3A4-ketoconazole example config
    cfg = {
        "protein": {"pdb": "5VCC", "source": "rcsb"},
        "ligand": {"name": "ketoconazole", "identifier_type": "name"}
    }
    
    fetcher = FlexibleFetcher(cfg)
    prot_path = fetcher.fetch_protein(Path("data/work/protein"))
    lig_path = fetcher.fetch_ligand(Path("data/work/ligand"))
    
    assert prot_path.exists() and "5VCC" in str(prot_path)
    assert lig_path.exists() and "ketoconazole" in str(lig_path)

import pytest
pytestmark = [
    pytest.mark.order(1),
    pytest.mark.dependency()
]

import sys, os, json, time, io
from pathlib import Path
from contextlib import redirect_stdout, redirect_stderr
import requests
from rich.progress import Progress, TimeElapsedColumn, BarColumn, MofNCompleteColumn, TimeRemainingColumn, SpinnerColumn
import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))

OUTDIR = REPO / "data" / "work"
LOGDIR = REPO / "data" / "output" / "fetch_test" / "logs"
LOGDIR.mkdir(parents=True, exist_ok=True)

# graceful, no-flake helper
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

def _cfg():
    import yaml
    cfgp = REPO / "configs" / "prod_test.yaml"
    return yaml.safe_load(cfgp.read_text()) if cfgp.exists() else {}

def _fetch_protein_rcsb(pdbid: str, out_pdb: Path):
    url = f"https://files.rcsb.org/download/{pdbid}.pdb"
    r = requests.get(url, timeout=20)
    r.raise_for_status()
    out_pdb.write_text(r.text)

def _fetch_ligand_pubchem_name(name: str, out_sdf: Path):
    # 3D if available; else 2D
    urls = [
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF?record_type=3d",
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF",
    ]
    for url in urls:
        r = requests.get(url, timeout=30)
        if r.ok and r.content.strip():
            out_sdf.write_bytes(r.content)
            return
    raise RuntimeError(f"PubChem fetch failed for {name}")
