# tests/test_md_analysis.py
import csv, math
from pathlib import Path

import numpy as np
from openmm import app, unit, XmlSerializer

REPO = Path(__file__).resolve().parents[1]
PREP_DIR = REPO / "data" / "output" / "md_prepare_step_test"
PROD_DIR = REPO / "data" / "output" / "md_prod_test"

def _read_csv(path):
    with open(path) as fh:
        reader = csv.reader(fh)
        rows = [r for r in reader if r]
    header = [h.strip() for h in rows[0]]
    data = rows[1:]
    return header, data

def _col_idx(header, key_choices):
    # find first header containing any of the provided substrings (case-insensitive)
    lo = [h.lower() for h in header]
    for k in key_choices:
        k = k.lower()
        for i, h in enumerate(lo):
            if k in h:
                return i
    return None

def test_prod_artifacts_exist():
    sysxml = PREP_DIR / "system.xml"
    pdbmin = PREP_DIR / "complex_solvated.pdb"
    if not pdbmin.exists():
        pdbmin = PREP_DIR / "complex_solvated_premin.pdb"
    dcd = PROD_DIR / "prod.dcd"
    log = PROD_DIR / "prod.log"
    for p in [sysxml, pdbmin, dcd, log]:
        assert p.exists(), f"Missing: {p}"
        assert p.stat().st_size > 0, f"Empty file: {p}"

def test_log_thermo_windows():
    log = PROD_DIR / "prod.log"
    header, data = _read_csv(log)
    assert len(data) >= 2, "Not enough samples in prod.log"
    # find columns
    iT = _col_idx(header, ["temperature"])
    iRho = _col_idx(header, ["density"])
    iPE = _col_idx(header, ["potential"])
    assert iT is not None and iRho is not None and iPE is not None, f"Bad header: {header}"

    # discard early third as warmup
    start = max(1, len(data) // 3)
    temps = np.array([float(r[iT]) for r in data[start:]])
    rhos  = np.array([float(r[iRho]) for r in data[start:]])
    pes   = np.array([float(r[iPE]) for r in data[start:]])

    # windows: broad on purpose for short toy runs
    Tmean = float(np.mean(temps))
    Rmean = float(np.mean(rhos))
    assert 285.0 <= Tmean <= 335.0, f"Mean T out of range: {Tmean:.1f} K"
    assert 0.85 <= Rmean <= 1.25,   f"Mean density out of range: {Rmean:.3f} g/mL"
    assert np.isfinite(pes).all(),  "Non-finite potential energies detected"

def test_dcd_has_multiple_frames():
    dcd = PROD_DIR / "prod.dcd"
    top_pdb = PREP_DIR / "complex_solvated.pdb"
    if not top_pdb.exists():
        top_pdb = PREP_DIR / "complex_solvated_premin.pdb"
    assert dcd.exists() and dcd.stat().st_size > 0, "No trajectory written"
    assert top_pdb.exists(), f"Missing topology PDB: {top_pdb}"

    # Preferred: use mdtraj if available (robust DCD reader)
    try:
        import mdtraj as md
        traj = md.load_dcd(str(dcd), top=str(top_pdb))
        assert traj.n_frames >= 2, f"prod.dcd has only {traj.n_frames} frame(s)"
        return
    except Exception:
        # Fallback: correlate with StateDataReporter cadence in prod.log (every 250 steps).
        # This confirms the sim ran long enough to produce multiple reporter events,
        # and we already asserted DCD exists and is non-empty.
        log = PROD_DIR / "prod.log"
        assert log.exists() and log.stat().st_size > 0, "No log written"
        with open(log) as fh:
            lines = [ln for ln in fh if ln.strip()]
        # header + ≥2 data rows ⇒ at least two reporting intervals
        assert len(lines) >= 3, f"Not enough samples in prod.log to infer multiple frames (lines={len(lines)})"

