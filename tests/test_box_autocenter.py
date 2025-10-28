# tests/test_box_autocenter.py
from pathlib import Path
import math

REPO = Path(__file__).resolve().parents[1]

# Import the helper directly from your wrapper
from ovmpk.docking.smina_wrapper import autocenter_box_from_pdb

def _mean_fe_from_pdb(pdb: Path):
    """Parse HEM/FE positions (PDB columns) and return mean (x,y,z)."""
    fes = []
    with pdb.open("r") as fh:
        for line in fh:
            if not (line.startswith("HETATM") or line.startswith("ATOM")):
                continue
            resname = line[17:20].strip()
            atom = line[12:16].strip()
            if atom.upper().startswith("FE") and resname in {"HEM", "HEC", "HEA", "HEB"}:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                fes.append((x, y, z))
    if not fes:
        return None
    n = len(fes)
    return (sum(x for x,_,_ in fes)/n, sum(y for _,y,_ in fes)/n, sum(z for *_,z in fes)/n)

def _dist(a, b):
    return math.sqrt(sum((ai - bi)**2 for ai, bi in zip(a, b)))

def test_autocenter_matches_heme_fe_centroid():
    # Use latest prepared receptor from your pipeline
    pdb_dir = REPO / "data" / "work" / "protein_prep"
    pdbs = sorted(pdb_dir.glob("*.pdb"), key=lambda p: p.stat().st_mtime, reverse=True)
    assert pdbs, f"No PDBs found in {pdb_dir}"
    pdb = pdbs[0]

    expected = _mean_fe_from_pdb(pdb)
    assert expected is not None, "No HEM/FE found in receptor PDB — this test assumes a heme target."

    center = autocenter_box_from_pdb(pdb, auto_center="heme_fe", fallback_center=(-999,-999,-999))
    assert center is not None and len(center) == 3
    assert _dist(center, expected) < 1e-3, f"Center {center} != FE mean {expected}"

def test_autocenter_fallback_when_no_metal(tmp_path):
    # Create a tiny non-heme PDB (alpha carbon only)
    fake = tmp_path / "mini.pdb"
    fake.write_text(
        "ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00  0.00           C  \n"
        "END\n"
    )
    fb = (1.0, 2.0, 3.0)
    center = autocenter_box_from_pdb(fake, auto_center="heme_fe", fallback_center=fb)
    assert tuple(center) == fb, "Fallback center not respected when no heme present."
