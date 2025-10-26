from pathlib import Path
import os
import requests

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"

def _dry_file(path: Path, contents: str = "PDB"):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(contents)

def fetch(gene: str, cfg: dict):
    dry = bool(os.getenv("OVM_DRY_RUN", ""))
    base = Path("data/input/proteins")
    base.mkdir(parents=True, exist_ok=True)
    # Minimal mapping for MVP; real logic will search PDB by gene and policy
    pdb_id = "1TQN"  # CYP3A4 apo (example); for dry-run we just create a file
    outp = base / f"{pdb_id}_apo.pdb"
    if dry:
        _dry_file(outp, contents="REMARK dry-run PDB")
    else:
        if not outp.exists():
            r = requests.get(RCSB_URL.format(pdb_id=pdb_id), timeout=30)
            r.raise_for_status()
            outp.write_text(r.text)
    return {"apo": outp}
