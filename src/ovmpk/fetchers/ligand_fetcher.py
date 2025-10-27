from pathlib import Path
import os, urllib.parse, requests

PUBCHEM_SDF = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF?record_type=3d"

def fetch(name: str, cfg: dict):
    dry = bool(os.getenv("OVM_DRY_RUN", ""))
    outp = Path("data/input/ligands") / f"{name}.sdf"
    outp.parent.mkdir(parents=True, exist_ok=True)
    if dry:
        outp.write_text("SDF dry-run")
        return outp
    if not outp.exists():
        url = PUBCHEM_SDF.format(name=urllib.parse.quote(name))
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        outp.write_bytes(r.content)
    return outp
