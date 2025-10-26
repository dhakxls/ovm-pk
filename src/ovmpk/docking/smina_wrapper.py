import os, shutil
from pathlib import Path

def run(protein_pdb: Path, ligand_sdf: Path, cfg: dict):
    # If smina is available and dry-run is off, you would invoke it here.
    # For MVP/CI: return a dummy "poses" list.
    dry = bool(os.getenv("OVM_DRY_RUN", "")) or shutil.which("smina") is None
    out = Path("results") / "poses.sdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    if dry:
        out.write_text("SDF pose - dry run")
        return [out.as_posix()]
    # TODO: real smina command invocation
    return [out.as_posix()]
