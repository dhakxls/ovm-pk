from __future__ import annotations
import subprocess, shutil
from pathlib import Path

def _need(bin_name: str) -> str:
    p = shutil.which(bin_name)
    if not p:
        raise RuntimeError(f"Required binary not found on PATH: {bin_name}")
    return p

def _run(cmd: list[str], log_path: Path|None=None) -> None:
    res = subprocess.run(cmd, text=True, capture_output=True)
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            "CMD: " + " ".join(cmd) + "\nRC: " + str(res.returncode)
            + "\n--- STDOUT ---\n" + res.stdout + "\n--- STDERR ---\n" + res.stderr + "\n"
        )
    if res.returncode != 0:
        raise RuntimeError(f"Command failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr}")

def to_pdbqt_receptor(inp_pdb: Path, out_pdbqt: Path) -> None:
    _need("obabel")
    # Try richer flags first, then degrade.
    logsuf = out_pdbqt.with_suffix(".obabel_receptor.log")
    for args in (["-h","-p","-a"], ["-h","-a"], ["-h"]):
        try:
            _run(["obabel","-ipdb",str(inp_pdb),"-opdbqt","-O",str(out_pdbqt), *args],
                 log_path=logsuf)
            if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
                return
        except Exception:
            continue
    raise RuntimeError(f"Failed to convert receptor PDB→PDBQT: {inp_pdb} -> {out_pdbqt}")

def to_pdbqt_ligand(inp_sdf: Path, out_pdbqt: Path) -> None:
    _need("obabel")
    log = out_pdbqt.with_suffix(".obabel_ligand.log")
    _run(["obabel","-isdf",str(inp_sdf),"-opdbqt","-O",str(out_pdbqt),"-h","--partialcharge","gasteiger"],
         log_path=log)
    if not (out_pdbqt.exists() and out_pdbqt.stat().st_size>0):
        raise RuntimeError(f"Failed to convert ligand SDF→PDBQT: {inp_sdf} -> {out_pdbqt}")
