from __future__ import annotations
import shutil, subprocess, re, json
from pathlib import Path
from typing import Optional, Dict, List, Tuple

def _need(bin_name: str) -> str:
    p = shutil.which(bin_name)
    if not p:
        raise RuntimeError(f"{bin_name} not found on PATH. Install GNINA or add it to PATH.")
    return p

def _run(cmd: list[str], log_path: Optional[Path]=None) -> subprocess.CompletedProcess:
    res = subprocess.run(cmd, text=True, capture_output=True)
    if log_path:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        log_path.write_text(
            "CMD: " + " ".join(cmd) + "\nRC: " + str(res.returncode)
            + "\n--- STDOUT ---\n" + res.stdout + "\n--- STDERR ---\n" + res.stderr + "\n"
        )
    if res.returncode != 0:
        raise RuntimeError(f"Command failed ({res.returncode}): {' '.join(cmd)}\n{res.stderr}")
    return res

def gnina_available() -> Tuple[bool, Dict[str,str]]:
    try:
        _need("gnina")
        cp = subprocess.run(["gnina","--version"], text=True, capture_output=True)
        ver = (cp.stdout or cp.stderr).strip()
        info = {"version": ver}
        # Heuristic GPU presence
        info["gpu_hint"] = "yes" if shutil.which("nvidia-smi") else "no"
        return True, info
    except Exception as e:
        return False, {"error": str(e)}

def _device_flags(device: str) -> list[str]:
    device = (device or "auto").lower()
    if device == "gpu":
        return ["--gpu"]
    if device == "cpu":
        return ["--cpu"]
    # auto: try GPU if available (gnina accepts --gpu without devices? keep conservative)
    return ["--gpu"] if shutil.which("nvidia-smi") else ["--cpu"]

def dock(
    receptor_pdbqt: Path,
    ligand_pdbqt: Path,
    center: Tuple[float,float,float],
    size: Tuple[float,float,float],
    num_modes: int = 20,
    exhaustiveness: int = 8,
    seed: int = 42,
    cpu_threads: int = 1,
    device: str = "auto",
    model: Optional[str] = None,
    out_sdf: Optional[Path] = None,
    log_path: Optional[Path] = None,
) -> Path:
    _need("gnina")
    assert receptor_pdbqt.exists(), f"Missing receptor: {receptor_pdbqt}"
    assert ligand_pdbqt.exists(), f"Missing ligand: {ligand_pdbqt}"
    out_sdf = out_sdf or ligand_pdbqt.with_name(ligand_pdbqt.stem + "_gnina.sdf")
    cmd = [
        "gnina",
        "-r", str(receptor_pdbqt),
        "-l", str(ligand_pdbqt),
        "--center_x", f"{center[0]:.3f}",
        "--center_y", f"{center[1]:.3f}",
        "--center_z", f"{center[2]:.3f}",
        "--size_x", f"{size[0]:.3f}",
        "--size_y", f"{size[1]:.3f}",
        "--size_z", f"{size[2]:.3f}",
        "--num_modes", str(num_modes),
        "--exhaustiveness", str(exhaustiveness),
        "--seed", str(seed),
        "--out", str(out_sdf),
        "--cpu", str(cpu_threads),
    ] + _device_flags(device)
    if model and model != "default":
        cmd += ["--model", model]
    _run(cmd, log_path=log_path)
    if not out_sdf.exists() or out_sdf.stat().st_size == 0:
        raise RuntimeError(f"GNINA docking produced no output: {out_sdf}")
    return out_sdf

def rescore(
    poses_sdf: Path,
    receptor_pdbqt: Path,
    device: str = "auto",
    model: Optional[str] = None,
    out_sdf: Optional[Path] = None,
    log_path: Optional[Path] = None,
) -> List[Dict]:
    _need("gnina")
    assert receptor_pdbqt.exists(), f"Missing receptor: {receptor_pdbqt}"
    assert poses_sdf.exists(), f"Missing poses SDF: {poses_sdf}"
    out_sdf = out_sdf or poses_sdf.with_name(poses_sdf.stem + "_gnina_rescored.sdf")
    cmd = [
        "gnina",
        "-r", str(receptor_pdbqt),
        "-l", str(poses_sdf),
        "--score_only",
        "--out", str(out_sdf),
    ] + _device_flags(device)
    if model and model != "default":
        cmd += ["--model", model]
    _run(cmd, log_path=log_path)
    return parse_scores_from_sdf(out_sdf)

def parse_scores_from_sdf(sdf_path: Path) -> List[Dict]:
    scores: List[Dict] = []
    if not sdf_path.exists():
        return scores
    pose = 0
    current: Dict[str, float] = {}
    def flush():
        nonlocal pose, current, scores
        if current:
            current["pose"] = pose
            scores.append(current)
            current = {}
    with open(sdf_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.strip() == "$$$$":
                flush()
                pose += 1
                continue
            # SDF data fields are like:
            # >  <CNNscore>
            # 0.732
            m = re.match(r">\s*<([^>]+)>", line)
            if m:
                key = m.group(1).strip()
                val = fh.readline().strip()
                # Map known fields
                kmap = {
                    "CNNscore": "CNNscore",
                    "CNNaffinity": "CNNaffinity",
                    "minimizedAffinity": "vina_affinity",
                    "minimizedAffinity_kcal_per_mol": "vina_affinity",
                    "CNNaffinity_kcal_per_mol": "CNNaffinity",
                }
                norm = kmap.get(key)
                if norm:
                    try:
                        current[norm] = float(val)
                    except Exception:
                        pass
    flush()
    return scores
