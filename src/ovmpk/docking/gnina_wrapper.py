from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

try:
    from ovmpk.utils.logging import get_logger  # type: ignore
    logger = get_logger("gnina_wrapper")
except Exception:  # pragma: no cover
    import logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("gnina_wrapper")


def has_gnina() -> bool:
    """Return True if the 'gnina' binary is available on PATH."""
    return shutil.which("gnina") is not None


def build_gnina_rescore_cmd(
    receptor_pdbqt: Path | str,
    poses_path: Path | str,
    cnn_scoring: str = "ad4",
    cnn_models: Optional[str] = None,
    extra_args: Optional[Iterable[str]] = None,
) -> List[str]:
    """
    Build a GNINA command that *rescoring-only* of an existing pose file.

    Notes:
    - We pass --score_only so GNINA evaluates and prints CNN/Vina scores
      for the supplied poses rather than docking.
    - We intentionally avoid GPU flags here (GNINA auto-selects if built with CUDA).
    """
    receptor_pdbqt = str(receptor_pdbqt)
    poses_path = str(poses_path)

    cmd = [
        "gnina",
        "-r", receptor_pdbqt,
        "-l", poses_path,
        "--score_only",
        "--cnn_scoring", cnn_scoring,
        "--verbosity", "1",
    ]
    if cnn_models:
        cmd += ["--cnn_models", cnn_models]
    if extra_args:
        cmd += list(extra_args)
    return cmd


def rescore(
    receptor_pdbqt: Path | str,
    poses_path: Path | str,
    out_csv: Path | str,
    cnn_scoring: str = "ad4",
    cnn_models: Optional[str] = None,
    extra_args: Optional[Iterable[str]] = None,
    dry_run: bool = False,
) -> Tuple[int, List[str]]:
    """
    Run GNINA in score-only mode and write a simple CSV with raw stdout lines.
    Returns (returncode, command_list). If dry_run=True, no subprocess is executed.
    """
    cmd = build_gnina_rescore_cmd(
        receptor_pdbqt=receptor_pdbqt,
        poses_path=poses_path,
        cnn_scoring=cnn_scoring,
        cnn_models=cnn_models,
        extra_args=extra_args,
    )

    if dry_run:
        logger.info("GNINA dry-run command: %s", " ".join(cmd))
        return 0, cmd

    if not has_gnina():
        raise FileNotFoundError("gnina not found on PATH. Install gnina or adjust PATH.")

    logger.info("Running GNINA rescoring:\n  %s", " ".join(cmd))
    proc = subprocess.run(cmd, text=True, capture_output=True)
    rc = proc.returncode
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)

    # Minimal CSV: keep raw stdout for auditing; downstream parse can happen elsewhere
    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("line\n")
        for line in (proc.stdout or "").splitlines():
            # Keep everything; parsing (e.g., 'CNNscore', 'CNNaffinity') can be added later
            f.write(f"{line.replace(',', ';')}\n")

    if rc != 0:
        logger.error("GNINA failed (rc=%d). Stderr tail:\n%s", rc, "\n".join(proc.stderr.splitlines()[-40:]))

    return rc, cmd
