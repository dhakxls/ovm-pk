from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Tuple

from ovmpk.docking.gnina_wrapper import has_gnina, rescore
from ovmpk.utils.logging import get_logger

logger = get_logger("gnina_engine")


class GninaEngine:
    """
    Lightweight "engine" that rescored existing SDF/PDBQT poses using GNINA's CNN.
    Intended to be used *after* pose generation with smina or another engine.
    """

    def __init__(self, cfg: Dict):
        self.cfg = cfg or {}
        self.enabled = bool(self.cfg.get("docking", {}).get("gnina", {}).get("enabled", False))
        self.cnn_scoring = self.cfg.get("docking", {}).get("gnina", {}).get("cnn_scoring", "ad4")
        self.cnn_models = self.cfg.get("docking", {}).get("gnina", {}).get("cnn_models", None)
        self.extra_args = self.cfg.get("docking", {}).get("gnina", {}).get("extra_args", [])

    def available(self) -> bool:
        return self.enabled and has_gnina()

    def rescore_file(
        self,
        receptor_pdbqt: Path,
        poses_path: Path,
        out_csv: Path,
        dry_run: bool = False,
    ) -> Tuple[int, list]:
        """
        Rescore poses via GNINA; writes a CSV containing GNINA stdout lines.
        """
        if not self.enabled:
            logger.info("GNINA rescoring is disabled in config; skipping.")
            return 0, []

        rc, cmd = rescore(
            receptor_pdbqt=receptor_pdbqt,
            poses_path=poses_path,
            out_csv=out_csv,
            cnn_scoring=self.cnn_scoring,
            cnn_models=self.cnn_models,
            extra_args=self.extra_args,
            dry_run=dry_run,
        )
        return rc, cmd
