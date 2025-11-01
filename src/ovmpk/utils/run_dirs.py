"""Utilities for managing per-run directory layouts."""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path
from typing import Optional

_RUN_ROOT_ENV = "OVMPK_RUN_ROOT"
_RUN_ID_ENV = "OVMPK_RUN_ID"

_run_dir_cache: Optional[Path] = None


def _compute_run_dir() -> Path:
    root = Path(os.environ.get(_RUN_ROOT_ENV, "test_run"))
    run_id = os.environ.get(_RUN_ID_ENV)
    if not run_id:
        run_id = datetime.now().strftime("run_%Y%m%d_%H%M%S")
        os.environ[_RUN_ID_ENV] = run_id
    run_dir = root / run_id
    run_dir.mkdir(parents=True, exist_ok=True)
    return run_dir


def set_run_context(run_root: Path, run_id: str) -> Path:
    """Persist the active run root/name for downstream helpers."""

    global _run_dir_cache
    os.environ[_RUN_ROOT_ENV] = str(run_root)
    os.environ[_RUN_ID_ENV] = run_id
    _run_dir_cache = (run_root / run_id)
    _run_dir_cache.mkdir(parents=True, exist_ok=True)
    return _run_dir_cache


def get_run_dir(create: bool = True) -> Path:
    """Return the active run directory (creating if necessary)."""

    global _run_dir_cache
    if _run_dir_cache is None:
        run_dir = _compute_run_dir()
        if not create and not run_dir.exists():
            return run_dir
        _run_dir_cache = run_dir
    elif create:
        _run_dir_cache.mkdir(parents=True, exist_ok=True)
    return _run_dir_cache


def stage_dir(name: str, create: bool = True) -> Path:
    """Return the path for a named stage within the current run."""

    run_dir = get_run_dir(create=create)
    path = run_dir / name
    if create:
        path.mkdir(parents=True, exist_ok=True)
    return path
