"""Input configuration utilities for automated ligand/enzyme curation."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Mapping, Optional

import yaml


def _default_run_name() -> str:
    return datetime.now().strftime("run_%Y%m%d_%H%M%S")


@dataclass
class AutoCurateConfig:
    """Normalized configuration for the auto-curation workflow."""

    protein_input: str
    ligand_input: str
    run_root: Path = Path("test_run")
    run_name: Optional[str] = None
    fetch: Dict[str, Any] = field(default_factory=dict)
    prep: Dict[str, Any] = field(default_factory=dict)
    ligand: Dict[str, Any] = field(default_factory=dict)
    metal: Dict[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.run_name:
            self.run_name = _default_run_name()
        self.run_root = Path(self.run_root)

    @property
    def run_dir(self) -> Path:
        return self.run_root / Path(self.run_name)


def _coerce_mapping(config: Any, *, source: str) -> Mapping[str, Any]:
    if isinstance(config, Mapping):
        return config
    raise TypeError(f"Expected mapping for {source}, got {type(config)!r}")


def load_config(source: Optional[Any] = None, **overrides: Any) -> AutoCurateConfig:
    """Load auto-curation config from YAML file, mapping, or kwargs."""

    data: Dict[str, Any] = {}

    if isinstance(source, (str, Path)):
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(f"Auto-curation config not found: {path}")
        data = yaml.safe_load(path.read_text()) or {}
    elif source is not None:
        data = dict(_coerce_mapping(source, source="config"))

    if overrides:
        data.update(overrides)

    required = {"protein_input", "ligand_input"}
    missing = [key for key in required if key not in data]
    if missing:
        raise ValueError(f"Auto-curation config missing required keys: {', '.join(missing)}")

    return AutoCurateConfig(
        protein_input=data["protein_input"],
        ligand_input=data["ligand_input"],
        run_root=Path(data.get("run_root", "test_run")),
        run_name=data.get("run_name"),
        fetch=dict(data.get("fetch", {})),
        prep=dict(data.get("prep", {})),
        ligand=dict(data.get("ligand", {})),
        metal=dict(data.get("metal", {})),
    )
