"""Ligand parameter generation and force-field bundling utilities."""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional


DEFAULT_TEMPLATE = Path("forcefields/shahrokh_heme_ic6_unl.ffxml")


@dataclass
class LigandParameterBundle:
    """Paths to generated force-field artifacts for a ligand."""

    ffxml: Path
    metadata: Path


def build_forcefield_bundle(
    ligand_prepped: Path,
    output_dir: Path,
    config: Optional[Dict[str, Any]] = None,
) -> LigandParameterBundle:
    """Generate or assemble force-field XML for the prepared ligand.

    Current implementation copies an existing template (e.g., Shahrokh heme +
    ligand ffxml) into the run assets directory so that downstream stages can
    use a run-scoped copy. Future iterations should invoke AmberTools/ParmEd to
    derive ligand parameters automatically when the tooling is available.
    """

    cfg = dict(config or {})

    source_template = Path(cfg.get("ffxml_template", DEFAULT_TEMPLATE))
    if not source_template.exists():
        raise FileNotFoundError(f"Force-field template not found: {source_template}")

    output_dir.mkdir(parents=True, exist_ok=True)
    bundle_name = cfg.get("bundle_name") or f"{ligand_prepped.stem}_combined.ffxml"
    ffxml_path = output_dir / bundle_name
    shutil.copyfile(source_template, ffxml_path)

    metadata = {
        "ligand_source": str(ligand_prepped),
        "template_ffxml": str(source_template),
        "notes": cfg.get("notes", "Copied existing template; no reparameterization performed."),
    }

    metadata_path = output_dir / "ffxml_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    return LigandParameterBundle(ffxml=ffxml_path, metadata=metadata_path)
