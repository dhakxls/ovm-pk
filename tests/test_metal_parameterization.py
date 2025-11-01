"""Tests for Stage 4 metal parameterization utilities."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Any

from src.ovmpk.analysis.metal_parameterization import (
    CoordinationReport,
    apply_metal_parameterization,
    identify_coordination,
    validate_geometry,
)


def test_identify_coordination(stage2_prepared_assets: Dict[str, Any]) -> None:
    protein_prepped: Path = stage2_prepared_assets["protein_prepped"]
    report = identify_coordination(protein_prepped)
    assert isinstance(report, CoordinationReport)
    assert report.bonds, "Expected metal coordination bonds to be detected"
    assert report.all_within_tolerance


def test_apply_metal_parameterization(stage2_prepared_assets: Dict[str, Any], tmp_path: Path) -> None:
    protein_prepped: Path = stage2_prepared_assets["protein_prepped"]
    output_dir = stage2_prepared_assets["run_dir"] / "metal_parameters"
    report = apply_metal_parameterization(protein_prepped, output_dir)

    assert report.bonded_xml_path and report.bonded_xml_path.exists()
    assert report.summary_path and report.summary_path.exists()
    assert report.all_within_tolerance

    validate_geometry(report)

    metadata = report.summary_path.read_text()
    assert report.metal_site in metadata
