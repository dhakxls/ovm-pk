"""Stage 3 tests: system setup built on Stage 2 prepared outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Dict, Any

import pytest

from src.ovmpk.prep import ligand_prep, protein_prep
from tests.conftest import run_stage3_system_setup

HAS_PDBFIXER = protein_prep.HAS_PDBFIXER
HAS_RDKIT = ligand_prep.HAS_RDKIT


@pytest.mark.skipif(not HAS_PDBFIXER, reason="PDBFixer is required for protein preparation")
@pytest.mark.skipif(not HAS_RDKIT, reason="RDKit is required for ligand preparation")
def test_system_setup_generates_solvated_neutral_system(
    stage2_prepared_assets: Dict[str, Any], tmp_path: Path
) -> None:
    """Validate system setup adds solvent, ions, and applies target force fields."""

    result = run_stage3_system_setup(stage2_prepared_assets, tmp_path)

    assert result["total_atoms"] > result["initial_atoms"], "Expected solvent to increase atom count"
    assert result["water_resnames"], "Expected water residues after solvation"
    assert result["ion_resnames"], "Expected neutralizing ions after system setup"

    assert result["system"].getNumParticles() == result["total_atoms"]
    assert result["system"].getNumForces() > 0

    solvated_pdb = result["solvated_pdb"]
    system_xml = result["system_xml"]
    assert solvated_pdb.exists() and solvated_pdb.stat().st_size > 0
    assert system_xml.exists() and system_xml.stat().st_size > 0
