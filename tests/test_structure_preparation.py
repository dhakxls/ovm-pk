"""Stage 2 tests: structure preparation built on Stage 1 fetching outputs."""
from __future__ import annotations

import os
from pathlib import Path

import pytest

# Ensure Stage 1 fixture is registered for this module
from tests.test_data_fetching import test_config  # noqa: F401

from src.ovmpk.prep import ligand_prep, protein_prep


@pytest.fixture(scope="module")
def stage1_assets(test_config: dict) -> dict:
    """Reuse Stage 1 outputs (protein/ligand downloads) for downstream prep tests."""
    override_dir = os.environ.get("OVMPK_STAGE1_RUN_DIR")
    run_dir: Path = Path(override_dir) if override_dir else test_config["run_dir"]
    protein_file = Path(run_dir / f"{test_config['pdb_id']}.pdb")
    ligand_file = Path(run_dir / f"{test_config['ligand_name']}.sdf")

    if not protein_file.exists():
        protein_file = Path(test_config["fetcher"].fetch_protein(test_config["pdb_id"]))
    if not ligand_file.exists():
        ligand_file = Path(test_config["fetcher"].fetch_ligand(test_config["ligand_name"]))

    return {
        **test_config,
        "protein_path": protein_file,
        "ligand_path": ligand_file,
    }


@pytest.mark.skipif(not protein_prep.HAS_PDBFIXER, reason="PDBFixer is required for protein preparation")
def test_prepare_protein_creates_fixed_structure(stage1_assets: dict) -> None:
    """Validate protein preparation generates a cleaned PDB suitable for downstream steps."""
    input_path: Path = stage1_assets["protein_path"]
    cfg = {
        "prep": {
            "protein": {
                "ph": 7.4,
                "output_suffix": "_stage2_prepped",
                "run_pdbfixer": True,
            }
        }
    }

    expected_output = protein_prep.WORK_DIR / f"{input_path.stem}_stage2_prepped.pdb"
    expected_output.unlink(missing_ok=True)

    prepared_path = protein_prep.prepare({"apo": input_path}, cfg)

    assert prepared_path.exists()
    assert prepared_path == expected_output
    assert prepared_path.stat().st_size > 0

    from openmm import app

    parsed = app.PDBFile(str(prepared_path))
    assert parsed.topology.getNumAtoms() > 0
    assert any(
        res.name in {"HEM", "HEC", "HEM"} for chain in parsed.topology.chains() for res in chain.residues()
    ), "Expected heme cofactor to be retained after preparation"


@pytest.mark.skipif(not ligand_prep.HAS_RDKIT, reason="RDKit is required for ligand preparation")
def test_prepare_ligand_generates_protonated_sdf(stage1_assets: dict) -> None:
    """Validate ligand preparation adds hydrogens/charges and writes new SDF."""
    input_path: Path = stage1_assets["ligand_path"]
    cfg = {
        "prep": {
            "ligand": {
                "ph": 7.4,
                "charge_method": "gasteiger",
                "output_suffix": "_stage2_prepped",
                "force_3d": True,
            }
        }
    }

    expected_output = ligand_prep.WORK_DIR / f"{input_path.stem}_stage2_prepped.sdf"
    expected_output.unlink(missing_ok=True)

    prepared_path = ligand_prep.prepare(input_path, cfg)

    assert prepared_path.exists()
    assert prepared_path == expected_output
    assert prepared_path.stat().st_size > 0

    content = prepared_path.read_text()
    assert "V3000" in content
    assert "$$$$" in content

    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(prepared_path), removeHs=False)
    mol = supplier[0]
    assert mol is not None
    assert mol.GetNumAtoms() >= stage1_assets.get("expected_ligand_atom_floor", 20)
""
