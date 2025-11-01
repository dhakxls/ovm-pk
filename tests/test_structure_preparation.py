"""Stage 2 tests: structure preparation built on Stage 1 fetching outputs."""
from __future__ import annotations

from pathlib import Path
from typing import Dict

import pytest

from src.ovmpk.prep import ligand_prep, protein_prep


@pytest.mark.skipif(not protein_prep.HAS_PDBFIXER, reason="PDBFixer is required for protein preparation")
def test_prepare_protein_creates_fixed_structure(stage1_assets: Dict[str, Path]) -> None:
    """Validate protein preparation generates a cleaned PDB suitable for downstream steps."""
    input_path: Path = stage1_assets["protein_path"]
    output_dir = stage1_assets["run_dir"] / "protein_prep"
    cfg = {
        "prep": {
            "protein": {
                "ph": 7.4,
                "output_suffix": "_stage2_prepped",
                "run_pdbfixer": True,
                "output_dir": str(output_dir),
            }
        }
    }

    expected_output = output_dir / f"{input_path.stem}_stage2_prepped.pdb"
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
def test_prepare_ligand_generates_protonated_sdf(stage1_assets: Dict[str, Path]) -> None:
    """Validate ligand preparation adds hydrogens/charges and writes new SDF."""
    input_path: Path = stage1_assets["ligand_path"]
    output_dir = stage1_assets["run_dir"] / "ligand_prep"
    cfg = {
        "prep": {
            "ligand": {
                "ph": 7.4,
                "charge_method": "gasteiger",
                "output_suffix": "_stage2_prepped",
                "force_3d": True,
                "use_dimorphite": True,
                "output_dir": str(output_dir),
            }
        }
    }

    expected_output = output_dir / f"{input_path.stem}_stage2_prepped.sdf"
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
