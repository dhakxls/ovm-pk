"""Stage 3 tests: system setup built on Stage 2 prepared outputs."""
from __future__ import annotations

from pathlib import Path

import pytest

from tests.test_data_fetching import test_config  # noqa: F401
from src.ovmpk.prep import ligand_prep, protein_prep

try:
    from openmm import XmlSerializer, unit
    from openmm import app
    HAS_OPENMM = True
except ImportError:  # pragma: no cover - environment dependent
    HAS_OPENMM = False

HAS_PDBFIXER = protein_prep.HAS_PDBFIXER
HAS_RDKIT = ligand_prep.HAS_RDKIT


@pytest.fixture(scope="module")
def stage2_prepared_assets(test_config: dict) -> dict:
    """Run Stage 2 preparation to produce inputs for system setup."""
    # Ensure Stage 1 outputs are available
    protein_raw = Path(test_config["fetcher"].fetch_protein(test_config["pdb_id"]))
    ligand_raw = Path(test_config["fetcher"].fetch_ligand(test_config["ligand_name"]))

    prot_cfg = {
        "prep": {
            "protein": {
                "ph": 7.4,
                "output_suffix": "_stage2_prepped",
                "run_pdbfixer": True,
            }
        }
    }
    ligand_cfg = {
        "prep": {
            "ligand": {
                "ph": 7.4,
                "charge_method": "gasteiger",
                "output_suffix": "_stage2_prepped",
                "force_3d": True,
            }
        }
    }

    protein_prepped = protein_prep.prepare({"apo": protein_raw}, prot_cfg)
    ligand_prepped = ligand_prep.prepare(ligand_raw, ligand_cfg)

    return {
        **test_config,
        "protein_prepped": protein_prepped,
        "ligand_prepped": ligand_prepped,
    }


@pytest.mark.skipif(not HAS_OPENMM, reason="OpenMM is required for system setup tests")
@pytest.mark.skipif(not HAS_PDBFIXER, reason="PDBFixer is required for protein preparation")
@pytest.mark.skipif(not HAS_RDKIT, reason="RDKit is required for ligand preparation")
def test_system_setup_generates_solvated_neutral_system(stage2_prepared_assets: dict, tmp_path) -> None:
    """Validate system setup adds solvent, ions, and applies target force fields."""
    protein_path: Path = stage2_prepared_assets["protein_prepped"]
    ligand_path: Path = stage2_prepared_assets["ligand_prepped"]

    # Load protein
    protein_pdb = app.PDBFile(str(protein_path))

    # Convert prepared ligand SDF to PDB for merging
    from rdkit import Chem

    ligand_mol = Chem.SDMolSupplier(str(ligand_path), removeHs=False)[0]
    assert ligand_mol is not None, "Ligand could not be loaded from prepared SDF"

    ligand_pdb_path = tmp_path / f"{ligand_path.stem}.pdb"
    Chem.MolToPDBFile(ligand_mol, str(ligand_pdb_path))
    ligand_pdb = app.PDBFile(str(ligand_pdb_path))

    modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
    modeller.add(ligand_pdb.topology, ligand_pdb.positions)
    initial_atoms = modeller.topology.getNumAtoms()

    # Load force fields (protein, solvent, lipids)
    forcefield_files = [
        "amber14/protein.ff14SB.xml",
        "amber14/tip3p_standard.xml",
        "amber14/lipid17.xml",
        str(Path("forcefields/shahrokh_heme_ic6_unl.ffxml")),
    ]
    try:
        forcefield = app.ForceField(*forcefield_files)
    except Exception as exc:  # pragma: no cover - environment specific availability
        pytest.skip(f"Required force field XML not available: {exc}")

    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
        neutralize=True,
    )

    total_atoms = modeller.topology.getNumAtoms()
    assert total_atoms > initial_atoms, "Expected solvent to increase atom count"

    water_resnames = {res.name for res in modeller.topology.residues() if res.name in {"HOH", "WAT", "SOL"}}
    assert water_resnames, "Expected water residues after solvation"

    ion_resnames = {res.name for res in modeller.topology.residues() if res.name.upper() in {"NA", "CL", "K", "CA"}}
    assert ion_resnames, "Expected neutralizing ions after system setup"

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        constraints=app.HBonds,
    )
    assert system.getNumParticles() == total_atoms
    assert system.getNumForces() > 0

    output_dir = Path("data/work/system_setup")
    output_dir.mkdir(parents=True, exist_ok=True)

    solvated_pdb = output_dir / f"{protein_path.stem}_stage3_solvated.pdb"
    with open(solvated_pdb, "w") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)

    system_xml = output_dir / f"{protein_path.stem}_stage3_system.xml"
    system_xml.write_text(XmlSerializer.serialize(system))

    assert solvated_pdb.exists()
    assert system_xml.exists()
    assert solvated_pdb.stat().st_size > 0
    assert system_xml.stat().st_size > 0
