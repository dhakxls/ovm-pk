"""Shared pytest configuration, fixtures, and helpers."""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path
from typing import Dict, Any

import pytest

from src.ovmpk.config.env_defaults import (
    ensure_defaults,
    load_protein_ligand_inputs,
    resolve_default_inputs,
)
from src.ovmpk.fetchers.core import FlexibleFetcher
from src.ovmpk.fetchers.pdb_resolver import resolve_pdb
from src.ovmpk.prep import ligand_prep, protein_prep
from src.ovmpk.utils.run_dirs import set_run_context, stage_dir


def pytest_addoption(parser):
    defaults = ensure_defaults(
        load_protein_ligand_inputs(Path("protein_ligand_env.md")),
        {"protein": "", "ligand": ""},
    )

    parser.addoption(
        "--protein",
        action="store",
        default=defaults.get("protein", ""),
        help="PDB ID or protein name (defaults to value in protein_ligand_env.md)",
    )
    parser.addoption(
        "--ligand",
        action="store",
        default=defaults.get("ligand", ""),
        help="Ligand name (defaults to value in protein_ligand_env.md)",
    )


@pytest.fixture(scope="module")
def test_config(request) -> Dict[str, Any]:
    """Resolve Stage 1 defaults and instantiate fetchers for downstream tests."""

    inputs = resolve_default_inputs(
        Path("protein_ligand_env.md"),
        {
            "target_species_id": 9606,
            "exp_method": "X-RAY DIFFRACTION",
            "max_resolution": 2.5,
            "results_limit": 25,
        },
    )

    protein_token = request.config.getoption("--protein") or inputs["protein_token"]
    ligand_name = request.config.getoption("--ligand") or inputs["ligand"]

    resolver_preferences = inputs["resolver_preferences"]
    if request.config.getoption("--protein"):
        resolver_preferences = dict(resolver_preferences)
        resolver_preferences.update(
            {
                "target_species_id": 9606,
                "exp_method": "X-RAY DIFFRACTION",
                "max_resolution": 2.5,
                "results_limit": 25,
            }
        )

    if not protein_token or not ligand_name:
        pytest.exit("protein_ligand_env.md must define both 'protein' and 'ligand'.")

    pdb_id = resolve_pdb(protein_token, resolver_preferences)

    run_dir_override = os.environ.get("OVMPK_STAGE1_RUN_DIR")
    if run_dir_override:
        run_dir = Path(run_dir_override)
    else:
        run_dir = Path(f"test_run/run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    run_dir.mkdir(parents=True, exist_ok=True)
    set_run_context(run_dir.parent, run_dir.name)

    stage1_dir = stage_dir("stage1")
    protein_stage1_dir = stage_dir("stage1/protein")
    ligand_stage1_dir = stage_dir("stage1/ligand")

    fetcher = FlexibleFetcher(
        {
            "protein": {
                "save_dir": str(protein_stage1_dir),
                "resolver": resolver_preferences,
            },
            "ligand": {
                "save_dir": str(ligand_stage1_dir),
                "identifier_type": "name",
            },
        }
    )

    return {
        "protein_token": protein_token,
        "pdb_id": pdb_id,
        "ligand_name": ligand_name,
        "fetcher": fetcher,
        "run_dir": run_dir,
        "stage1_dir": stage1_dir,
        "resolver_preferences": resolver_preferences,
    }


@pytest.fixture(scope="module")
def stage1_assets(test_config: Dict[str, Any]) -> Dict[str, Any]:
    """Ensure Stage 1 outputs are available for downstream stages."""

    stage1_dir: Path = test_config["stage1_dir"]
    protein_path = stage1_dir / "protein" / f"{test_config['pdb_id']}.pdb"
    ligand_path = stage1_dir / "ligand" / f"{test_config['ligand_name']}.sdf"

    if not protein_path.exists():
        protein_path = Path(test_config["fetcher"].fetch_protein(test_config["pdb_id"]))
    if not ligand_path.exists():
        ligand_path = Path(test_config["fetcher"].fetch_ligand(test_config["ligand_name"]))

    return {
        **test_config,
        "protein_path": protein_path,
        "ligand_path": ligand_path,
    }


@pytest.fixture(scope="module")
def stage2_prepared_assets(stage1_assets: Dict[str, Any]) -> Dict[str, Any]:
    """Run Stage 2 prep once for reuse in Stage 3/4 tests."""

    protein_input = stage1_assets["protein_path"]
    ligand_input = stage1_assets["ligand_path"]
    run_dir: Path = stage1_assets["run_dir"]

    protein_out_dir = stage_dir("protein_prep")
    ligand_out_dir = stage_dir("ligand_prep")

    prot_cfg = {
        "prep": {
            "protein": {
                "ph": 7.4,
                "output_suffix": "_stage2_prepped",
                "run_pdbfixer": True,
                "output_dir": str(protein_out_dir),
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
                "use_dimorphite": False,
                "output_dir": str(ligand_out_dir),
            }
        }
    }

    protein_prepped = protein_prep.prepare({"apo": protein_input}, prot_cfg)
    ligand_prepped = ligand_prep.prepare(ligand_input, ligand_cfg)

    return {
        **stage1_assets,
        "protein_prepped": Path(protein_prepped),
        "ligand_prepped": Path(ligand_prepped),
    }


def run_stage3_system_setup(stage2_prepared_assets: Dict[str, Any], tmp_path: Path) -> Dict[str, Any]:
    """Execute Stage 3 system setup and return metadata about the solvated system."""

    try:
        from openmm import XmlSerializer, app, unit
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise pytest.skip(f"OpenMM not available: {exc}")

    try:
        from rdkit import Chem
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise pytest.skip(f"RDKit not available: {exc}")

    protein_path: Path = stage2_prepared_assets["protein_prepped"]
    ligand_path: Path = stage2_prepared_assets["ligand_prepped"]

    protein_pdb = app.PDBFile(str(protein_path))
    ligand_mol = Chem.SDMolSupplier(str(ligand_path), removeHs=False)[0]
    if ligand_mol is None:
        raise AssertionError("Ligand could not be loaded from prepared SDF")

    ligand_pdb_path = tmp_path / f"{ligand_path.stem}.pdb"
    Chem.MolToPDBFile(ligand_mol, str(ligand_pdb_path))
    ligand_pdb = app.PDBFile(str(ligand_pdb_path))

    modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
    modeller.add(ligand_pdb.topology, ligand_pdb.positions)
    initial_atoms = modeller.topology.getNumAtoms()

    forcefield_files = [
        "amber14/protein.ff14SB.xml",
        "amber14/tip3p_standard.xml",
        "amber14/lipid17.xml",
        str(Path("forcefields/shahrokh_heme_ic6_unl.ffxml")),
    ]
    forcefield = app.ForceField(*forcefield_files)

    modeller.addSolvent(
        forcefield,
        model="tip3p",
        padding=1.0 * unit.nanometer,
        ionicStrength=0.15 * unit.molar,
        neutralize=True,
    )

    total_atoms = modeller.topology.getNumAtoms()

    water_resnames = {
        res.name for res in modeller.topology.residues() if res.name in {"HOH", "WAT", "SOL"}
    }
    ion_resnames = {
        res.name for res in modeller.topology.residues() if res.name.upper() in {"NA", "CL", "K", "CA"}
    }

    system = forcefield.createSystem(
        modeller.topology,
        nonbondedMethod=app.PME,
        constraints=app.HBonds,
    )

    output_dir = stage2_prepared_assets["run_dir"] / "system_setup"
    output_dir.mkdir(parents=True, exist_ok=True)

    solvated_pdb = output_dir / f"{protein_path.stem}_stage3_solvated.pdb"
    with open(solvated_pdb, "w", encoding="utf-8") as fh:
        app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)

    system_xml = output_dir / f"{protein_path.stem}_stage3_system.xml"
    system_xml.write_text(XmlSerializer.serialize(system), encoding="utf-8")

    return {
        "protein_path": protein_path,
        "ligand_path": ligand_path,
        "solvated_pdb": solvated_pdb,
        "system_xml": system_xml,
        "output_dir": output_dir,
        "initial_atoms": initial_atoms,
        "total_atoms": total_atoms,
        "water_resnames": water_resnames,
        "ion_resnames": ion_resnames,
        "system": system,
    }
