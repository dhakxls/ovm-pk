"""Tests for automated ligand/enzyme curation workflow scaffolding."""

from pathlib import Path
from typing import Dict

import pytest

from src.ovmpk.curation.input_manager import AutoCurateConfig, load_config
from src.ovmpk.curation.workflow import AutoCurator


def test_load_config_from_mapping(tmp_path: Path) -> None:
    config_data: Dict[str, str] = {
        "protein_input": "5VCC",
        "ligand_input": "ketoconazole",
        "run_root": str(tmp_path),
        "run_name": "run_test",
    }

    config = load_config(config_data)
    assert isinstance(config, AutoCurateConfig)
    assert config.run_dir == Path(tmp_path / "run_test")
    assert config.protein_input == "5VCC"
    assert config.ligand_input == "ketoconazole"


def test_autocurator_prepares_directories(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    config = load_config(
        {
            "protein_input": "5VCC",
            "ligand_input": "ketoconazole",
            "run_root": str(tmp_path),
            "run_name": "run_test",
        }
    )

    class DummyFetcher:
        def __init__(self, _cfg):
            pass

        def fetch_protein(self, protein_input: str) -> Path:
            assert protein_input == "5VCC"
            path = tmp_path / "protein.pdb"
            path.write_text("ATOM")
            return path

        def fetch_ligand(self, ligand_input: str) -> Path:
            assert ligand_input == "ketoconazole"
            path = tmp_path / "ligand.sdf"
            path.write_text("$$$$")
            return path

    monkeypatch.setattr("src.ovmpk.curation.workflow.FlexibleFetcher", DummyFetcher)
    monkeypatch.setattr("src.ovmpk.curation.workflow.protein_prep.prepare", lambda paths, cfg: tmp_path / "protein_prepped.pdb")
    monkeypatch.setattr("src.ovmpk.curation.workflow.ligand_prep.prepare", lambda path, cfg: tmp_path / "ligand_prepped.sdf")
    monkeypatch.setattr("src.ovmpk.curation.workflow.AutoCurator._setup_system", lambda self, prot, lig: (None, None))
    monkeypatch.setattr("src.ovmpk.curation.workflow.apply_metal_parameterization", lambda protein, outdir, profile=None: None)

    curator = AutoCurator(config)
    result = curator.run()

    assert result.run_dir.exists()
    assert (result.run_dir / "protein_prep").exists()
    assert (result.run_dir / "ligand_prep").exists()
    assert result.ffxml_bundle is not None
    assert result.ffxml_bundle.ffxml.exists()
