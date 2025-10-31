"""Stage 1 fetching tests that pull defaults from configuration env file."""

from datetime import datetime
from pathlib import Path

import pytest

from src.ovmpk.config.env_defaults import resolve_default_inputs
from src.ovmpk.fetchers.core import FlexibleFetcher
from src.ovmpk.fetchers.pdb_resolver import resolve_pdb


@pytest.fixture(scope="module")
def test_config(request):
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

    if request.config.getoption("--protein"):
        resolver_preferences.update({
            "target_species_id": 9606,
            "exp_method": "X-RAY DIFFRACTION",
            "max_resolution": 2.5,
            "results_limit": 25,
        })

    if not protein_token or not ligand_name:
        pytest.exit("protein_ligand_env.md must define both 'protein' and 'ligand'.")

    pdb_id = resolve_pdb(protein_token, resolver_preferences)

    # Create run directory
    run_dir = Path(f"test_run/run_{datetime.now().strftime('%Y%m%d_%H%M%S')}")
    run_dir.mkdir(parents=True, exist_ok=True)

    return {
        'protein_token': protein_token,
        'pdb_id': pdb_id,
        'ligand_name': ligand_name,
        'fetcher': FlexibleFetcher({
            'protein': {
                'save_dir': str(run_dir),
                'resolver': resolver_preferences,
            },
            'ligand': {
                'save_dir': str(run_dir),
                'identifier_type': 'name',
            }
        }),
        'run_dir': run_dir
    }

def test_fetch_protein(test_config):
    """Test protein download"""
    protein_path = test_config['fetcher'].fetch_protein(test_config['pdb_id'])
    assert protein_path.exists()
    assert protein_path.suffix == ".pdb"
    assert "HEADER" in protein_path.read_text()[:100]

@pytest.mark.slow
def test_fetch_ligand(test_config):
    """Test ligand download - validates PubChem SDF format"""
    ligand_path = test_config['fetcher'].fetch_ligand(test_config['ligand_name'])
    
    # Basic file checks
    assert ligand_path.exists()
    assert ligand_path.suffix == ".sdf"
    assert ligand_path.stat().st_size > 0
    
    # PubChem-specific SDF validation
    content = ligand_path.read_text()
    lines = content.splitlines()
    
    # Should contain both V2000 marker and $$$$ separator
    assert any("V2000" in line for line in lines)  # Structure block
    assert "$$$$" in content[-100:]  # File separator

def pytest_addoption(parser):
    """Register custom command-line options using environment defaults."""
    defaults = ensure_defaults(
        load_protein_ligand_inputs(Path("protein_ligand_env.md")),
        {"protein": "", "ligand": ""},
    )

    parser.addoption(
        "--protein",
        action="store",
        default=defaults.get("protein", ""),
        help="Protein input (defaults to value in protein_ligand_env.md)",
    )
    parser.addoption(
        "--ligand",
        action="store",
        default=defaults.get("ligand", ""),
        help="Ligand name (defaults to value in protein_ligand_env.md)",
    )
