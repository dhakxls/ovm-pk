"""Stage 1 fetching tests that reuse fixtures from tests.conftest."""


def test_fetch_protein(test_config):
    """Test protein download"""
    protein_path = test_config['fetcher'].fetch_protein(test_config['pdb_id'])
    assert protein_path.exists()
    assert protein_path.suffix == ".pdb"
    assert "HEADER" in protein_path.read_text()[:100]

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
    assert any("V2000" in line for line in lines)
    assert "$$$$" in content[-100:]
