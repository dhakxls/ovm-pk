"""Operational robustness tests."""
import pytest
from ovmpk import Pipeline

@pytest.mark.order(40)
def test_error_handling():
    """Test pipeline error handling."""
    with pytest.raises(RuntimeError):
        Pipeline({
            "protein": {"pdb": "INVALID"},
            "ligand": {"name": "ketoconazole"}
        }).run()
