"""Experimental validation tests."""
import pytest
from ovmpk import Pipeline

@pytest.mark.order(20)
def test_dg_prediction():
    """Test ΔG prediction against experimental value."""
    config = {
        "protein": {"pdb": "5VCC"},
        "ligand": {"name": "ketoconazole"},
        "experimental": {"dg_kcal": -10.3},
        "docking": {
            "box": {
                "auto_center": "heme_fe",
                "fallback_center": [-15.846, -23.032, -11.293],
                "size": [12.0, 12.0, 12.0]
            }
        }
    }
    results = Pipeline(config).run()
    assert "physics_scores" in results
    assert abs(results["physics_scores"]["total"] - (-10.3)) <= 1.5
