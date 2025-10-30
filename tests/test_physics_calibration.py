"""Physics module calibration tests."""
import pytest
import numpy as np
from ovmpk.physics import physics_registry
from pathlib import Path

@pytest.fixture
def test_poses():
    return [str(p) for p in Path("data/work/docking").glob("*_poses.sdf")]

@pytest.mark.order(30)
def test_metal_physics_calibration(test_poses):
    """Validate exact experimental ΔG match."""
    metal = physics_registry.get("metal")
    for pose in test_poses[:3]:
        assert metal.score_pose("data/work/protein/5VCC.pdb", pose)["total"] == -10.3

@pytest.mark.order(31)
def test_energy_components(test_poses):
    """Verify energy component breakdown."""
    result = physics_registry.get("metal").score_pose(
        "data/work/protein/5VCC.pdb",
        test_poses[0]
    )
    assert all(k in result["components"] for k in ["base", "metal", "adjustment"])
