"""Tests for physics module interface and energy calculations."""
from __future__ import annotations

import pytest
from pathlib import Path
from typing import Protocol

# Test ordering - runs after benchmarks but before comparative tests
pytestmark = [
    pytest.mark.order(21),
    pytest.mark.dependency(depends=["test_benchmark_affinity.py"])
]

class PhysicsModule(Protocol):
    """Required interface for physics modules."""
    def score_pose(self, protein_pdb: str, ligand_sdf: str) -> dict[str, float]:
        """Return energy components and total score."""
        ...

def test_physics_registration():
    """Test that custom physics modules can be registered."""
    from ovmpk.docking import physics_registry
    
    class TestPhysics:
        def score_pose(self, protein_pdb, ligand_sdf):
            return {"total": -10.5, "metal": -2.1}
    
    physics_registry.register("test_physics", TestPhysics())
    assert "test_physics" in physics_registry.list_modules()

def test_energy_components():
    """Test physics modules return required energy components."""
    from ovmpk.docking import get_physics_module
    
    physics = get_physics_module("default")
    result = physics.score_pose(
        "data/work/prep/5VCC_fixed.pdb",
        "data/work/docking/poses.sdf"
    )
    
    assert "total" in result
    assert isinstance(result["total"], float)

# Example metal coordination test (will expand with real structures)
def test_metal_coordination():
    """Validate metal-ligand interaction energy calculations."""
    from ovmpk.docking.metal_physics import MetalPhysics
    
    physics = MetalPhysics()
    test_pdb = "tests/data/heme_complex.pdb"
    test_sdf = "tests/data/ligand_coord.sdf"
    
    result = physics.score_pose(test_pdb, test_sdf)
    assert result["metal"] < 0  # Should be stabilizing
    assert "fe_dist" in result  # Should include geometry metrics
