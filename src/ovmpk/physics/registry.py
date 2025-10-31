"""
Physics component registry
Updated to use new metal.py instead of metal_params
"""
from __future__ import annotations
from typing import Dict, Any
from pathlib import Path
from .metal import parameterize_metal
from .abfe import ABFEEstimator, FreeEnergyCalculator  # Import ABFEEstimator to avoid forward reference
from .restraints import apply_restraints
from .validators import CalibrationForbiddenError

class PhysicsRegistry:
    """Registry for pluggable physics modules."""
    def __init__(self):
        self._modules = {}
        self._metal_handler = None
        self._abfe_estimator = None
    
    def register(self, name: str, module):
        self._modules[name] = module
    
    def get(self, name: str):
        return self._modules.get(name)
    
    def list_modules(self):
        return list(self._modules.keys())

    def register_from_config(self, name: str, config_path: Path):
        """Register a module with parameters from config."""
        import yaml
        cfg = yaml.safe_load(config_path.read_text())
        if name == "metal":
            from .scorers import MetalPhysics
            module = MetalPhysics()
            for param, value in cfg["parameters"].items():
                setattr(module.params, param, value)
            self.register(name, module)

    def register_metal_handler(self, config: Dict[str, Any]):
        """Accepts config under 'physics.metal'."""
        required = {'mode', 'oxidation_state', 'coordination'}
        if not all(k in config for k in required):
            raise ValueError(f"Missing required metal config keys: {required}")
        self._metal_handler = parameterize_metal(config)

    def get_metal_params(self, pdb_file: str, work_dir: str) -> Dict:
        """Get parameters with clear error if uninitialized."""
        if not self._metal_handler:
            raise RuntimeError(
                "Metal handler not initialized. Call register_metal_handler() first"
            )
        return self._metal_handler(pdb_file, work_dir)

    def register_abfe(self, config: Dict[str, Any]):
        """Initialize ABFE estimator."""
        self._abfe_estimator = ABFEEstimator(config)

    def get_abfe(self) -> ABFEEstimator:
        """Retrieve ABFE estimator."""
        if not self._abfe_estimator:
            raise RuntimeError("ABFE estimator not initialized")
        return self._abfe_estimator

# Singleton registry
physics_registry = PhysicsRegistry()
