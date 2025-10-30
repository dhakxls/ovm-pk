"""Physics module registry system."""
from typing import Dict, Any
from pathlib import Path

class PhysicsRegistry:
    """Registry for pluggable physics modules."""
    def __init__(self):
        self._modules = {}
    
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

physics_registry = PhysicsRegistry()
