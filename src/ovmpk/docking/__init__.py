"""Docking module interface."""
from .smina_wrapper import run, autocenter_box_from_pdb, SminaScorer
from .metal_physics import MetalPhysics

__all__ = ["run", "autocenter_box_from_pdb", "SminaScorer", "MetalPhysics"]

# Physics module registry
class PhysicsRegistry:
    def __init__(self):
        self._modules = {}
    
    def register(self, name: str, module):
        self._modules[name] = module
    
    def get(self, name: str):
        return self._modules.get(name)
    
    def list_modules(self):
        return list(self._modules.keys())

physics_registry = PhysicsRegistry()

# Register default modules
from .smina_wrapper import SminaScorer
physics_registry.register("smina", SminaScorer())

# Add metal physics module
from .metal_physics import MetalPhysics
physics_registry.register("metal", MetalPhysics())