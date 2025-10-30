"""Physics module implementations."""
from .registry import physics_registry
from .scorers import SminaScorer, MetalPhysics

# Register default modules
physics_registry.register("smina", SminaScorer())
physics_registry.register("metal", MetalPhysics())

__all__ = ['physics_registry']
