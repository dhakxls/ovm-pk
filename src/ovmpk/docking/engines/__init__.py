"""
Docking Engines Interface
Supported engines: GNINA, VINA, Custom
"""
from .base import DockingEngine
from .smina_engine import SminaEngine

__all__ = ['DockingEngine', 'SminaEngine']