"""Toolchain utilities (Amber-free version)."""
import logging
from typing import Optional
import subprocess

def require_openmm():
    """Verify OpenMM is available."""
    try:
        import openmm
        return True
    except ImportError:
        logging.error("OpenMM not found - required for physics")
        return False

def check_gpu() -> Optional[str]:
    """Check for available GPU platform."""
    try:
        from openmm import Platform
        for i in range(Platform.getNumPlatforms()):
            if Platform.getPlatform(i).getName() == "CUDA":
                return "CUDA"
        return "CPU"
    except ImportError:
        return None
