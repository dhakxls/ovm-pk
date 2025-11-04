"""
OVM-PK Command Line Interface
----------------------------
This module provides command-line tools for running the OVM-PK pipeline.
Command-line interface for OVM-PK.

This package provides command-line tools for running protein-ligand analysis
using the OVM-PK pipeline.
"""

from .run_analysis import run_interactive_analysis, main

__all__ = ['run_interactive_analysis', 'main']
