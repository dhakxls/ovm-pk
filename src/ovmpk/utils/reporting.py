"""Benchmark reporting utilities."""
import yaml
import math
from pathlib import Path

def dg_to_ki(dg_kcal, T=298.15):
    """Convert ΔG to Ki in nM"""
    R = 0.0019872041  # kcal/(mol·K)
    return math.exp(dg_kcal / (R * T)) * 1e9

def load_benchmarks(config_path):
    """Load benchmark systems from YAML"""
    with open(config_path) as f:
        return yaml.safe_load(f)["systems"]
