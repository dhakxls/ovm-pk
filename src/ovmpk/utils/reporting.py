"""Benchmark reporting utilities."""
import yaml
import math
import json
from pathlib import Path
from typing import Dict
import numpy as np

def dg_to_ki(dg_kcal, T=298.15):
    """Convert ΔG to Ki in nM"""
    R = 0.0019872041  # kcal/(mol·K)
    return math.exp(dg_kcal / (R * T)) * 1e9

def load_benchmarks(config_path):
    """Load benchmark systems from YAML"""
    with open(config_path) as f:
        return yaml.safe_load(f)["systems"]

def generate_report(results: Dict, config: Dict) -> str:
    """Create report with temperature-aware metrics."""
    temp_k = config['conditions']['temperature_K']
    base_report = {
        "deltaG_kcal": results['physics_scores']['total'],
        "Ki_nM": results['physics_scores']['ki_nm'],
        "temperature_K": temp_k,
        "conditions": {
            "pH": config['conditions']['pH'],
            "ionic_strength_M": config['conditions']['ionic_strength_M']
        }
    }
    report = {
        **base_report,
        "analysis": {
            "min_overlap": float(np.min(results['overlaps']['matrix'])),
            "converged": results['converged'],
            "neff": int(np.min(results['overlaps']['neff']))
        }
    }
    return json.dumps(report, indent=2)
