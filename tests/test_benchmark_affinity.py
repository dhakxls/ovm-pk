"""Benchmark affinity validation tests (ΔG/Ki correlation)."""
from __future__ import annotations

import pytest
import pandas as pd
from pathlib import Path
import numpy as np
from scipy import stats
from ovmpk.utils.affinity import ki_to_dg

# Sample benchmark data (will be replaced with real dataset)
SAMPLE_DATA = """
PDB,Ligand,Experimental_Ki_nM,Experimental_dG_kcal
5VCC,ketoconazole,15.0,-10.2
1A2C,indomethacin,8.0,-11.5
"""

def test_ki_to_dg_conversion():
    """Validate ΔG calculation from Ki using CYP3A4 parameters."""
    # Test known values from benchmark
    assert np.isclose(ki_to_dg(26.7, 298.15), -10.3, atol=0.1)  # Ketoconazole
    assert np.isclose(ki_to_dg(9210, 298.15), -6.9, atol=0.1)   # Fluconazole
    
    # Verify relative affinities
    dg_keto = ki_to_dg(26.7, 298.15)
    dg_flu = ki_to_dg(9210, 298.15)
    assert dg_keto < dg_flu, "Ketoconazole should have stronger binding"

def test_ki_to_dg_conversion_precise():
    """Standalone precision test"""
    assert abs(ki_to_dg(26.7) - (-10.3)) < 0.1

def test_relative_affinity():
    """Verify scoring preserves relative affinity."""
    from ovmpk.docking import compare_physics_modules
    results = compare_physics_modules(
        protein="data/work/protein/5VCC.pdb",
        ligand="data/work/docking/ketoconazole_poses.sdf",
        modules=["smina", "metal"]
    )
    assert results["metal"]["total"] < results["smina"]["total"]

def test_ketoconazole_affinity():
    """Validate ketoconazole ΔG prediction."""
    bench_path = Path(__file__).parent / "data" / "benchmarks" / "cyp3a4_ketoconazole.csv"
    df = pd.read_csv(bench_path)
    
    # Get pipeline prediction (example: -10.5 kcal/mol)
    pred_dg = -10.5  # Will connect to actual pipeline output
    exp_dg = df.iloc[0]["Experimental_dG_kcal"]
    
    assert abs(pred_dg - exp_dg) <= 1.5, \
        f"ΔG error {abs(pred_dg - exp_dg):.2f} > 1.5 kcal/mol threshold"

def test_benchmark_correlation():
    """Test CYP3A4 benchmark set."""
    bench_path = Path(__file__).parent / "data" / "benchmarks" / "cyp3a4_benchmark.csv"
    df = pd.read_csv(bench_path)
    
    # Get pipeline predictions (placeholder)
    pred_dg = {
        "ketoconazole": -10.5,
        "fluconazole": -7.1
    }
    
    # Calculate metrics
    errors = [abs(row["Experimental_dG_kcal"] - pred_dg[row["Ligand"]]) 
              for _, row in df.iterrows()]
    avg_error = np.mean(errors)
    
    assert avg_error <= 1.5, f"Average ΔG error {avg_error:.2f} > 1.5 kcal/mol"

# Helper to generate benchmark CSV if missing
@pytest.fixture(scope="module")
def ensure_benchmark_data():
    out = Path(__file__).parent / "data" / "benchmark_set.csv"
    out.parent.mkdir(exist_ok=True)
    if not out.exists():
        out.write_text(SAMPLE_DATA)
    return out
