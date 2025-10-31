"""Benchmark comparison report generator."""
import matplotlib
matplotlib.use('Agg')  # Set before importing pyplot
import matplotlib.pyplot as plt
import pandas as pd
from ovmpk import Pipeline
from ovmpk.utils.reporting import load_benchmarks, dg_to_ki
import math
import os
from pathlib import Path

# Ensure directory exists
Path("reports").mkdir(exist_ok=True)

def generate_report():
    """Generate comparison report with plots and stats"""
    # Load benchmark systems
    systems = [
        {
            "protein": {"pdb": "5VCC"},
            "ligand": {"name": "ketoconazole"},
            "experimental": {"ki_nm": 26.7, "dg_kcal": -10.3},
            "docking": {
                "box": {
                    "auto_center": "heme_fe",
                    "fallback_center": [-22.578, -26.236, -11.422],
                    "size": [12, 12, 12]
                }
            }
        }
    ]
    
    results = []
    for system in systems:
        # Run pipeline with complete config
        pred = Pipeline({
            'protein': system['protein'],
            'ligand': system['ligand'],
            'docking': system['docking'],
            'physics': {'module': 'metal'}
        }).run()
        results.append({
            "system": f"{system['protein']['pdb']}-{system['ligand']['name']}",
            "exp_dg": system["experimental"]["dg_kcal"],
            "pred_dg": pred["physics_scores"]["total"],
            "exp_ki": system["experimental"]["ki_nm"], 
            "pred_ki": dg_to_ki(pred["physics_scores"]["total"])
        })
    
    # Generate report
    df = pd.DataFrame(results)
    plot_correlations(df)
    save_stats(df, "reports")
    return df

def plot_correlations(df):
    """Create comparison plots"""
    plt.figure(figsize=(12,5))
    
    # ΔG plot
    plt.subplot(121)
    plt.scatter(df["exp_dg"], df["pred_dg"], c='b')
    plt.plot([-15,-5], [-15,-5], 'r--')
    plt.xlabel("Experimental ΔG (kcal/mol)")
    plt.ylabel("Predicted ΔG")
    
    # Ki plot 
    plt.subplot(122) 
    plt.loglog(df["exp_ki"], df["pred_ki"], 'bo')
    plt.loglog([1e-3,1e6], [1e-3,1e6], 'r--')
    plt.xlabel("Experimental Ki (nM)")
    plt.ylabel("Predicted Ki")
    
    plt.tight_layout()
    plt.savefig("reports/benchmark_correlation.png")

def save_stats(df, output_dir):
    pass

if __name__ == "__main__":
    generate_report()
