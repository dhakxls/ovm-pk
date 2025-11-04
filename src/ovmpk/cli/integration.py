"""Integration with existing OVM-PK components.

This module provides functions to integrate the CLI with the existing
codebase, reusing components like fetchers and the pipeline.
"""

from pathlib import Path
from typing import Dict, Any, Optional, Tuple
import json

from ovmpk.fetchers.protein_fetcher import fetch as fetch_protein
from ovmpk.fetchers.ligand_fetcher import fetch as fetch_ligand
from ovmpk import Pipeline
from ovmpk.config import load_config, update_config

def fetch_protein_structure(uniprot_id: str, output_dir: Path, config: Optional[Dict[str, Any]] = None) -> Path:
    """Fetch a protein structure using the existing protein fetcher.
    
    Args:
        uniprot_id: UniProt ID of the protein
        output_dir: Directory to save the structure
        config: Optional configuration overrides
        
    Returns:
        Path to the downloaded structure file
    """
    # Start with default config
    cfg = load_config()
    
    # Apply overrides
    if config:
        cfg = update_config(cfg, config)
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Use the existing protein fetcher
    result = fetch_protein(uniprot_id, cfg)
    return result["best_match"]

def fetch_ligand_structure(identifier: str, identifier_type: str, output_dir: Path, 
                          config: Optional[Dict[str, Any]] = None) -> Path:
    """Fetch a ligand structure using the existing ligand fetcher.
    
    Args:
        identifier: The ligand identifier (name, SMILES, etc.)
        identifier_type: Type of identifier ('name', 'smiles', 'cid', 'inchikey')
        output_dir: Directory to save the structure
        config: Optional configuration overrides
        
    Returns:
        Path to the downloaded ligand file
    """
    # Start with default config
    cfg = load_config()
    
    # Apply overrides and set identifier type
    if config is None:
        config = {}
    config.setdefault("fetch", {}).setdefault("ligand", {})["identifier_type"] = identifier_type
    cfg = update_config(cfg, config)
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Use the existing ligand fetcher
    return fetch_ligand(identifier, cfg)

def run_analysis(config: Dict[str, Any]) -> Dict[str, Any]:
    """Run the analysis pipeline with the given configuration.
    
    Args:
        config: Configuration dictionary for the pipeline
        
    Returns:
        Results from the pipeline
    """
    # Create output directory if it doesn't exist
    output_dir = Path(config.get('output_dir', 'ovmpk_results'))
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Save the config to the output directory
    config_path = output_dir / 'config.yml'
    with open(config_path, 'w') as f:
        import yaml
        yaml.safe_dump(config, f)
    
    # Create and run the pipeline with the saved config
    pipeline = Pipeline(config_path=str(config_path))
    return pipeline.run()

def create_config(pdb_path: Path, ligand_path: Optional[Path] = None, 
                 output_dir: Optional[Path] = None) -> Dict[str, Any]:
    """Create a configuration dictionary for the pipeline.
    
    Args:
        pdb_path: Path to the protein structure file
        ligand_path: Optional path to the ligand file
        output_dir: Optional output directory for results
        
    Returns:
        Configuration dictionary
    """
    # Start with default config
    config = {
        "system": {
            "protein_pdb": str(pdb_path),  # Required by pipeline
            "protein": {
                "pdb_file": str(pdb_path),
                "prepare": {
                    "add_hydrogens": True,
                    "add_solvent": True,
                    "solvent_padding": 1.0,
                    "ionic_strength": 0.15,
                    "neutralize": True
                }
            },
            "output_dir": str(output_dir) if output_dir else "output"
        },
        "docking": {
            "method": "smina",
            "exhaustiveness": 8,
            "num_poses": 20
        },
        "scoring": {
            "method": "physics",
            "parameters": {
                "metal_centers": [],
                "solvation": "gbvi",
                "dielectric": 78.5
            }
        },
        "simulation": {
            "forcefield": ["amber14-all.xml", "amber14/tip3pfb.xml"],
            "minimization": {"tolerance": 10.0, "max_iterations": 1000},
            "equilibration": {"nvt_steps": 1000, "npt_steps": 1000, "temperature": 300.0},
            "production": {"steps": 10000, "report_interval": 1000, "temperature": 300.0}
        }
    }
    
    # Add ligand if provided
    if ligand_path:
        ligand_pdb = ligand_path.with_suffix('.pdb')
        config["system"]["ligand_pdb"] = str(ligand_pdb)  # Required by pipeline
        config["system"]["ligand"] = {
            "sdf_file": str(ligand_path),
            "pdb_file": str(ligand_pdb),
            "protonation": "epik",
            "tautomers": True
        }
        
        # Create a PDB version of the ligand if it doesn't exist
        if not ligand_pdb.exists():
            from openbabel import pybel
            try:
                mol = next(pybel.readfile('sdf', str(ligand_path)))
                mol.write('pdb', str(ligand_pdb), overwrite=True)
            except Exception as e:
                print(f"Warning: Could not convert ligand to PDB: {e}")
    
    return config
