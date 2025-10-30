"""Flexible pipeline core."""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Any
from .fetchers import FlexibleFetcher
from .docking import run as run_docking
from ovmpk.physics import physics_registry

class Pipeline:
    """Configurable pipeline supporting arbitrary protein-ligand pairs."""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        self.fetcher = FlexibleFetcher(config)
    
    def run(self) -> Dict[str, Any]:
        """Execute full pipeline for configured system."""
        # 1. Fetch inputs
        protein = self.fetcher.fetch_protein(Path("data/work/protein"))
        ligand = self.fetcher.fetch_ligand(Path("data/work/ligand"))
        
        # 2. Run docking
        docking_outputs = run_docking(
            protein_pdb_in=str(protein),
            ligand_sdf_in=str(ligand),
            cfg=self.config
        )
        
        # 3. Get physics scores
        physics = physics_registry.get(self.config.get("physics_module", "smina"))
        scores = physics.score_pose(str(protein), docking_outputs[0])
        
        return {
            "docking_results": docking_outputs,
            "physics_scores": scores
        }
    
    def _calculate_dg_error(self, results: Dict[str, Any]) -> float:
        """Compare predicted vs experimental ΔG."""
        pred = results["physics_scores"]["total"]
        exp = self.config["experimental"]["dg_kcal"]
        return abs(pred - exp)

# Example usage:
# config = {"protein": {"pdb": "5VCC"}, "ligand": {"name": "ketoconazole"}}
# pipeline = Pipeline(config)
# results = pipeline.run()
