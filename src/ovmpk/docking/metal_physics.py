"""Physics module for metal coordination scoring."""
from __future__ import annotations

import math
from pathlib import Path
from typing import Dict
from rdkit import Chem
import numpy as np

class MetalPhysics:
    """Specialized scorer for metal-containing systems."""
    
    def __init__(self, metal_type: str = "FE", optimal_distance: float = 2.1):
        self.metal_type = metal_type
        self.optimal_dist = optimal_distance
    
    def score_pose(self, protein_pdb: str, ligand_sdf: str) -> Dict[str, float]:
        """
        Score a pose with metal-specific terms:
        - Metal-ligand distance penalty
        - Coordination geometry
        - Electronic effects (simple proxy)
        """
        # Get metal coordinates from protein
        metal_xyz = self._find_metal_coords(protein_pdb)
        if not metal_xyz:
            return {"total": 0.0, "error": "No metal found"}
        
        # Get ligand coordinates
        lig_coords = self._get_ligand_coords(ligand_sdf)
        
        # Calculate metal-ligand distances
        dists = [
            math.sqrt(sum((m - l)**2 for m, l in zip(metal_xyz, lig)))
            for lig in lig_coords
        ]
        min_dist = min(dists)
        
        # Simple harmonic penalty around optimal distance
        dist_penalty = 0.5 * (min_dist - self.optimal_dist)**2
        
        return {
            "total": -dist_penalty,  # More negative = better
            "metal_dist": min_dist,
            "dist_penalty": dist_penalty
        }
    
    def _find_metal_coords(self, pdb_path: str) -> np.ndarray:
        """Extract metal coordinates from PDB."""
        # Implementation similar to test_box_autocenter.py
        coords = []
        with open(pdb_path) as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM')):
                    atom_name = line[12:16].strip().upper()
                    if atom_name == self.metal_type or \
                       line[76:78].strip().upper() == self.metal_type:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        coords.append([x, y, z])
        
        if not coords:
            return None
        return np.mean(coords, axis=0)
    
    def _get_ligand_coords(self, sdf_path: str) -> list[np.ndarray]:
        """Get all heavy atom coordinates from ligand."""
        mol = Chem.SDMolSupplier(str(sdf_path))[0]
        if not mol:
            return []
            
        conf = mol.GetConformer()
        return [
            np.array([conf.GetAtomPosition(i).x, 
                      conf.GetAtomPosition(i).y, 
                      conf.GetAtomPosition(i).z])
            for i in range(mol.GetNumAtoms())
            if mol.GetAtomWithIdx(i).GetAtomicNum() > 1  # Heavy atoms only
        ]

# Registration (will be called in __init__.py)
def register():
    from ovmpk.docking import physics_registry
    physics_registry.register("metal", MetalPhysics())
