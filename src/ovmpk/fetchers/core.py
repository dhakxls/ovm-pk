"""Generalized protein/ligand fetcher with configurable sources."""
from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict, Any
import requests

class FlexibleFetcher:
    """Handles multiple retrieval methods for proteins/ligands."""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
    
    def fetch_protein(self, out_dir: Path) -> Path:
        """Get protein structure from configured source."""
        out_dir.mkdir(parents=True, exist_ok=True)  # Ensure dir exists
        src = self.config.get("protein", {}).get("source", "rcsb")
        
        if src == "rcsb":
            return self._fetch_rcsb_pdb(
                self.config["protein"]["pdb"], 
                out_dir
            )
        elif src == "custom":
            return Path(self.config["protein"]["local_path"])
        # Add other sources (alphafold, etc.)
        
    def fetch_ligand(self, out_dir: Path) -> Path:
        """Get ligand from configured source."""
        out_dir.mkdir(parents=True, exist_ok=True)  # Ensure dir exists
        id_type = self.config["ligand"].get("identifier_type", "name")
        
        if id_type == "name":
            return self._fetch_by_name(
                self.config["ligand"]["name"],
                out_dir
            )
        elif id_type == "smiles":
            return self._generate_from_smiles(
                self.config["ligand"]["smiles"],
                out_dir
            )
        # Add other ID types
    
    # Example implementations
    def _fetch_rcsb_pdb(self, pdb_id: str, out_dir: Path) -> Path:
        """Fetch from RCSB PDB."""
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        out_path = out_dir / f"{pdb_id}.pdb"
        out_path.write_text(requests.get(url).text)
        return out_path
    
    def _fetch_by_name(self, name: str, out_dir: Path) -> Path:
        """Fetch ligand by name from PubChem."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF"
        out_path = out_dir / f"{name}.sdf"
        out_path.write_bytes(requests.get(url).content)
        return out_path
    
    def _generate_from_smiles(self, smiles: str, out_dir: Path) -> Path:
        """Generate 3D structure from SMILES."""
        out_path = out_dir / "ligand.sdf"
        # Would use RDKit or OpenBabel here
        return out_path
