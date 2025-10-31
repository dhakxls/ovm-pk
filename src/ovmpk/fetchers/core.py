"""Generalized protein/ligand fetcher with configurable sources."""
from __future__ import annotations
from pathlib import Path
from typing import Optional, Dict, Any
import requests
from .pdb_resolver import resolve_pdb

class FlexibleFetcher:
    """Handles multiple retrieval methods for proteins/ligands."""
    
    def __init__(self, config: Dict[str, Any]):
        self.config = config
        # Force use of configured directories
        self.protein_dir = Path(config['protein']['save_dir'])
        self.ligand_dir = Path(config['ligand']['save_dir'])
        self.protein_dir.mkdir(parents=True, exist_ok=True)
        self.ligand_dir.mkdir(parents=True, exist_ok=True)
    
    def fetch_protein(self, protein_input: str) -> Path:
        """Smart protein fetcher with auto-resolution"""
        resolver_cfg = self.config.get('protein', {}).get('resolver', {})
        pdb_id = resolve_pdb(protein_input, resolver_cfg)
        
        # Download from RCSB
        return self._download_rcsb(pdb_id, self.protein_dir/f"{pdb_id}.pdb")
    
    def fetch_ligand(self, name: str) -> Path:
        """Get ligand from configured source."""
        id_type = self.config["ligand"].get("identifier_type", "name")
        
        if id_type == "name":
            return self._fetch_pubchem_ligand(
                name,
                self.ligand_dir / f"{name}.sdf"
            )
        elif id_type == "smiles":
            return self._generate_from_smiles(
                self.config["ligand"]["smiles"],
                self.ligand_dir / "ligand.sdf"
            )
        # Add other ID types
    
    # Example implementations
    def _fetch_rcsb_pdb(self, pdb_id: str, out_path: Path) -> Path:
        """Fetch from RCSB PDB."""
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        out_path.write_text(requests.get(url).text)
        return out_path
    
    def _fetch_pubchem_ligand(self, name: str, out_path: Path) -> Path:
        """Fetch ligand by name from PubChem."""
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/SDF"
        out_path.write_bytes(requests.get(url).content)
        return out_path
    
    def _generate_from_smiles(self, smiles: str, out_path: Path) -> Path:
        """Generate 3D structure from SMILES."""
        # Would use RDKit or OpenBabel here
        return out_path
    
    def _download_rcsb(self, pdb_id: str, out_path: Path) -> Path:
        """Fetch from RCSB PDB."""
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        out_path.write_text(requests.get(url).text)
        return out_path

def looks_like_pdb_id(protein_input: str) -> bool:
    # implement logic to check if protein_input looks like a PDB ID
    pass
