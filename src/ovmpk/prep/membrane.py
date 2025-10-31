"""Membrane preparation using OpenMM (Amber-free)."""
import logging
from pathlib import Path
from openmm.app import Modeller, PDBFile
from openmmforcefields.generators import SystemGenerator

def prepare_membrane(pdb_path: Path, out_dir: Path, config: dict):
    """Prepare membrane system using OpenMM force fields."""
    try:
        # Load protein
        pdb = PDBFile(str(pdb_path))
        
        # Create system generator
        generator = SystemGenerator(
            forcefields=[config['physics']['forcefields']['protein'],
                        config['physics']['forcefields']['lipids']],
            small_molecule_forcefield=config['physics']['forcefields']['ligand'],
            molecules=[],
            cache=None
        )
        
        # Build membrane
        modeller = Modeller(pdb.topology, pdb.positions)
        modeller.addMembrane(
            lipidType='POPC',
            ionicStrength=config['conditions']['ionic_strength_M'],
            positiveIon='Na+',
            negativeIon='Cl-'
        )
        
        # Save output
        output_path = out_dir / "membrane.pdb"
        with open(output_path, 'w') as f:
            PDBFile.writeFile(modeller.topology, modeller.positions, f)
            
        return output_path
        
    except Exception as e:
        logging.error(f"Membrane preparation failed: {str(e)}")
        raise
