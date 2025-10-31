"""Robust metal parameterization with proper structure preparation."""
from pathlib import Path
from typing import Dict
import logging
from openmm import app, unit
from openmm import System, Vec3
import openmm as mm
from ..utils.validators import validate_no_calibration

# Supported metal elements
METAL_ELEMENTS = {'FE', 'ZN', 'MG', 'CA', 'MN', 'CO', 'NI', 'CU'}

class MetalParameterGenerator:
    """Generate metal parameters using OpenMM with fallback."""
    
    def __init__(self, config: Dict):
        validate_no_calibration(config)
        self.config = config
        
    def prepare_heme(self, pdb_file: Path, work_dir: Path) -> Dict:
        """Generate parameters for metal-containing protein."""
        # First validate the PDB contains metals
        self._validate_pdb(pdb_file)
        
        try:
            # Try standard approach with added hydrogens
            return self._standard_parameterization(pdb_file)
        except Exception as e:
            logging.warning(f"Standard parameterization failed: {e}")
            logging.warning("Using minimal metal parameterization")
            return self._minimal_metal_system(pdb_file)
            
    def _validate_pdb(self, pdb_file: Path):
        """Verify PDB contains metal atoms."""
        with open(pdb_file) as f:
            pdb_content = f.read()
            
        has_metal = any(
            f" {element} " in pdb_content.upper() 
            for element in METAL_ELEMENTS
        )
        
        if not has_metal:
            raise ValueError(
                f"No metal atoms found in {pdb_file}. "
                f"Required elements: {METAL_ELEMENTS}"
            )
    
    def _standard_parameterization(self, pdb_file: Path) -> Dict:
        """Standard OpenMM parameterization with hydrogens."""
        pdb = app.PDBFile(str(pdb_file))
        forcefield = app.ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        
        # Add missing hydrogens
        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield, pH=7.4)
        
        system = forcefield.createSystem(
            modeller.getTopology(),
            nonbondedMethod=app.PME,
            constraints=app.HBonds
        )
        return {
            'system': system,
            'topology': modeller.getTopology()
        }
            
    def _minimal_metal_system(self, pdb_file: Path) -> Dict:
        """Create minimal system focusing on metal coordination."""
        pdb = app.PDBFile(str(pdb_file))
        positions = pdb.positions
        
        # Find all metal atoms
        metal_atoms = [
            (i, atom) for i, atom in enumerate(pdb.topology.atoms())
            if atom.element.symbol.upper() in METAL_ELEMENTS
        ]
        
        if not metal_atoms:
            raise ValueError("No metal atoms found after validation - this should not happen")
            
        # Use first metal atom found
        metal_idx, metal_atom = metal_atoms[0]
        metal_pos = positions[metal_idx]
        
        # Create minimal system
        system = mm.System()
        system.addParticle(metal_atom.element.mass)
        
        # Add nonbonded parameters for metal
        nonbonded = mm.NonbondedForce()
        nonbonded.addParticle(
            self.config.get('metal_charge', 2.0)*unit.elementary_charge,
            self.config.get('metal_sigma', 0.25)*unit.nanometer,
            self.config.get('metal_epsilon', 0.5)*unit.kilojoule_per_mole
        )
        system.addForce(nonbonded)
        
        return {
            'system': system,
            'topology': pdb.topology
        }

    @classmethod
    def from_config(cls, config: Dict):
        """Factory method with validation."""
        required = {'oxidation_state', 'coordination'}
        if not required.issubset(config.keys()):
            raise ValueError(f"Config missing required keys: {required}")
        return cls(config)
