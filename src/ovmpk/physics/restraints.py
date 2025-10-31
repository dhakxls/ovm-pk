"""Boresch-style restraints for ABFE."""
from typing import Dict, List, Tuple
import numpy as np
from openmm import CustomExternalForce, app

def create_boresch_restraint(
    system, 
    protein_atoms: List[int], 
    ligand_atoms: List[int],
    reference_positions: np.ndarray
) -> Tuple[Dict, CustomExternalForce]:
    """Set up Boresch restraints between protein and ligand.
    
    Returns:
        Tuple of (restraint_params, force_object)
    """
    # Anchor points: protein (r1,r2,r3) and ligand (l1,l2,l3)
    params = {
        "r1": protein_atoms[0], "r2": protein_atoms[1], "r3": protein_atoms[2],
        "l1": ligand_atoms[0], "l2": ligand_atoms[1], "l3": ligand_atoms[2],
        "distance": 0.3,  # nm
        "angle1": 2.0,    # radians
        "angle2": 2.0,
        "k_distance": 100,  # kJ/mol/nm²
        "k_angle": 100      # kJ/mol/rad²
    }
    
    force = CustomExternalForce(
        """
        distance_restraint + angle1_restraint + angle2_restraint;
        distance_restraint = 0.5*k_distance*(distance(r1,l1)-d0)^2;
        angle1_restraint = 0.5*k_angle*(angle(r2,r1,l1)-a1)^2;
        angle2_restraint = 0.5*k_angle*(angle(r1,l1,l2)-a2)^2;
        """
    )
    
    return params, force

"""
Restraints implementation for OpenMM systems
"""
from openmm import app, CustomBondForce

def apply_restraints(system, config=None):
    """Apply harmonic restraints to heavy atoms"""
    if config is None:
        config = {"force_constant": 100}  # kJ/mol/nm^2
        
    restraint_force = CustomBondForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint_force.addPerBondParameter("k")
    restraint_force.addPerBondParameter("x0")
    restraint_force.addPerBondParameter("y0") 
    restraint_force.addPerBondParameter("z0")
    
    system.addForce(restraint_force)
    return system
