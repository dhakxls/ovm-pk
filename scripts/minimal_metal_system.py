"""Minimal working metal parameterization."""
from openmm import app, System, Vec3
from openmm import unit
import openmm as mm

# Create minimal topology
topology = app.Topology()
chain = topology.addChain()
residue = topology.addResidue("METAL", chain)

# Add atoms (Fe + 4 N)
fe = topology.addAtom("FE", app.Element.getBySymbol('Fe'), residue)
n1 = topology.addAtom("N1", app.Element.getBySymbol('N'), residue)
n2 = topology.addAtom("N2", app.Element.getBySymbol('N'), residue) 
n3 = topology.addAtom("N3", app.Element.getBySymbol('N'), residue)
n4 = topology.addAtom("N4", app.Element.getBySymbol('N'), residue)

# Create positions (square planar)
positions = [
    Vec3(0, 0, 0),  # Fe
    Vec3(2, 0, 0),  # N1
    Vec3(0, 2, 0),  # N2
    Vec3(-2, 0, 0), # N3
    Vec3(0, -2, 0)  # N4
]

# Create empty System
system = mm.System()

# Add particles
system.addParticle(55.845*unit.amu)  # Fe
for _ in range(4):
    system.addParticle(14.007*unit.amu)  # N

# Add bonds
bond_force = mm.HarmonicBondForce()
bond_force.addBond(0, 1, 2.1*unit.angstroms, 100*unit.kilojoules_per_mole/unit.nanometer**2)
bond_force.addBond(0, 2, 2.1*unit.angstroms, 100*unit.kilojoules_per_mole/unit.nanometer**2)
bond_force.addBond(0, 3, 2.1*unit.angstroms, 100*unit.kilojoules_per_mole/unit.nanometer**2)
bond_force.addBond(0, 4, 2.1*unit.angstroms, 100*unit.kilojoules_per_mole/unit.nanometer**2)
system.addForce(bond_force)

# Add nonbonded interactions
nonbonded = mm.NonbondedForce()
nonbonded.addParticle(2.0*unit.elementary_charge, 0.25*unit.nanometer, 0.5*unit.kilojoule_per_mole)  # Fe
for _ in range(4):
    nonbonded.addParticle(-0.5*unit.elementary_charge, 0.31*unit.nanometer, 0.71*unit.kilojoule_per_mole)  # N
system.addForce(nonbonded)

print("Successfully created minimal metal coordination system")
print(f"System has {system.getNumParticles()} particles")
print(f"System has {system.getNumForces()} forces")
