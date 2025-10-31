"""Complete pipeline execution with robust structure preparation."""
import os
import sys
import requests
from pathlib import Path
from openmm.app import PDBFile, Modeller, ForceField
from pdbfixer import PDBFixer
from simtk import unit

# Project paths
PROJECT_ROOT = Path(__file__).parent.parent
INPUT_DIR = PROJECT_ROOT / "data/input/proteins"
WORK_DIR = PROJECT_ROOT / "data/work/protein"

# Ensure directories exist
INPUT_DIR.mkdir(parents=True, exist_ok=True)
WORK_DIR.mkdir(parents=True, exist_ok=True)

# Download PDB if needed
PDB_ID = "5VCC"
input_pdb = INPUT_DIR / f"{PDB_ID}.pdb"
work_pdb = WORK_DIR / f"{PDB_ID}_fixed.pdb"

if not input_pdb.exists():
    print(f"Downloading {PDB_ID}.pdb from RCSB...")
    url = f"https://files.rcsb.org/download/{PDB_ID}.pdb"
    response = requests.get(url)
    response.raise_for_status()
    input_pdb.write_text(response.text)
    print(f"Saved to {input_pdb}")

# Step 1: Use PDBFixer for comprehensive fixes
print("Running PDBFixer...")
fixer = PDBFixer(str(input_pdb))

# Preserve heme group (FE atom)
fixer.removeHeterogens(keepWater=False)
for residue in fixer.topology.residues():
    if any(atom.element.symbol == 'FE' for atom in residue.atoms()):
        fixer.heterogens.append(residue)

# Fix structure
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Save intermediate fixed PDB
intermediate_pdb = WORK_DIR / "5VCC_intermediate.pdb"
PDBFile.writeFile(fixer.topology, fixer.positions, open(intermediate_pdb, 'w'))

# Step 2: Use OpenMM Modeller for hydrogen addition
print("Adding hydrogens with Modeller...")
pdb = PDBFile(str(intermediate_pdb))
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forcefield, pH=7.4)

# Save final prepared PDB
print(f"Saving prepared PDB to {work_pdb}")
PDBFile.writeFile(modeller.topology, modeller.positions, open(work_pdb, 'w'))

# Verify file
print("Verifying PDB file...")
assert work_pdb.exists(), f"{work_pdb} does not exist"
PDBFile(str(work_pdb))  # Test OpenMM can read it
print("PDB file validated")

# Import and run pipeline after path setup
sys.path.append(str(PROJECT_ROOT))
from scripts.run_with_logging import CONFIG, Pipeline

print("Starting pipeline...")
Pipeline(CONFIG, work_dir=str(WORK_DIR.parent)).run()
print("Pipeline completed successfully")
