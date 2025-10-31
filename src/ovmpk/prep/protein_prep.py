from pathlib import Path
import os
from typing import Dict, Any, Iterable, Mapping

# Define a work directory for protein prep outputs
WORK_DIR = Path("data/work/protein_prep")

# Try importing PDBFixer and OpenMM
try:
    import pdbfixer
    from openmm import app
    HAS_PDBFIXER = True
except ImportError:
    HAS_PDBFIXER = False

def prepare(paths: Dict[str, Path], cfg: Dict[str, Any]) -> Path:
    """
    Prepares a protein PDB file for docking or simulation using PDBFixer.

    Steps (if PDBFixer is available):
    1. Reads the input PDB file (typically the 'apo' structure).
    2. Finds missing residues and atoms.
    3. Adds missing heavy atoms.
    4. Adds missing hydrogens based on a specified pH.
    5. Importantly, keeps heterogens (like Heme) in the structure.
    6. Writes the processed structure to a new PDB file in the work directory.

    Args:
        paths: Dictionary containing input paths, expects "apo" key.
        cfg: Configuration dictionary, expects a 'prep.protein' section.

    Returns:
        Path to the prepared PDB file (either the fixed one or the original).
    """
    input_pdb_path = paths.get("apo")
    if input_pdb_path is None or not input_pdb_path.exists():
        raise FileNotFoundError(f"Input PDB file not found in paths dictionary or path invalid: {input_pdb_path}")

    # Get config parameters
    prep_cfg = cfg.get("prep", {}).get("protein", {})
    target_ph = float(prep_cfg.get("ph", 7.4))
    output_suffix = prep_cfg.get("output_suffix", f"_fixed_ph{target_ph}")
    run_fixer = prep_cfg.get("run_pdbfixer", True) # Option to disable fixer

    WORK_DIR.mkdir(parents=True, exist_ok=True)
    outp_pdb = WORK_DIR / f"{input_pdb_path.stem}{output_suffix}.pdb"

    if HAS_PDBFIXER and run_fixer:
        print(f"[info] Running PDBFixer on {input_pdb_path} (target pH: {target_ph})...")
        try:
            fixer = pdbfixer.PDBFixer(filename=str(input_pdb_path))

            # Find missing elements but keep heterogens like Heme
            fixer.findMissingResidues()
            fixer.findNonstandardResidues() # Identify non-standard ones
            fixer.findMissingAtoms()

            # Add missing heavy atoms, DO NOT remove heterogens
            fixer.addMissingAtoms()

            # Add missing hydrogens at the target pH
            fixer.addMissingHydrogens(target_ph)

            # Write the fixed PDB file
            with open(outp_pdb, 'w') as f:
                app.PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True) # keepIds helps maintain residue/atom numbering

            heterogen_cfg = _resolve_heterogen_config(prep_cfg.get("heterogens"))
            _normalize_heterogens(outp_pdb, heterogen_cfg)

            print(f"[info] PDBFixer complete. Output: {outp_pdb}")
            return outp_pdb

        except Exception as e:
            print(f"[warn] PDBFixer failed: {e}. Returning original PDB path: {input_pdb_path}")
            # Fallback to original PDB if fixer fails
            return input_pdb_path
    else:
        if not run_fixer:
            print("[info] PDBFixer step explicitly disabled in config.")
        else:
            print("[warn] PDBFixer library not found. Skipping protein preparation/fixing step.")
        # Return the original PDB path if fixer isn't run
        return input_pdb_path

DEFAULT_HETEROGEN_CONFIG = {
    "allowed_residues": {"HEM"},
    "rename_map": {
        "HEM": {
            "HAA": "HAA1",
            "HAAA": "HAA2",
            "HAD": "HAD1",
            "HADA": "HAD2",
            "HBA": "HBA1",
            "HBAA": "HBA2",
            "HBB": "HBB1",
            "HBBA": "HBB2",
            "HBC": "HBC1",
            "HBCA": "HBC2",
            "HBD": "HBD1",
            "HBDA": "HBD2",
            "HMA": "HMA1",
            "HMAA": "HMA2",
            "HMAB": "HMA3",
            "HMB": "HMB1",
            "HMBA": "HMB2",
            "HMBB": "HMB3",
            "HMC": "HMC1",
            "HMCA": "HMC2",
            "HMCB": "HMC3",
            "HMD": "HMD1",
            "HMDA": "HMD2",
            "HMDB": "HMD3",
        }
    },
    "drop_atoms": {
        "HEM": {"H2A", "H2D"},
    },
}


def _resolve_heterogen_config(heterogen_cfg: Dict[str, Any] | None) -> Dict[str, Any]:
    """Merge user-provided heterogen config with defaults."""

    allowed = set(DEFAULT_HETEROGEN_CONFIG["allowed_residues"])
    rename_map: Dict[str, Dict[str, str]] = {
        res: dict(mapping) for res, mapping in DEFAULT_HETEROGEN_CONFIG["rename_map"].items()
    }
    drop_atoms: Dict[str, set[str]] = {
        res: set(atoms) for res, atoms in DEFAULT_HETEROGEN_CONFIG["drop_atoms"].items()
    }

    if heterogen_cfg:
        user_allowed = heterogen_cfg.get("allowed_residues")
        if user_allowed is not None:
            allowed = {res.upper() for res in user_allowed}

        user_rename = heterogen_cfg.get("rename_map")
        if isinstance(user_rename, Mapping):
            for res, mapping in user_rename.items():
                if not isinstance(mapping, Mapping):
                    continue
                rename_map[res.upper()] = {name: new for name, new in mapping.items() if new}

        user_drop = heterogen_cfg.get("drop_atoms")
        if isinstance(user_drop, Mapping):
            for res, atoms in user_drop.items():
                drop_atoms[res.upper()] = {atom for atom in atoms}

    return {
        "allowed_residues": allowed,
        "rename_map": rename_map,
        "drop_atoms": drop_atoms,
    }


def _normalize_heterogens(pdb_path: Path, heterogen_cfg: Dict[str, Any]) -> None:
    """Normalize heterogen residue naming based on configuration."""

    allowed = heterogen_cfg.get("allowed_residues", set())
    rename_map = heterogen_cfg.get("rename_map", {})
    drop_atoms = heterogen_cfg.get("drop_atoms", {})

    try:
        lines = pdb_path.read_text().splitlines()
    except FileNotFoundError:
        return

    updated_lines = []
    dropped_serials = set()

    for line in lines:
        if len(line) < 20:
            updated_lines.append(line)
            continue

        record = line[0:6].strip()
        if record in {"ATOM", "HETATM"}:
            resname = line[17:20].strip().upper()
            if record == "HETATM" and allowed and resname not in allowed:
                dropped_serials.add(line[6:11].strip())
                continue

            residue_rename = rename_map.get(resname, {})
            residue_drop = drop_atoms.get(resname, set())

            atom_name = line[12:16].strip()
            if atom_name in residue_drop:
                dropped_serials.add(line[6:11].strip())
                continue

            new_name = residue_rename.get(atom_name)
            if new_name:
                line = f"{line[:12]}{new_name:>4}{line[16:]}"

            updated_lines.append(line)
            continue

        if record == "ANISOU" and line[6:11].strip() in dropped_serials:
            continue

        updated_lines.append(line)

    pdb_path.write_text("\n".join(updated_lines) + "\n")

def prepare_protein(pdb_file: Path, config: Dict) -> Path:
    """Process protein based on config."""
    if config['conditions']['membrane_mode'] == 'POPC':
        from .membrane import MembraneBuilder
        builder = MembraneBuilder(config)
        return builder.build(pdb_file)
    else:
        return _prepare_soluble(pdb_file, config)

def _prepare_soluble(pdb_file: Path, config: Dict) -> Path:
    """Standard soluble protein prep."""
    # Existing implementation
    return pdb_file.parent / "solvated.pdb"
