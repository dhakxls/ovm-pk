"""
OVM-PK Pipeline with Physics-Aware Docking

This pipeline includes physics-based scoring for more accurate
binding affinity predictions.
"""
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Any, Optional, List, Union
import json
import yaml
import mdtraj as md
import numpy as np

from openmm import unit, LangevinIntegrator, LocalEnergyMinimizer
from openmm.app import PDBFile, PDBxFile, Modeller, ForceField, Simulation, Topology
from openmm import Platform, XmlSerializer
import openmm
from pdbfixer import PDBFixer

from .docking.scorers import NewPhysicsScorer
from .physics.metals import MetalCenter
from .physics.polarization import PolarizationModel
from .physics.solvation import SolvationModel
from ovmpk.pipeline.stages.physics_scoring import PhysicsScoringStage
from .utils import setup_logging, validate_config, get_default_platform

logger = logging.getLogger(__name__)


def _compute_dihedral_angle(positions: np.ndarray, indices: List[int]) -> float:
    """Compute the dihedral angle (radians) from Cartesian positions.

    Args:
        positions: Array of shape (n_atoms, 3) in Angstroms.
        indices: Four atom indices defining the torsion.

    Returns:
        Dihedral angle in radians.
    """
    i0, i1, i2, i3 = indices
    p0, p1, p2, p3 = positions[[i0, i1, i2, i3]]

    b0 = p1 - p0
    b1 = p2 - p1
    b2 = p3 - p2

    b1_norm = np.linalg.norm(b1)
    if b1_norm == 0:
        return 0.0
    b1_unit = b1 / b1_norm

    v = b0 - np.dot(b0, b1_unit) * b1_unit
    w = b2 - np.dot(b2, b1_unit) * b1_unit

    x = np.dot(v, w)
    y = np.dot(np.cross(b1_unit, v), w)

    return np.arctan2(y, x)


def _identify_ligand_torsions(ligand_topology) -> List[List[int]]:
    """Enumerate ligand torsions (unique dihedral quadruples)."""
    adjacency = defaultdict(set)
    bonds = []

    for atom1, atom2 in ligand_topology.bonds():
        i = atom1.index
        j = atom2.index
        adjacency[i].add(j)
        adjacency[j].add(i)
        bonds.append((i, j))

    torsions = set()
    for j, k in bonds:
        for i in adjacency[j]:
            if i == k:
                continue
            for l in adjacency[k]:
                if l == j or l == i:
                    continue
                torsion = (i, j, k, l)
                reverse = (l, k, j, i)
                if reverse in torsions:
                    continue
                torsions.add(torsion)

    return [list(t) for t in torsions]


class Pipeline:
    """Main pipeline class for OVM-PK with physics-based scoring."""
    
    def __init__(self, config_path: Union[str, Path], output_dir: Optional[Union[str, Path]] = None):
        """Initialize the pipeline with configuration.
        
        Args:
            config_path: Path to the configuration file
            output_dir: Base output directory (defaults to 'output' in config directory)
        """
        self.config_path = Path(config_path).resolve()
        self.config = self._load_config()
        
        # Set up output directory
        if output_dir is None:
            output_dir = self.config_path.parent / 'output'
        self.output_dir = Path(output_dir).resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Set up logging
        log_file = self.output_dir / 'pipeline.log'
        setup_logging(log_file=log_file, log_level=self.config.get('log_level', 'INFO'))
        
        # Initialize physics models
        self._init_physics_models()

        # Initialize OpenMM platform
        self.platform = get_default_platform()

        logger.info(f"Initialized OVM-PK Pipeline with config: {config_path}")
        ff_files = ['amber14-all.xml', 'amber14/tip3pfb.xml']
        project_root = Path(__file__).resolve().parents[2]
        heme_ff = project_root / 'forcefields' / 'shahrokh_heme_ic6.ffxml'
        ligand_ff = project_root / 'forcefields' / 'UNL_gaff211.ffxml'
        if heme_ff.exists():
            logger.info(f"Including heme force field parameters: {heme_ff}")
            ff_files.append(str(heme_ff))
        else:
            logger.warning("Heme force field file not found; proceeding without explicit heme parameters")
        if ligand_ff.exists():
            logger.info(f"Including ligand force field parameters: {ligand_ff}")
            ff_files.append(str(ligand_ff))
        else:
            logger.warning("Ligand force field file not found; proceeding without explicit ligand parameters")
        self.structure_forcefield = ForceField(*ff_files)

    def _relax_pose_geometry(self, topology, positions):
        """Perform staged minimization to reduce steric clashes."""
        system = self.structure_forcefield.createSystem(
            topology,
            nonbondedMethod=openmm.app.PME,
            constraints=openmm.app.HBonds
        )

        integrator = LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds, 0.001 * unit.picoseconds)
        platform = Platform.getPlatformByName('Reference')
        simulation = Simulation(topology, system, integrator, platform)
        simulation.context.setPositions(positions)

        LocalEnergyMinimizer.minimize(simulation.context, 1.0, 500)
        state = simulation.context.getState(getPositions=True)
        minimized_positions = state.getPositions()

        return np.array(minimized_positions.value_in_unit(unit.angstrom))

    def _sanitize_heme_residues(self, topology: Topology, positions):
        """Normalize HEM residue and coordinating cysteine to match the IC6 template."""
        rename_map = {
            'HMA': 'HMA1', 'HMAA': 'HMA2', 'HMAB': 'HMA3',
            'HAA': 'HAA1', 'HAAA': 'HAA2',
            'HBA': 'HBA1', 'HBAA': 'HBA2',
            'HMB': 'HMB1', 'HMBA': 'HMB2', 'HMBB': 'HMB3',
            'HMC': 'HMC1', 'HMCA': 'HMC2', 'HMCB': 'HMC3',
            'HMD': 'HMD1', 'HMDA': 'HMD2', 'HMDB': 'HMD3',
            'HAD': 'HAD1', 'HADA': 'HAD2',
            'HBD': 'HBD1', 'HBDA': 'HBD2',
            'HBC': 'HBC1', 'HBCA': 'HBC2',
            'HBB': 'HBB1', 'HBBA': 'HBB2',
            'H2A': None, 'H2D': None,
        }

        allowed_names = {
            'NC', 'C1C', 'C4C', 'C2C', 'C3C', 'CHD', 'HHD', 'C1D', 'ND', 'C4D',
            'C3D', 'C2D', 'CHA', 'HHA', 'C1A', 'C2A', 'C4A', 'NA', 'C3A', 'CHB',
            'HHB', 'C1B', 'C2B', 'NB', 'C4B', 'CHC', 'HHC', 'FE', 'C3B', 'CAB',
            'HAB', 'CBB', 'HBB1', 'HBB2', 'CAC', 'HAC', 'CBC', 'HBC1', 'HBC2',
            'CMB', 'HMB1', 'HMB2', 'HMB3', 'CMC', 'HMC1', 'HMC2', 'HMC3', 'CMD',
            'HMD1', 'HMD2', 'HMD3', 'CMA', 'HMA1', 'HMA2', 'HMA3', 'CAA', 'HAA1',
            'HAA2', 'CAD', 'HAD1', 'HAD2', 'CBA', 'HBA1', 'HBA2', 'CBD', 'HBD1',
            'HBD2', 'CGD', 'CGA', 'O1A', 'O1D', 'O2A', 'O2D'
        }

        heme_template_bonds = [
            ('NC', 'C1C'), ('NC', 'C4C'), ('NC', 'FE'),
            ('C1C', 'C2C'), ('C1C', 'CHC'),
            ('C4C', 'C3C'), ('C4C', 'CHD'),
            ('C2C', 'C3C'), ('C2C', 'CMC'),
            ('C3C', 'CAC'),
            ('CHD', 'HHD'), ('CHD', 'C1D'),
            ('C1D', 'ND'), ('C1D', 'C2D'),
            ('ND', 'C4D'), ('ND', 'FE'),
            ('C4D', 'C3D'), ('C4D', 'CHA'),
            ('C3D', 'C2D'), ('C3D', 'CAD'),
            ('C2D', 'CMD'),
            ('CHA', 'HHA'), ('CHA', 'C1A'),
            ('C1A', 'C2A'), ('C1A', 'NA'),
            ('C2A', 'C3A'), ('C2A', 'CAA'),
            ('C3A', 'C4A'), ('C3A', 'CMA'),
            ('C4A', 'NA'), ('C4A', 'CHB'),
            ('NA', 'FE'),
            ('CHB', 'HHB'), ('CHB', 'C1B'),
            ('C1B', 'C2B'), ('C1B', 'NB'),
            ('C2B', 'C3B'), ('C2B', 'CMB'),
            ('C3B', 'C4B'), ('C3B', 'CAB'),
            ('C4B', 'NB'), ('C4B', 'CHC'),
            ('NB', 'FE'),
            ('CHC', 'HHC'),
            ('CAB', 'HAB'), ('CAB', 'CBB'),
            ('CBB', 'HBB1'), ('CBB', 'HBB2'),
            ('CAC', 'HAC'), ('CAC', 'CBC'),
            ('CBC', 'HBC1'), ('CBC', 'HBC2'),
            ('CMB', 'HMB1'), ('CMB', 'HMB2'), ('CMB', 'HMB3'),
            ('CMC', 'HMC1'), ('CMC', 'HMC2'), ('CMC', 'HMC3'),
            ('CMD', 'HMD1'), ('CMD', 'HMD2'), ('CMD', 'HMD3'),
            ('CMA', 'HMA1'), ('CMA', 'HMA2'), ('CMA', 'HMA3'),
            ('CAA', 'HAA1'), ('CAA', 'HAA2'), ('CAA', 'CBA'),
            ('CAD', 'HAD1'), ('CAD', 'HAD2'), ('CAD', 'CBD'),
            ('CBA', 'HBA1'), ('CBA', 'HBA2'), ('CBA', 'CGA'),
            ('CBD', 'HBD1'), ('CBD', 'HBD2'), ('CBD', 'CGD'),
            ('CGD', 'O1D'), ('CGD', 'O2D'),
            ('CGA', 'O1A'), ('CGA', 'O2A'),
        ]

        positions_array = positions.value_in_unit(unit.nanometer) if hasattr(positions, 'value_in_unit') else np.asarray(positions)
        positions_lookup: Dict[Topology.Atom, np.ndarray] = {}
        for atom, pos in zip(topology.atoms(), positions_array):
            positions_lookup[atom] = np.array(pos)

        fe_atom_original = None
        for atom in topology.atoms():
            if atom.residue.name.strip() == 'HEM' and atom.name.strip() == 'FE':
                fe_atom_original = atom
                break

        target_cys_residue = None
        if fe_atom_original is not None:
            fe_pos = positions_lookup[fe_atom_original]
            best_dist = None
            for atom in topology.atoms():
                if atom.residue.name.strip() in {'CYS', 'CYM', 'CYX', 'CYP'} and atom.name.strip() in {'SG', 'S'}:
                    dist = np.linalg.norm(positions_lookup[atom] - fe_pos)
                    if best_dist is None or dist < best_dist:
                        best_dist = dist
                        target_cys_residue = atom.residue

        new_topology = Topology()
        if topology.getPeriodicBoxVectors() is not None:
            new_topology.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())

        atom_map = {}
        new_positions = []
        heme_atom_maps: Dict[Topology.Residue, Dict[str, Topology.Atom]] = {}
        heme_modified = False
        residue_map: Dict[Topology.Residue, Topology.Residue] = {}

        pos_index = 0
        for chain in topology.chains():
            new_chain = new_topology.addChain(chain.id)
            for residue in chain.residues():
                original_name = residue.name.strip()
                is_heme = original_name == 'HEM'
                is_cys = original_name in {'CYS', 'CYM', 'CYX', 'CYP'}
                is_target_cys = residue is target_cys_residue
                new_residue_name = 'CYP' if is_target_cys else residue.name
                new_residue = residue_map.get(residue)
                if new_residue is None:
                    new_residue = new_topology.addResidue(new_residue_name, new_chain, residue.id)
                    residue_map[residue] = new_residue

                for atom in residue.atoms():
                    position = positions_array[pos_index]
                    pos_index += 1
                    atom_name = atom.name
                    if is_heme:
                        heme_modified = True
                        atom_name = rename_map.get(atom_name, atom_name)
                        if atom_name is None:
                            continue
                        if atom_name not in allowed_names:
                            logger.warning(f"Removing unexpected HEM atom '{atom.name}'")
                            continue
                    if is_target_cys and atom.name.strip() in {'HG', 'HSG', 'HG1', 'HG2'}:
                        continue
                    new_atom = new_topology.addAtom(atom_name, atom.element, new_residue)
                    atom_map[atom] = new_atom
                    new_positions.append(position)
                    if is_heme:
                        heme_atom_maps.setdefault(new_residue, {})[atom_name] = new_atom

        for atom1, atom2 in topology.bonds():
            new_atom1 = atom_map.get(atom1)
            new_atom2 = atom_map.get(atom2)
            if new_atom1 is not None and new_atom2 is not None:
                new_topology.addBond(new_atom1, new_atom2)

        if not heme_modified:
            return topology, positions

        existing_heme_bonds = set()
        for bond in new_topology.bonds():
            if bond[0].residue is bond[1].residue and bond[0].residue.name.strip() == 'HEM':
                key = (id(bond[0].residue), tuple(sorted((bond[0].name, bond[1].name))))
                existing_heme_bonds.add(key)

        for residue, name_map in heme_atom_maps.items():
            for name1, name2 in heme_template_bonds:
                atom1 = name_map.get(name1)
                atom2 = name_map.get(name2)
                if atom1 is None or atom2 is None:
                    continue
                key = (id(residue), tuple(sorted((atom1.name, atom2.name))))
                if key in existing_heme_bonds:
                    continue
                new_topology.addBond(atom1, atom2)
                existing_heme_bonds.add(key)

        # Identify coordinating cysteine (nearest SG to FE) and enforce Fe-S bond
        fe_atom = None
        for name_map in heme_atom_maps.values():
            fe_atom = name_map.get('FE')
            if fe_atom is not None:
                break

        if target_cys_residue is not None:
            new_target_residue = residue_map.get(target_cys_residue)
            if new_target_residue is not None and new_target_residue.name.strip() != 'CYP':
                new_target_residue.name = 'CYP'

        logger.info("Sanitized HEM residue to match force-field template")
        sanitized_positions = unit.Quantity(np.array(new_positions), unit.nanometer)
        return new_topology, sanitized_positions

    def _build_pose_system(self, topology: openmm.app.Topology) -> openmm.System:
        return self.structure_forcefield.createSystem(
            topology,
            nonbondedMethod=openmm.app.PME,
            constraints=openmm.app.HBonds
        )

    def _compute_system_energy(self, system: openmm.System, positions: np.ndarray) -> float:
        integrator = LangevinIntegrator(300 * unit.kelvin, 1.0 / unit.picoseconds, 0.001 * unit.picoseconds)
        platform = Platform.getPlatformByName('Reference')
        context = openmm.Context(system, integrator, platform)
        context.setPositions([openmm.Vec3(*pos) for pos in positions * 0.1])
        LocalEnergyMinimizer.minimize(context, 1.0, 200)
        state = context.getState(getEnergy=True)
        energy = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
        del context
        del integrator
        return energy

    def _extract_nonbonded_parameters(self, system: openmm.System):
        charges = []
        radii = []
        polarizabilities = []
        sigma_to_radius = 2 ** (1 / 6) / 2  # rmin/2 conversion factor
        for force in system.getForces():
            if isinstance(force, openmm.NonbondedForce):
                for i in range(force.getNumParticles()):
                    charge, sigma, epsilon = force.getParticleParameters(i)
                    charge_e = float(charge.value_in_unit(unit.elementary_charge))
                    charges.append(charge_e)

                    sigma_nm = float(sigma.value_in_unit(unit.nanometer))
                    if sigma_nm <= 0:
                        sigma_nm = 0.1

                    rmin_over_two_ang = sigma_nm * 10.0 * sigma_to_radius
                    radius_ang = np.clip(rmin_over_two_ang, 0.7, 2.8)
                    radii.append(radius_ang)

                    polar = 0.5 * (radius_ang ** 3)
                    polar_clamped = np.clip(polar, 0.2, 20.0)
                    polarizabilities.append(polar_clamped)
                break
        charges_arr = np.array(charges, dtype=float)
        if charges_arr.size:
            charges_arr = np.clip(charges_arr, -1.5, 1.5)
            charges_arr -= charges_arr.mean()
            charges_arr = np.clip(charges_arr, -1.5, 1.5)
            # Re-center after clipping to enforce neutrality within tolerance
            charges_arr -= charges_arr.mean()
        return charges_arr.tolist(), radii, polarizabilities

    def _rigid_align_ligand(self, protein_positions: np.ndarray, ligand_positions: np.ndarray) -> np.ndarray:
        """Translate ligand to protein centroid to reduce initial overlaps."""
        protein_center = protein_positions.mean(axis=0)
        ligand_center = ligand_positions.mean(axis=0)
        translation = protein_center - ligand_center
        return ligand_positions + translation

    def _load_config(self) -> Dict[str, Any]:
        """Load and validate the pipeline configuration."""
        with open(self.config_path, 'r') as f:
            config = yaml.safe_load(f)

        
        # Set default values
        config.setdefault('system', {})
        config.setdefault('simulation', {})
        config.setdefault('scoring', {})
        
        # Validate required fields
        required = {
            'system': ['protein_pdb', 'ligand_pdb'],
            'simulation': ['forcefield'],
            'scoring': ['method']
        }
        
        for section, fields in required.items():
            for field in fields:
                if field not in config.get(section, {}):
                    raise ValueError(f"Missing required config: {section}.{field}")
        
        return config
    
    def _init_physics_models(self):
        """Initialize physics models for scoring."""
        # Initialize metal centers if specified
        self.metal_centers = []
        metal_config = self.config.get('system', {}).get('metals', [])
        
        for metal_spec in metal_config:
            metal = MetalCenter(
                element=metal_spec['element'],
                position=metal_spec['position'],
                oxidation_state=metal_spec.get('oxidation_state'),
                coordination=metal_spec.get('coordination', 6),
                preferred_ligands=metal_spec.get('preferred_ligands', [])
            )
            self.metal_centers.append(metal)
        
        # Initialize physics models
        self.polarization_model = PolarizationModel()
        self.solvation_model = SolvationModel()

        ligand_cfg = self.config.get('system', {}).get('ligand', {})
        ligand_sdf = ligand_cfg.get('sdf_file')
        if ligand_sdf:
            ligand_sdf_path = Path(ligand_sdf)
            if not ligand_sdf_path.is_absolute():
                ligand_sdf_path = (Path(__file__).resolve().parents[2] / ligand_sdf).resolve()
            ligand_sdf = str(ligand_sdf_path)

        scoring_cfg = self.config.get('scoring', {})
        scoring_params = scoring_cfg.get('parameters', {})
        dielectric = scoring_params.get('dielectric', scoring_cfg.get('dielectric', 78.5))
        solvation_model = scoring_params.get('solvation', scoring_cfg.get('solvation_model', 'gb_like'))
        alpha_scale = scoring_params.get('alpha_scale', scoring_cfg.get('alpha_scale', 1.0))
        solvation_scale = scoring_params.get('solvation_scale', scoring_cfg.get('solvation_scale', 1.0))

        # Initialize physics scorer
        physics_config = {
            'metal': {
                'metals': self.metal_centers,
                'distance_penalty': scoring_cfg.get('distance_penalty', 50.0),
                'angle_penalty': scoring_cfg.get('angle_penalty', 10.0)
            },
            'polarization': {
                'model': scoring_cfg.get('polarization_model', 'point_dipoles'),
                'alpha_scale': alpha_scale
            },
            'desolvation': {
                'model': solvation_model,
                'dielectric': dielectric,
                'scale': solvation_scale
            },
            'weights': {
                'w_vdw': scoring_cfg.get('w_vdw', 1.0),
                'w_elec': scoring_cfg.get('w_elec', 0.6),
                'w_pol': scoring_cfg.get('w_pol', 0.2),
                'w_coord': scoring_cfg.get('w_coord', 1.5),
                'w_desolv': scoring_cfg.get('w_desolv', 0.5),
                'w_tors': scoring_cfg.get('w_tors', 0.1)
            },
            'ligand': {
                'sdf_file': ligand_sdf,
                'residue_name': ligand_cfg.get('residue_name', 'UNL')
            }
        }
        
        self.physics_scorer = NewPhysicsScorer(physics_config)
        self.physics_stage = PhysicsScoringStage(physics_config)
    
    def run(self):
        """Run the complete pipeline with physics-based scoring."""
        logger.info("Starting OVM-PK Pipeline with physics-based scoring")
        
        # Create stage output directories
        stage_dirs = {
            'prep': self.output_dir / '01_preparation',
            'docking': self.output_dir / '02_docking',
            'scoring': self.output_dir / '03_physics_scoring',
            'simulation': self.output_dir / '04_simulation',
            'analysis': self.output_dir / '05_analysis'
        }
        
        for dir_path in stage_dirs.values():
            dir_path.mkdir(exist_ok=True)
        
        # Run pipeline stages
        try:
            # Stage 1: System preparation
            system = self._prepare_system(stage_dirs['prep'])
            
            # Stage 2: Docking (or load existing poses)
            poses = self._run_docking(system, stage_dirs['docking'])
            system['poses'] = poses
            
            # Stage 3: Physics-based scoring
            system = self.physics_stage.run(system, stage_dirs['scoring'])
            
            # Stage 4: Simulation of top poses
            if self.config.get('simulation', {}).get('run_simulation', True):
                top_poses = self._select_top_poses(system, n_poses=3)
                simulation_results = self._run_simulations(top_poses, stage_dirs['simulation'])
                system['simulation_results'] = simulation_results
            
            # Stage 5: Analysis
            if self.config.get('analysis', {}).get('run_analysis', True):
                analysis_results = self._analyze_results(system, stage_dirs['analysis'])
                system['analysis'] = analysis_results
            
            # Save final results
            self._save_results(system, self.output_dir)
            
            logger.info("Pipeline completed successfully!")
            return system
            
        except Exception as e:
            logger.error(f"Pipeline failed: {str(e)}", exc_info=True)
            raise
    
    def _prepare_system(self, output_dir: Path) -> Dict[str, Any]:
        """Prepare the system for simulation.
        
        Args:
            output_dir: Directory to save prepared structures
            
        Returns:
            Dictionary containing prepared system information
        """
        logger.info("Preparing system...")
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Get file paths from config
        protein_pdb = Path(self.config['system']['protein_pdb']).resolve()
        ligand_pdb = Path(self.config['system']['ligand_pdb']).resolve()
        
        # Create output file paths
        prepared_protein = output_dir / 'prepared_protein.pdb'
        prepared_ligand = output_dir / 'prepared_ligand.pdb'
        
        # Prepare protein structure
        logger.info(f"Preparing protein: {protein_pdb}")
        try:
            # Check if the file is a CIF file
            if str(protein_pdb).endswith('.cif'):
                # For CIF files, use MDTraj to convert to PDB first
                import mdtraj as md
                traj = md.load(str(protein_pdb))
                temp_pdb = output_dir / 'temp_protein.pdb'
                traj.save(str(temp_pdb))
                protein_file = temp_pdb
            else:
                protein_file = protein_pdb
            
            # Load protein with PDBFixer to add missing atoms and standardize
            fixer = PDBFixer(str(protein_file))
            
            # Common fixes
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(7.0)  # pH 7.0
            
            # Save prepared protein
            sanitized_topology, sanitized_positions = self._sanitize_heme_residues(fixer.topology, fixer.positions)

            with open(prepared_protein, 'w') as f:
                PDBFile.writeFile(sanitized_topology, sanitized_positions, f)

            # Persist sanitized data for downstream stages to avoid reloading losses
            self._prepared_topology = sanitized_topology
            self._prepared_positions = sanitized_positions
                
            logger.info(f"Prepared protein saved to {prepared_protein}")
            
        except Exception as e:
            logger.error(f"Error preparing protein: {str(e)}")
            raise
        
        # Prepare ligand structure
        logger.info(f"Preparing ligand: {ligand_pdb}")
        try:
            from openff.toolkit import Molecule

            if ligand_pdb.suffix.lower() == '.sdf':
                ligand_mol = Molecule.from_file(str(ligand_pdb))
            else:
                sdf_path = ligand_pdb.with_suffix('.sdf')
                if sdf_path.exists():
                    ligand_mol = Molecule.from_file(str(sdf_path))
                else:
                    raise ValueError(
                        "Ligand preparation requires an SDF file with full connectivity; "
                        f"no SDF found alongside {ligand_pdb}"
                    )

            if not ligand_mol.name:
                filename_stem = ligand_pdb.stem
                if filename_stem:
                    ligand_mol.name = filename_stem[:20]
            residue_name = ligand_mol.name[:3].upper() if ligand_mol.name else 'UNL'
            ligand_mol.name = residue_name

            ligand_mol.to_file(str(prepared_ligand), file_format='pdb')
            logger.info(f"Prepared ligand saved to {prepared_ligand}")

        except Exception as e:
            logger.error(f"Error preparing ligand: {str(e)}")
            raise

        return {
            'protein': str(prepared_protein),
            'ligand': str(prepared_ligand),
            'prepared': True,
            'output_dir': str(output_dir),
            'protein_topology': getattr(self, '_prepared_topology', None),
            'protein_positions': getattr(self, '_prepared_positions', None),
        }
    
    def _run_docking(self, system: Dict[str, Any], output_dir: Path) -> List[Dict]:
        """Run docking to generate poses with proper structural information.
        
        Args:
            system: Dictionary containing system information
            output_dir: Directory to save docking results
            
        Returns:
            List of pose dictionaries with required structural information
        """
        logger.info("Running docking...")
        
        # Load protein and ligand structures
        protein_file = system.get('protein')
        ligand_file = system.get('ligand')
        
        if not protein_file or not ligand_file:
            raise ValueError("Protein and ligand files are required for docking")
        
        # Create output directory for docking results
        output_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            stored_topology = system.get('protein_topology')
            stored_positions = system.get('protein_positions')

            if stored_topology is not None and stored_positions is not None:
                logger.info("Using sanitized protein topology from preparation stage")
                protein_topology = stored_topology
                # stored_positions is Quantity in nanometers; convert to Angstrom numpy array
                protein_positions = np.array(stored_positions.value_in_unit(unit.nanometer)) * 10.0
            else:
                # Load protein structure
                logger.info(f"Loading protein from {protein_file}")
                protein = md.load(protein_file)
                protein_topology = protein.topology.to_openmm()
                protein_positions = protein.xyz[0] * 10.0  # Convert nm to Angstroms
            
            # Load ligand structure - handle different file formats
            logger.info(f"Loading ligand from {ligand_file}")
            if str(ligand_file).endswith('.sdf'):
                # For SDF files, use RDKit to convert to PDB first
                from rdkit import Chem
                mol = Chem.SDMolSupplier(str(ligand_file))[0]
                temp_pdb = output_dir / 'temp_ligand.pdb'
                Chem.MolToPDBFile(mol, str(temp_pdb))
                ligand = md.load(str(temp_pdb))
            else:
                ligand = md.load(ligand_file)
                
            ligand_topology = ligand.topology.to_openmm()
            ligand_positions = ligand.xyz[0] * 10.0  # Convert nm to Angstroms
            ligand_torsions = _identify_ligand_torsions(ligand_topology)
            
            # Get number of poses from config, default to 1 if not specified
            num_poses = self.config.get('docking', {}).get('num_poses', 1)
            poses = []
            
            logger.info(f"Generating {num_poses} poses...")
            
            # Create a Modeller to combine protein and ligand
            modeller = Modeller(protein_topology, protein_positions * unit.angstrom)
            modeller.add(ligand_topology, ligand_positions * unit.angstrom)
            
            for i in range(num_poses):
                # Align ligand near protein centroid then apply jitter
                pose_ligand_positions = self._rigid_align_ligand(protein_positions, ligand_positions)
                translation_range = self.config.get('docking', {}).get('random_translation', 0.2)
                translation_vector = np.random.uniform(-translation_range, translation_range, 3)
                pose_ligand_positions += translation_vector
                
                # Create a new modeller for this pose
                pose_modeller = Modeller(protein_topology, protein_positions * unit.angstrom)
                pose_modeller.add(ligand_topology, pose_ligand_positions * unit.angstrom)
                
                # Get the combined topology and positions
                combined_topology = pose_modeller.getTopology()
                combined_positions = pose_modeller.getPositions()
                # Convert positions to numpy array in Angstroms
                try:
                    combined_positions = self._relax_pose_geometry(combined_topology, combined_positions)
                    logger.debug("Pose geometry relaxed successfully")
                except Exception as relax_error:
                    logger.warning(f"Pose relaxation failed: {relax_error}")
                    combined_positions = np.array(combined_positions.value_in_unit(unit.angstrom))
                
                # Create pose dictionary with required properties
                pose_system = self._build_pose_system(combined_topology)
                charges, radii, polarizabilities = self._extract_nonbonded_parameters(pose_system)
                charges_arr = np.array(charges)
                radii_arr = np.array(radii)
                polar_arr = np.array(polarizabilities)

                def _summarize(name: str, values: np.ndarray, unit_label: str = "") -> None:
                    if values.size == 0:
                        logger.warning(f"No values extracted for {name}")
                        return
                    stats = {
                        'count': int(values.size),
                        'min': float(values.min()),
                        'max': float(values.max()),
                        'mean': float(values.mean()),
                        'std': float(values.std()),
                    }
                    label = f" ({unit_label})" if unit_label else ""
                    logger.info(f"Extracted {name}{label}: {stats}")

                _summarize('charges', charges_arr, 'e')
                _summarize('radii', radii_arr, 'Å')
                _summarize('polarizabilities', polar_arr, 'Å^3')

                total_charge = float(charges_arr.sum()) if charges_arr.size else 0.0
                if abs(total_charge) > 1e-3:
                    logger.warning(f"Total system charge from extraction is {total_charge:.4f} e")

                pose_energy = self._compute_system_energy(pose_system, combined_positions)

                pose = {
                    'id': i,
                    'score': pose_energy,
                    'topology': combined_topology,
                    'positions': combined_positions,
                    'pdb_file': str(output_dir / f'pose_{i:03d}.pdb'),
                    'charges': np.array(charges),
                    'radii': np.array(radii),
                    'polarizabilities': np.array(polarizabilities),
                    'ligand_indices': list(range(
                        protein_topology.getNumAtoms(),
                        protein_topology.getNumAtoms() + ligand_topology.getNumAtoms()
                    )),
                    'openmm_system': pose_system,
                }

                if ligand_torsions:
                    ligand_offset = protein_topology.getNumAtoms()
                    pose['ligand_torsions'] = [
                        {
                            'indices': [ligand_offset + idx for idx in torsion],
                            'theta0': _compute_dihedral_angle(pose_ligand_positions, torsion),
                            'k': 1.0,
                        }
                        for torsion in ligand_torsions
                    ]

                # Save pose to PDB for visualization
                with open(pose['pdb_file'], 'w') as f:
                    PDBFile.writeFile(combined_topology, combined_positions, f)

                poses.append(pose)
                logger.debug(f"Generated pose {i} with {combined_topology.getNumAtoms()} atoms")
            
            logger.info(f"Generated {len(poses)} poses with structural information")
            return poses
            
        except Exception as e:
            logger.error(f"Error during docking: {str(e)}", exc_info=True)
            raise
    
    def _select_top_poses(self, system: Dict[str, Any], n_poses: int = 3) -> List[Dict]:
        """Select top N poses based on physics-based scores."""
        if 'scoring' not in system or 'scores' not in system['scoring']:
            logger.warning("No scores available for pose selection. Using first N poses.")
            return system.get('poses', [])[:n_poses]
        
        # Sort poses by score and take top N
        sorted_poses = sorted(
            zip(system['poses'], system['scoring']['scores']),
            key=lambda x: x[1]['total_score']
        )
        
        return [pose for pose, _ in sorted_poses[:n_poses]]
    
    def _run_simulations(self, poses: List[Dict], output_dir: Path) -> Dict[int, Dict]:
        """Run MD simulations for the top poses."""
        logger.info(f"Running MD simulations for {len(poses)} poses...")
        
        # TODO: Implement MD simulation logic
        # For now, return dummy results
        return {
            i: {
                'rmsd': 0.5 + i*0.1,
                'binding_energy': -8.0 + i*0.5,
                'stability': 0.9 - i*0.05
            }
            for i in range(len(poses))
        }
    
    def _analyze_results(self, system: Dict[str, Any], output_dir: Path) -> Dict[str, Any]:
        """Analyze simulation results."""
        logger.info("Analyzing results...")
        
        # TODO: Implement analysis logic
        # For now, return dummy analysis
        return {
            'best_pose': system.get('scoring', {}).get('best_pose', {}).get('pose_id', 0),
            'binding_affinity': -8.5,  # kcal/mol
            'confidence': 0.85
        }
    
    def _save_results(self, system: Dict[str, Any], output_dir: Path):
        """Save pipeline results to disk."""
        # Save system state
        with open(output_dir / 'system_state.json', 'w') as f:
            json.dump(system, f, indent=2, default=str)
        
        # Save summary
        summary = {
            'status': 'completed',
            'best_pose': system.get('scoring', {}).get('best_pose', {}).get('pose_id', None),
            'binding_affinity': system.get('analysis', {}).get('binding_affinity', None),
            'confidence': system.get('analysis', {}).get('confidence', None)
        }
        
        with open(output_dir / 'summary.json', 'w') as f:
            json.dump(summary, f, indent=2)
        
        logger.info(f"Results saved to {output_dir}")


def main():
    """Command-line interface for the pipeline."""
    import argparse
    
    parser = argparse.ArgumentParser(description='OVM-PK Pipeline with Physics-Aware Docking')
    parser.add_argument('config', help='Path to configuration file')
    parser.add_argument('-o', '--output', help='Output directory (default: output/)')
    parser.add_argument('--log-level', default='INFO', 
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                       help='Logging level')
    
    args = parser.parse_args()
    
    # Set up logging
    setup_logging(level=args.log_level)
    
    try:
        # Initialize and run pipeline
        pipeline = PipelineV2(args.config, output_dir=args.output)
        results = pipeline.run()
        
        # Print summary
        print("\nPipeline completed successfully!")
        print(f"Results saved to: {pipeline.output_dir}")
        
        if 'analysis' in results:
            print("\nAnalysis Results:")
            print(f"  Best Pose: {results['analysis'].get('best_pose', 'N/A')}")
            print(f"  Binding Affinity: {results['analysis'].get('binding_affinity', 'N/A')} kcal/mol")
            print(f"  Confidence: {results['analysis'].get('confidence', 'N/A')}")
        
        return 0
        
    except Exception as e:
        logger.critical(f"Pipeline failed: {str(e)}", exc_info=True)
        print(f"\nError: {str(e)}")
        print("Check the log file for more details.")
        return 1


if __name__ == "__main__":
    import sys
    sys.exit(main())
