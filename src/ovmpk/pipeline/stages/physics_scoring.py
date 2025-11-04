"""
Physics-based scoring stage for the OVM-PK pipeline.

This stage uses a physics-based scoring function to evaluate the binding affinity of
small molecules to a protein target. The scoring function is based on a combination
of molecular mechanics and continuum solvation models.

The stage takes a system dictionary as input, which contains the protein-ligand
complex, and outputs an updated system dictionary with the scoring results.
"""
from typing import Dict, Any, Optional
import logging
from pathlib import Path
import numpy as np
from openmm import unit
from ovmpk.docking.scorers import NewPhysicsScorer
from ovmpk.physics.metals import MetalCenter
from ovmpk.physics.polarization import PolarizationModel
from ovmpk.physics.solvation import SolvationModel
from ovmpk.utils.typing import ArrayLike

logger = logging.getLogger(__name__)

class PhysicsScoringStage:
    """Physics-based scoring stage for the pipeline."""
    
    def __init__(self, config: Dict[str, Any]):
        """Initialize the physics scoring stage.
        
        Args:
            config: Configuration dictionary for the scoring stage
        """
        self.config = config
        self.scorer = self._setup_scorer()
    
    def _setup_scorer(self) -> NewPhysicsScorer:
        """Set up the physics-based scorer with configuration."""
        # Extract metal configuration
        metal_config = self.config.get('metal', {})
        
        # Set up metal centers
        metals = []
        for metal_spec in metal_config.get('centers', []):
            metal = MetalCenter(
                element=metal_spec['element'],
                position=metal_spec['position'],
                oxidation_state=metal_spec.get('oxidation_state'),
                coordination=metal_spec.get('coordination', 6),
                preferred_ligands=metal_spec.get('preferred_ligands', [])
            )
            metals.append(metal)
        
        # Set up physics models
        weights_cfg = self.config.get('weights', {})
        pol_cfg = self.config.get('polarization', {})
        desolv_cfg = self.config.get('desolvation', {})
        physics_config = {
            'metal': {
                'metals': metals,
                'distance_penalty': metal_config.get('distance_penalty', 50.0),
                'angle_penalty': metal_config.get('angle_penalty', 10.0)
            },
            'polarization': {
                'model': pol_cfg.get('model', self.config.get('polarization_model', 'point_dipoles')),
                'alpha_scale': pol_cfg.get('alpha_scale', self.config.get('alpha_scale', 1.0))
            },
            'desolvation': {
                'model': desolv_cfg.get('model', self.config.get('solvation_model', 'gb_like')),
                'dielectric': desolv_cfg.get('dielectric', self.config.get('dielectric', 78.5)),
                'scale': desolv_cfg.get('scale', self.config.get('solvation_scale', 1.0))
            },
            'weights': {
                'w_vdw': weights_cfg.get('w_vdw', self.config.get('w_vdw', 1.0)),
                'w_elec': weights_cfg.get('w_elec', self.config.get('w_elec', 0.6)),
                'w_pol': weights_cfg.get('w_pol', self.config.get('w_pol', 0.2)),
                'w_coord': weights_cfg.get('w_coord', self.config.get('w_coord', 1.5)),
                'w_desolv': weights_cfg.get('w_desolv', self.config.get('w_desolv', 0.5)),
                'w_tors': weights_cfg.get('w_tors', self.config.get('w_tors', 0.1))
            }
        }
        
        return NewPhysicsScorer(physics_config)
    
    def _prepare_pose(self, pose: Dict[str, Any], system: Dict[str, Any]) -> Dict[str, Any]:
        """Prepare a pose with all required properties for physics scoring.
        
        Args:
            pose: The pose to prepare
            system: The system dictionary containing additional context
            
        Returns:
            Prepared pose with all required properties
        """
        from openmm.app import PDBFile
        from openmm import Vec3
        import numpy as np
        
        # Ensure pose is a dictionary
        if not isinstance(pose, dict):
            pose = {'raw': pose}
        
        # Add required properties if missing
        if 'topology' not in pose:
            if 'topology' in system:
                pose['topology'] = system['topology']
            elif 'protein_pdb' in system:
                pdb = PDBFile(system['protein_pdb'])
                pose['topology'] = pdb.topology
                pose['positions'] = pdb.positions
        
        # Ensure positions are set
        if 'positions' not in pose and 'topology' in pose:
            pose['positions'] = [Vec3(0, 0, 0) for _ in range(pose['topology'].getNumAtoms())]
        
        # Add system reference if not present
        if 'system' not in pose and 'system' in system:
            pose['system'] = system['system']
        
        # Add default charges and polarizabilities if missing
        if 'charges' not in pose and 'topology' in pose:
            pose['charges'] = [0.0] * pose['topology'].getNumAtoms()
        
        if 'radii' not in pose and 'topology' in pose:
            pose['radii'] = [1.5] * pose['topology'].getNumAtoms()
        
        if 'polarizabilities' not in pose and 'topology' in pose:
            pose['polarizabilities'] = [0.0] * pose['topology'].getNumAtoms()
        
        return pose
    
    def run(self, system: Dict[str, Any], output_dir: Path) -> Dict[str, Any]:
        """Run the physics-based scoring stage.
        
        Args:
            system: System dictionary containing topology, coordinates, etc.
            output_dir: Directory to write output files
            
        Returns:
            Updated system dictionary with scoring results
        """
        logger.info("Running physics-based scoring stage")
        
        # Extract poses from system
        poses = system.get('poses', [])
        if not poses:
            logger.warning("No poses found for scoring")
            return system
        
        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Score each pose
        scores = []
        for i, pose in enumerate(poses):
            # Prepare the pose with all required properties
            try:
                prepared_pose = self._prepare_pose(pose, system)
            except Exception as e:
                logger.error(f"Failed to prepare pose {i}: {str(e)}")
                continue
            logger.debug(f"Scoring pose {i+1}/{len(poses)}")
            try:
                score = self.scorer.score_pose(prepared_pose)
            except Exception as e:
                logger.error(f"Error scoring pose {i}: {str(e)}")
                continue
            scores.append({
                'pose_id': i,
                'total_score': score.total,
                'breakdown': score.breakdown
            })
        
        # Sort poses by score (lower is better)
        scores.sort(key=lambda x: x['total_score'])
        
        # Update system with scoring results
        system['scoring'] = {
            'method': 'physics_based',
            'scores': scores,
            'best_pose': scores[0] if scores else None
        }
        
        # Save scoring results
        self._save_results(system, output_dir)
        
        return system
    
    def _save_results(self, system: Dict[str, Any], output_dir: Path):
        """Save scoring results to files."""
        if 'scoring' not in system:
            return
        
        # Save detailed scores to JSON
        import json
        scores_file = output_dir / 'physics_scores.json'
        with open(scores_file, 'w') as f:
            json.dump(system['scoring'], f, indent=2, default=str)
        
        # Save summary to text file
        summary_file = output_dir / 'scoring_summary.txt'
        with open(summary_file, 'w') as f:
            f.write("Physics-based Scoring Results\n")
            f.write("=" * 50 + "\n\n")
            
            if not system['scoring'].get('scores'):
                f.write("No poses were scored.\n")
                return
                
            # Write top 5 poses
            f.write("Top 5 Poses:\n")
            f.write("-" * 50 + "\n")
            for i, score in enumerate(system['scoring']['scores'][:5]):
                f.write(f"{i+1}. Pose {score['pose_id'] + 1}: {score['total_score']:.2f} kcal/mol\n")
                
                # Write energy breakdown
                if 'breakdown' in score and 'terms' in score['breakdown']:
                    f.write("   Energy breakdown:\n")
                    for term, value in score['breakdown']['terms'].items():
                        f.write(f"   - {term}: {value:.2f} kcal/mol\n")
                f.write("\n")
            
            # Write best pose info
            best = system['scoring']['best_pose']
            if best:
                f.write("\nBest Pose:\n")
                f.write("-" * 50 + "\n")
                f.write(f"Pose ID: {best['pose_id'] + 1}\n")
                f.write(f"Total Score: {best['total_score']:.2f} kcal/mol\n\n")
                
                if 'breakdown' in best and 'metal_breakdown' in best['breakdown']:
                    f.write("Metal Coordination:\n")
                    for metal_id, metal_data in best['breakdown']['metal_breakdown'].items():
                        f.write(f"- {metal_id}: {metal_data.get('energy', 0):.2f} kcal/mol\n")
        
        logger.info(f"Scoring results saved to {output_dir}")
