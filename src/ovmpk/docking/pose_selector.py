from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import numpy as np # Add numpy for distance calculation

# Try importing RDKit
try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

def _get_heme_fe_coords(cfg: Dict[str, Any]) -> Optional[Tuple[float, float, float]]:
    """Extracts Heme Fe coordinates from the config if specified."""
    # Look for coordinates potentially defined in docking or prep sections
    docking_cfg = cfg.get("docking", {})
    box_cfg = docking_cfg.get("box", {})
    # Using the default center from config as the Fe coordinate source
    # This assumes the user sets the box center to the Fe coords for heme-focused docking
    center = box_cfg.get("center")
    if center and len(center) == 3:
        try:
            return (float(center[0]), float(center[1]), float(center[2]))
        except (ValueError, TypeError):
            print("[warn] Could not parse 'docking.box.center' as Fe coordinates.")
            return None

    # Add other potential config locations if needed later
    # prep_cfg = cfg.get("prep", {}).get("protein", {})
    # fe_coords = prep_cfg.get("fe_coords")
    # ...

    print("[warn] Heme Fe coordinates not found in config ('docking.box.center'). Cannot perform Fe-N distance scoring.")
    return None

def _calculate_min_fe_n_distance(mol: Chem.Mol, fe_xyz_np: np.ndarray) -> Tuple[float, Optional[np.ndarray]]:
    """Calculates the minimum distance between Fe and any N atom in the molecule."""
    min_dist = float("inf")
    closest_n_xyz = None
    if not mol:
        return min_dist, closest_n_xyz

    conf = mol.GetConformer()
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            pos = conf.GetAtomPosition(atom.GetIdx())
            n_xyz = np.array([pos.x, pos.y, pos.z])
            dist = np.linalg.norm(fe_xyz_np - n_xyz)
            if dist < min_dist:
                min_dist = dist
                closest_n_xyz = n_xyz
    return min_dist, closest_n_xyz

def pick(pose_sdf_files: List[str], cfg: Dict[str, Any]) -> Dict[str, Any]:
    """
    Selects the best pose from a list of SDF files based on configuration.

    Current Strategy:
    - If Fe coordinates are defined in config and RDKit is available:
        Select the pose (across all files/molecules) with the minimum Fe-N distance.
    - Otherwise (fallback):
        Select the first molecule from the first SDF file.

    Args:
        pose_sdf_files: A list of paths (as strings) to SDF files containing docked poses.
        cfg: The configuration dictionary.

    Returns:
        A dictionary containing information about the best pose found, e.g.,
        {'sdf_path': str, 'pose_index': int, 'score': float, 'selection_method': str}
        Returns an empty dictionary if no poses are found or RDKit is missing for scoring.
    """
    if not pose_sdf_files:
        print("[warn] No pose SDF files provided to pose selector.")
        return {}

    if not HAS_RDKIT:
        print("[warn] RDKit not found. Selecting first pose from first file as fallback.")
        return {
            'sdf_path': pose_sdf_files[0],
            'pose_index': 1, # 1-based index
            'score': None, # No score available
            'selection_method': 'first_pose_fallback_no_rdkit'
        }

    fe_coords = _get_heme_fe_coords(cfg)
    selection_method = 'first_pose_fallback_no_fe' # Default if Fe coords missing

    best_pose_info = {
        'sdf_path': pose_sdf_files[0], # Default to first
        'pose_index': 1,
        'score': float('inf'), # Use distance as score; lower is better
        'selection_method': selection_method
    }
    found_poses = False

    if fe_coords:
        fe_xyz_np = np.array(fe_coords)
        print(f"[info] Selecting pose based on minimum Fe-N distance from Fe at {fe_coords}")
        selection_method = 'min_FeN_distance'
        min_overall_dist = float('inf')

        for sdf_path_str in pose_sdf_files:
            sdf_path = Path(sdf_path_str)
            if not sdf_path.exists() or sdf_path.stat().st_size == 0:
                print(f"[warn] Pose SDF file not found or empty, skipping: {sdf_path}")
                continue

            try:
                # Ensure sanitize=True to handle potential RDKit issues gracefully
                suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=True)
                if suppl is None: # Check if supplier creation failed
                     print(f"[warn] Could not create RDKit supplier for SDF file: {sdf_path}")
                     continue

                for idx, mol in enumerate(suppl):
                    if mol is None:
                        # RDKit often fails silently on bad molecules, print index for debug
                        print(f"[warn] Could not read molecule index {idx} from {sdf_path}")
                        continue

                    found_poses = True # Mark that we found at least one processable pose
                    min_dist, _ = _calculate_min_fe_n_distance(mol, fe_xyz_np)

                    if min_dist < min_overall_dist:
                        min_overall_dist = min_dist
                        best_pose_info['sdf_path'] = str(sdf_path)
                        best_pose_info['pose_index'] = idx + 1 # Store 1-based index
                        best_pose_info['score'] = min_dist
                        best_pose_info['selection_method'] = selection_method

            except Exception as e:
                print(f"[error] Failed to process SDF file {sdf_path}: {e}")
                continue # Skip to next file on error

        if found_poses and selection_method == 'min_FeN_distance':
             # Check if a finite score was actually found
             if np.isfinite(best_pose_info['score']):
                 print(f"[info] Best pose found: Index {best_pose_info['pose_index']} in {best_pose_info['sdf_path']} with Fe-N distance = {best_pose_info['score']:.3f} Å")
             else:
                 # If no pose had a valid distance (e.g., no N atoms), fall back
                 print("[warn] No poses with valid Fe-N distance found. Falling back to first pose.")
                 selection_method = 'first_pose_fallback_no_valid_score'
                 # Revert to default first pose info (but keep path if first file was bad)
                 first_valid_sdf = next((s for s in pose_sdf_files if Path(s).exists() and Path(s).stat().st_size > 0), None)
                 if first_valid_sdf:
                     best_pose_info['sdf_path'] = first_valid_sdf
                     best_pose_info['pose_index'] = 1
                     best_pose_info['score'] = None
                 else: # Should not happen if smina_wrapper returned files, but safety check
                     print("[error] Could not find any valid SDF file for fallback.")
                     return {}
        elif not found_poses:
             print("[warn] No valid poses could be read from any SDF file.")
             return {} # Return empty if nothing was read

    else:
        # Fallback if Fe coords aren't defined - just confirm the first pose exists
        try:
            first_valid_sdf = next((s for s in pose_sdf_files if Path(s).exists() and Path(s).stat().st_size > 0), None)
            if not first_valid_sdf:
                 print("[error] Could not find any valid SDF file for fallback.")
                 return {}

            suppl = Chem.SDMolSupplier(first_valid_sdf, removeHs=False, sanitize=True)
            if suppl and suppl[0] is not None:
                found_poses = True
                best_pose_info['sdf_path'] = first_valid_sdf
                best_pose_info['pose_index'] = 1
                best_pose_info['score'] = None # No score calculated
                print(f"[info] Selecting first pose from {best_pose_info['sdf_path']} (Fe coords not specified in config).")
            else:
                 print(f"[warn] Could not read the first pose from {first_valid_sdf}.")
                 return {}
        except Exception as e:
            print(f"[error] Failed to process first SDF file {first_valid_sdf} for fallback: {e}")
            return {}


    # Update selection method in the final dict
    best_pose_info['selection_method'] = selection_method
    return best_pose_info

