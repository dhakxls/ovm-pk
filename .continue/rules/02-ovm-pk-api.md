# OVM-PK API Architecture

## Package Structure

### ovmpk.docking
- `BaseEngine`: Abstract base class for docking engines
  ```python
  class BaseEngine:
      def generate_poses(self, receptor, ligand, box, n_poses) -> List[Pose]
  ```
- `BaseScorer`: Abstract base class for scoring functions
  ```python
  class BaseScorer:
      def score_pose(self, receptor, pose, context) -> EnergyBreakdown
      def score_batch(self, receptor, poses, context) -> List[EnergyBreakdown]
  ```

### ovmpk.fetchers
- `protein_fetcher`: RCSB PDB fetching utilities
- `ligand_fetcher`: Ligand structure fetching

### ovmpk.prep
- `protein_prep`: Protein structure preparation
- `ligand_prep`: Ligand structure preparation

[Rest of API documentation continues...]
