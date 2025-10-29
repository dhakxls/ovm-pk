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

### ovmpk.utils
- `logging`: Structured logging utilities
- `toolchain`: Build and execution helpers
- `smiles_fetch`: SMILES format handling

## Configuration Schema

```yaml
docking:
  engine: smina             # or other engine
  scorer: baseline_vina     # or new_physics
  box:
    auto_center: heme_fe    # or explicit coords
    size: [12.0, 12.0, 12.0]

scoring:
  new_physics:
    metal:
      enable: true
      metals: [FE]
      preferred_bonds:
        - donor: N
          ideal: 2.1
          k: 50.0
    polarization:
      enable: true
      model: point_dipoles

benchmark:
  dataset_csv: data/bench/kiset.csv
  temperature_default: 298.15
  metrics_out: results/bench/metrics.csv
```

## Common Usage Patterns

### Fetching Structures
```python
from ovmpk.fetchers import fetch_protein, fetch_ligand

# Fetch protein structure
fetch_protein(pdbid="5VCC", out_dir=work_dir)

# Fetch ligand structure
fetch_ligand(name="ketoconazole", out_dir=work_dir)
```

### Structure Preparation
```python
from ovmpk.prep import prepare_protein, prepare_ligand

# Prepare protein
prepare_protein(input_pdb, output_pdb, ph=7.4)

# Prepare ligand
prepare_ligand(input_sdf, output_sdf, ph=7.4)
```

### Docking Operations
```python
from ovmpk.docking import BaseEngine, BaseScorer

# Using engine
poses = engine.generate_poses(receptor, ligand, box, n_poses=20)

# Scoring poses
scores = scorer.score_batch(receptor, poses, context)
```

## Testing Architecture

### Test Dependencies
1. `test_fetchers.py` - Base fetching tests
2. `test_prep.py` - Preparation tests (depends on fetchers)
3. `test_docking.py` - Docking tests (depends on prep)
4. Other tests follow similar dependency chain

### Test Fixtures
```python
@pytest.fixture(scope="session")
def shared_fetchers():
    """Ensure fetchers run before prep"""
    from test_fetchers import test_fetchers_end_to_end
    test_fetchers_end_to_end(True)
    return True
```

## File Locations

### Input/Output Structure
- Inputs: `data/input/{proteins,ligands}/`
- Work: `data/work/`
- Outputs: `data/output/run-{timestamp}/`
- Logs: `data/output/run-{timestamp}/logs/`

### Configuration Files
- Main configs: `configs/*.yaml`
- Test configs: `configs/*_test.yaml`
- Environment: `environment.yml`