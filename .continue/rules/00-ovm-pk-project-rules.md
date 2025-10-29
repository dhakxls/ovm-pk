# OVM-PK Project Rules

## Project Overview

OVM-PK is a physics-aware docking and validation pipeline for protein-ligand interactions, focusing on:
- Experimental affinity (Ki) and ΔG replication
- Pluggable docking physics modules
- Head-to-head comparison against established baselines

## Project Architecture

The project follows a modular structure:
- `src/ovmpk/` - Main package implementation
  - `docking/` - Docking engines and scoring modules
  - `fetchers/` - RCSB and ligand fetching utilities
  - `prep/` - Protein and ligand preparation tools
  - `utils/` - Common utilities and logging
- `scripts/` - Environment setup and MD runners
- `configs/` - YAML configuration files (default/fast_test/prod_test)
- `data/` - Data directory structure
  - `input/` - Input files for proteins and ligands
  - `work/` - Working directory for intermediate files
  - `output/` - Results and analysis outputs
- `tests/` - Test suite following dependency order:
  1. fetchers → 
  2. prep → 
  3. docking → 
  4. pose_selection → 
  5. md_* tests → 
  6. others

## Coding Standards

### Python Standards
- Use Python 3.10+ features
- Type hints required for all new code
- Follow existing import order: built-ins → third-party → local
- Doc strings required for all public functions and classes

### Testing Standards
- All new features must have corresponding tests
- Tests must be placed in appropriate dependency order
- Use pytest fixtures for shared resources
- Test files should match the module they test

### File Naming Conventions
- Module files: lowercase with underscores (e.g., `protein_prep.py`)
- Test files: prefix with `test_` (e.g., `test_protein_prep.py`)
- Config files: descriptive YAML names (e.g., `prod_test.yaml`)

## Data Management

### Directory Structure
- Input data goes in `data/input/{proteins,ligands}`
- Working files in `data/work/`
- Test outputs in `data/output/run-{timestamp}`
- Cache files in `data/cache/`

### File Formats
- Proteins: PDB, PDBQT formats
- Ligands: SDF, MOL2 formats
- Trajectories: DCD format
- Config: YAML format

## Configuration Standards

### YAML Configuration
- All runs must specify a config file
- Use `configs/prod_test.yaml` as template
- Required sections:
  - protein configuration
  - ligand parameters
  - docking settings
  - MD parameters

## Testing Workflow

### Test Organization
1. Fetchers tests must run first
2. Prep tests depend on fetcher results
3. Docking tests depend on prep
4. All MD tests follow

### Test Output
- Logs go to `data/output/run-{timestamp}/logs`
- Each test run gets unique timestamp
- Clean up old test outputs regularly

## Common Operations

### Environment Setup
- Use Mamba/Conda environment from `environment.yml`
- Run `scripts/create_or_update_env.sh` for setup
- Activate environment before any operation

### Running Tests
- Full suite: `pytest tests/`
- Single module: `pytest tests/test_module.py`
- With output: `pytest -v tests/`

### Code Generation Rules
- New modules go in appropriate `src/ovmpk/` subdirectory
- Match existing file structure
- Include appropriate imports and type hints
- Add corresponding test file

## Refactoring Guidelines

### Code Organization
- Keep modules focused and single-purpose
- Use appropriate subdirectory in `src/ovmpk/`
- Update imports when moving files
- Update tests to match new structure

### Dependency Management
- Add new dependencies to `environment.yml`
- Run `scripts/freeze_env.sh` after changes
- Update `runlogs/env.lock.yml` for reproducibility

## File Operations

### Temporary Files
- Use `data/work/` for intermediate files
- Clean up after operations complete
- Don't commit temporary files

### Output Management
- Use timestamps in output directories
- Clean up old output directories
- Keep latest successful run

## Version Control

### Git Guidelines
- Don't commit binary files
- Don't commit temporary files
- Follow existing .gitignore rules
- Keep commits focused and well-described

### Protected Files
- Don't modify `environment.yml` without approval
- Don't commit `data/` contents except .gitkeep
- Don't commit cache files or outputs