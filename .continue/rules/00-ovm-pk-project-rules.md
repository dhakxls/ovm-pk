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
- `tests/` - Test suite following dependency order

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

[Rest of the content continues...]
