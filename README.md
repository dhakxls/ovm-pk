# ovm-pk

A pragmatic protein–ligand pipeline (docking → MD → free energy).

## Quickstart

```bash
# (Recommended) create env via conda/mamba
mamba env create -f environment.yml
mamba activate ovmpk

# Install package in editable mode (optional)
pip install -e .

# Smoke run (uses dry-run stubs; no internet/tools required)
OVM_DRY_RUN=1 ovmpk run --config configs/fast_test.yaml --gene CYP3A4 --ligand ketoconazole
```

See `configs/*.yaml` for knobs (structure_policy, blind_mode, seeds, etc.).
