# OVM-PK: Physics-Aware Protein-Ligand Simulation Pipeline

A comprehensive molecular dynamics pipeline for protein-ligand systems with a focus on metalloproteins and physics-based scoring.

## 🚀 Features

* **Interactive CLI**: User-friendly command-line interface for easy analysis
* **Automated Structure Retrieval**: Fetches protein structures and ligand information automatically
* **Physics-Aware Docking**: Advanced scoring functions incorporating metal coordination, polarization, and solvation effects
* **End-to-End Pipeline**: From structure preparation to simulation and analysis
* **Modular Architecture**: Easily extensible with custom scorers and physics models
* **Reproducible**: Version-controlled configurations and containerized environments
* **MM/GBSA benchmarking**: Lightweight calculator with snapshot relaxation, ensemble sampling, and documentation for ΔG validation
* **High Performance**: GPU-accelerated simulations with OpenMM

## 📦 Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/ovm-pk.git
   cd ovm-pk
   ```

2. Set up the conda environment:
   ```bash
   mamba env create -f environment.yml
   mamba activate ovmpk
   pip install -e .
   ```

## 🏃‍♂️ Quick Start

### Interactive Mode (Recommended for New Users)

Run the interactive CLI to analyze protein-ligand interactions:

```bash
ovmpk run --interactive
```

Follow the prompts to:
1. Enter a protein name or UniProt ID (e.g., "CYP3A4" or "P08684")
2. Optionally enter a ligand name (e.g., "ketoconazole")
3. Let the pipeline handle the rest!

### Advanced Usage

For more control, you can also use a configuration file:

1. Prepare a configuration file (see `configs/template_ligand_enzyme.yaml`)
2. Run the pipeline:
   ```bash
   ovmpk run -c configs/your_config.yaml
   ```

### MM/GBSA Binding Energy Workflow

Once weight-scan results are generated, the CLI can compute MM/GBSA binding free energies directly:

```bash
# Baseline energies without additional relaxation
python -m ovmpk.cli.run_analysis mmgbsa \
  ovmpk_results/benchmarks/weight_scan_results.json \
  --out ovmpk_results/benchmarks/mmgbsa_results.json \
  --gbsa-model GBn2

# Relax complex/receptor/ligand snapshots, archive intermediates, then evaluate energies
python -m ovmpk.cli.run_analysis mmgbsa \
  ovmpk_results/benchmarks/weight_scan_results.json \
  --out ovmpk_results/benchmarks/mmgbsa_relaxed_results.json \
  --gbsa-model GBn2 \
  --minimize \
  --relax \
  --relax-output ovmpk_results/benchmarks/relaxed_snapshots \
  --relax-md-steps 1000 \
  --relax-random-seed 314159

# Launch an ensemble with MD sampling (reuses relaxed snapshots when provided)
python -m ovmpk.cli.run_analysis mmgbsa \
  ovmpk_results/benchmarks/weight_scan_results.json \
  --out ovmpk_results/benchmarks/mmgbsa_ensemble_results.json \
  --gbsa-model GBn2 \
  --minimize \
  --relax \
  --md-steps 5000 \
  --sample-interval 50 \
  --relax-random-seed 314159
```

Each JSON output captures complex/receptor/ligand energies (kcal/mol), ΔG_bind, and—when relaxation is enabled—per-state energy diagnostics plus paths to archived PDBs for reproducibility. See [`docs/pipeline_stages.md`](docs/pipeline_stages.md) for a deeper discussion of validation strategy, replica ensembles, and next-tier free-energy corrections.

## 🏗️ Project Structure

```
ovm-pk/
├── configs/              # Configuration files
├── data/                 # Input/output data (gitignored)
│   ├── input/            # Input structures and parameters
│   └── output/           # Pipeline outputs
├── docs/                 # Documentation
├── examples/             # Example scripts
├── forcefields/          # Force field parameters
├── scripts/              # Utility scripts
├── src/                  # Source code
│   └── ovmpk/            # Main package
│       ├── analysis/     # Analysis tools
│       ├── config/       # Configuration management
│       ├── docking/      # Docking engines and scorers
│       ├── fetchers/     # Data fetching utilities
│       ├── md/           # Molecular dynamics
│       ├── physics/      # Physics models
│       └── utils/        # Utility functions
└── tests/                # Test suite
```

## 📚 Documentation

For detailed documentation, please see the [docs](docs/) directory.

## 🤝 Contributing

Contributions are welcome! Please see our [Contributing Guidelines](CONTRIBUTING.md).

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📧 Contact

For questions or support, please open an issue or contact the maintainers.

## 🛠️ System Requirements

- **OS**: Linux (Ubuntu 22.04+ recommended, WSL2 supported)
- **CPU**: x86_64 with AVX2 support
- **GPU**: NVIDIA GPU with CUDA support (recommended) or CPU-only mode
- **Memory**: 16GB RAM minimum, 32GB+ recommended
- **Storage**: 50GB+ free space for simulations

### Dependencies

- **System**:
  ```bash
  sudo apt-get update
  sudo apt-get install -y build-essential cmake git curl wget \
      mesa-utils xauth x11-apps
  ```

- **NVIDIA (for GPU acceleration)**:
  - Install NVIDIA drivers and CUDA toolkit
  - Ensure `nvidia-smi` shows your GPU

### Python Environment

The project uses Mamba/Conda for dependency management. Install dependencies with:

```bash
mamba env create -f environment.yml
mamba activate ovmpk
pip install -e .
```

## 🔍 Testing

Run the test suite with:

```bash
pytest tests/ -v
```

## 📦 Data Management

Input and output data are stored in the `data/` directory, which is gitignored by default. The recommended structure is:

```
data/
├── input/
│   ├── proteins/     # PDB/mmCIF files
│   └── ligands/      # Ligand files (SDF, MOL2, PDBQT)
└── output/           # Pipeline outputs
    ├── 01_prep/      # Preprocessed structures
    ├── 02_docking/   # Docking results
    ├── 03_scoring/   # Physics-based scoring
    ├── 04_simulation/# MD trajectories
    └── 05_analysis/  # Analysis results
```
## 🚀 Quick Start with Example

1. **Prepare your environment**:
   ```bash
   # Clone the repository
   git clone https://github.com/yourusername/ovm-pk.git
   cd ovm-pk
   
   # Set up the environment
   mamba env create -f environment.yml
   mamba activate ovmpk
   pip install -e .
   ```

2. **Run an example simulation**:
   ```bash
   # Create input directories
   mkdir -p data/input/proteins data/input/ligands
   
   # Copy example files (replace with your own)
   cp examples/data/1abc.pdb data/input/proteins/
   cp examples/data/ligand.mol2 data/input/ligands/
   
   # Run the pipeline
   python -m ovmpk.pipeline_v2 configs/example_config.yaml -o data/output
   ```

3. **View results**:
   - Docking poses: `data/output/02_docking/`
   - Physics scores: `data/output/03_scoring/physics_scores.json`
   - Simulation trajectories: `data/output/04_simulation/`
   - Analysis results: `data/output/05_analysis/`

## 📊 Example Configuration

Here's a minimal configuration file example (`configs/example_config.yaml`):

```yaml
system:
  protein_pdb: data/input/proteins/1abc.pdb
  ligand_pdb: data/input/ligands/ligand.mol2
  metals:
    - element: FE
      position: [10.0, 20.0, 30.0]
      oxidation_state: 2
      coordination: 6

scoring:
  method: physics_based
  weights:
    w_vdw: 1.0
    w_elec: 0.6
    w_pol: 0.2
    w_coord: 1.5
    w_desolv: 0.5
    w_tors: 0.1

simulation:
  forcefield: amber14-all.xml
  water_model: tip3p.xml
  temperature: 300  # Kelvin
  pressure: 1.0     # bar
  nsteps: 1000000   # 2 ns with 2 fs timestep
```

## 🧪 Testing the Installation

To verify your installation, run the test suite:

```bash
# Run all tests
pytest tests/

# Run specific test module
pytest tests/test_physics_scoring.py -v

# Run with coverage report
pytest --cov=src/ovmpk tests/
```

## 🤝 Contributing

We welcome contributions! Please follow these steps:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add some amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📜 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Documentation

For detailed documentation, please see the [docs](docs/) directory.

## 📧 Contact

For questions or support, please open an issue on GitHub or contact the maintainers.

## 🚀 Next Steps

- Explore the [examples](examples/) directory for more usage examples
- Check out the [configuration guide](docs/configuration.md) for advanced options
- Read about [custom physics models](docs/custom_physics.md) to extend the functionality
- Join our [community discussions](https://github.com/yourusername/ovm-pk/discussions) to ask questions and share ideas
If you ever need to **rebuild**:

```bash
mamba deactivate || true
mamba env remove -n ovmpk || true
./scripts/setup_conda.sh
mamba activate ovmpk
pip install -e .
```

---

## Quick Start

```bash
# Install dependencies
mamba env create -f environment.yml
mamba activate ovm-pk

# Run the full pipeline
python scripts/run.py pipeline configs/default.yaml

# Or run individual stages
python scripts/run.py md prepare input.pdb -o system.xml
python scripts/run.py md run system.xml -o trajectory.dcd

# Get help
python scripts/run.py --help
python scripts/run.py pipeline --help
python scripts/run.py md --help
```

---

## Pipeline Overview

### 1) Fetch

* **Proteins** from RCSB (filters: method, resolution, required non-polymers e.g., HEM).
* **Ligands** by name/CID/SMILES (PubChem/ChEMBL resolution and sanity).

### 2) Prep

* **Protein**: light clean → **PDBFixer** (missing atoms + H @ pH) → PDB.
* **Ligand**: **dimorphite-dl** protonation at pH → SDF; PDBQT later via OpenBabel.

### 3) Dock (smina; GNINA optional)

* Receptor/ligand **PDBQT** via OpenBabel.
* **Boxing**:

  * Explicit center/size **or**
  * `autocenter_box_from_pdb(...)` to detect **heme Fe** and center the box.
* **Multi-seed** runs; outputs SDF poses + per-seed logs.

### 4) Pose selection

* Pick **best-scoring pose**, assert Fe–N proximity when relevant.

### 5) MD

* **Equilibration**: short NVT + NPT. **Production**: short run; DCD + `.chk`.
* **Resume**: continue from checkpoint. **Analysis**: basic integrity checks.

---

## Configs

See `configs/prod_test.yaml`. Key knobs (excerpt):

```yaml
docking:
  num_seeds: 5
  poses: 20
  exhaustiveness: 16
  seed: 42
  box:
    center: [-22.578, -26.236, -11.422]
    size: [12.0, 12.0, 12.0]
    auto_center: heme_fe
    fallback_center: [-15.846, -23.032, -11.293]
md:
  temperature: 310
  pressure: 1.0
  length_ns: 10
```

You can override ligand identity with `OVMPK_LIGAND_SMILES` for deterministic prep.

---

## Tests

Covers the full stack:

* **Fetchers** (RCSB + ligand resolver)
* **Prep** (protein & ligand)
* **Docking** (multi-seed smina; auto heme-Fe centering)
* **Pose selection** (Fe–N sanity)
* **MD** (prepare → NVT → NPT → prod → resume → analysis)

All tests are under `tests/` and are runnable individually or in suites.

---

## Visualization (PyMOL)

**Quick viewer** for the exact objects used in docking:

```pml
# Run inside PyMOL (File → Run or 'pymol -cq view.pml -- ...' style)
load data/work/docking/5VCC_fixed_ph7.4_cleaned_fixed.pdbqt, rec
load data/work/docking/ketoconazole_prepared_ph7.4.pdbqt, lig_docked  # if you saved PDBQT pose
# or the selected SDF pose:
# load data/output/pose_selection_test/selected_pose/ligand_best_pose_1_from_5VCC_fixed_ph7.4_cleaned_fixed_ketoconazole_prepared_ph7.4_s5573673_poses.sdf, lig_docked

hide everything
as surface, rec
set transparency, 0.25, rec
split_states lig_docked; disable lig_docked
as sticks, lig_docked_0001
util.cbaw lig_docked_0001
select pocket, byres (rec within 4 of lig_docked_0001)
show sticks, pocket
color tv_yellow, pocket
orient lig_docked_0001
```

> If your desktop has no X server, use `pymol -cq` to render images headlessly.

---

## GNINA (CNN) integration (optional)

* **Binary**: put `gnina` on `PATH` (official release or container).
* **Enable**: set `docking.engine: gnina` in a config, or keep `smina` and set `docking.rescore_with_gnina: true`.
* **Device**: `docking.gnina.device: auto|gpu|cpu` (auto prefers GPU).
* **Tests** (will skip if GNINA not found):

  ```bash
  pytest -q tests/test_gnina_presence.py tests/test_gnina_rescore.py
  ```

Artifacts will appear as `*_gnina*.{sdf,log}` under `data/work/docking/` with metrics alongside.

---

## Project Structure

```bash
ovm-pk/
├── configs/           # YAML configs for different pipeline stages
├── data/              # Input/Output data (gitignored)
├── docs/              # Documentation
├── examples/          # Example scripts and configurations
├── scripts/           # Pipeline scripts and tools
│   ├── tools/         # Specialized tools (e.g., autopocket)
│   └── utils/         # Utility scripts
├── src/               # Source code
│   └── ovmpk/         # Main package
│       ├── analysis/  # Analysis tools
│       ├── config/    # Configuration handling
│       ├── docking/   # Docking implementations
│       ├── fetch/     # Data fetching
│       ├── md/        # Molecular dynamics
│       ├── physics/   # Physics models and scoring
│       └── prep/      # Structure preparation
└── tests/             # Tests for each pipeline stage
```
tests/
  …                         # pytest suites across the pipeline
data/
  input/  work/  output/    # ephemeral artifacts (gitignored)
view.pml                    # handy PyMOL viewing macro
```

---

## 🧹 Housekeeping

* **Ignore rules**: large/binary artifacts, structures, logs are excluded by `.gitignore`.
* **Clean before re-runs**:

  ```bash
  rm -rf data/work/* data/output/* runlogs/*
  ```

---

## 🧭 Design Notes

* **Why smina?** Portable CPU search, multi-seed coverage.
* **Why OpenBabel PDBQT?** Stable PDB/SDF↔PDBQT conversions + Gasteiger charges.
* **Why OpenMM?** Modern GPU MD, straightforward XML & checkpoint I/O.
* **Why pytest?** Treat the whole pipeline like code: repeatable, observable, testable.

---

## 🆘 Troubleshooting

* **PyMOL can’t open PDBQT** → ensure you load with `load file.pdbqt, obj` (PyMOL understands PDBQT). If needed, `obabel in.pdbqt -O out.pdb` for viewing.
* **GPU not used** → check `./scripts/setup_wsl_gpu.sh` output and that OpenMM selects `CUDA` (see logs).
* **Headless** → use `pymol -cq …` or export images with `png out.png, ray=1`.

---

## 🙏 Acknowledgments

* **smina** (Vina-compatible) — Koes lab
* **OpenMM** — Eastman/Pande et al.
* **OpenBabel** — O’Boyle et al.
* **RDKit** — Landrum et al.
* **Open Force Field Toolkit** — OpenFF consortium
* **PDBFixer** / **MDTraj**
* **RCSB PDB APIs** for search & retrieval

---

## 📄 License

TBD. Until then, assume research/educational use only; check licenses of third-party tools listed above.

---

## 🔗 Changelog (high-level)

* Added one-shot bootstrap scripts (APT, WSL GPU check, mamba env, tool verifier)
* End-to-end pytest suite with rich UI and fail-fast
* Multi-seed docking; **auto-center box from heme Fe**
* Stepwise MD prep + short NVT/NPT/prod + resume + analysis
* `.gitignore` hardened for binary/trajectory outputs
* GNINA (optional) enablement & tests

