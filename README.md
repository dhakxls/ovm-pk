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

