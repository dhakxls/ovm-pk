# OVM-PK: Physics-Aware Docking & Validation Pipeline

## Overview
This document details the 5-stage pipeline for physics-aware docking and validation of protein-ligand systems, with special focus on metalloproteins. The pipeline is designed to be modular, with pluggable scoring functions and validation metrics to enable accurate prediction of binding affinities (ΔG) and comparison with experimental data.

---

## Stage 1: Data Fetching & Preparation
**Purpose**: Retrieve, prepare, and validate initial structures and experimental data

### Inputs
- PDB ID (e.g., "1TQN")
- Ligand ID (e.g., "KET" or SMILES string)
- Experimental binding data (Ki/Kd/ΔG°)
- Configuration (pH, temperature, force fields)

### Process
1. **Data Retrieval**
   - Fetch PDB/mmCIF from RCSB with biological assembly
   - Retrieve ligand from PubChem/ChEMBL or generate from SMILES
   - Download experimental binding data from public databases

2. **Structure Preparation**
   - Add missing atoms and predict protonation states (PROPKA)
   - Generate ligand 3D coordinates and protonation states
   - Handle metal coordination and cofactors

3. **Validation & Quality Control**
   - Validate electron density fit (if available)
   - Check stereochemistry and geometry
   - Verify experimental data consistency

### Outputs
- `protein_prepared.pdb`: Prepared protein structure
- `ligand_prepared.sdf`: Prepared ligand structure
- `experimental_data.json`: Binding affinity and conditions
- `validation_report.html`: Validation results

---

## Stage 2: System Setup & Parameterization
**Purpose**: Prepare the simulation system with appropriate force fields and parameters

### Inputs
- `protein_prepared.pdb` (from Stage 1)
- `ligand_prepared.sdf` (from Stage 1)
- Configuration (force fields, water model, ion concentration)

### Process
1. **System Assembly**
   - Add explicit solvent (TIP3P/OPC)
   - Add counterions for charge neutralization
   - Set up periodic boundary conditions

2. **Force Field Assignment**
   - Protein: AMBER/CHARMM/OPLS
   - Ligand: GAFF/CGenFF with AM1-BCC charges
   - Water: TIP3P/OPC
   - Special parameters for metals and cofactors

3. **Parameter Validation**
   - Check energy terms
   - Validate parameter assignments
   - Verify system neutrality

### Outputs
- `system.prmtop`: Topology file
- `system.inpcrd`: Coordinates
- `parameters/`: Force field files
- `system_validation.json`: Validation metrics

---

## Stage 3: Physics-Aware Docking
**Purpose**: Generate and score protein-ligand poses using physics-based scoring

### Inputs
- `system.prmtop` (from Stage 2)
- `system.inpcrd` (from Stage 2)
- Docking configuration (box size, exhaustiveness, etc.)

### Process
1. **Pose Generation**
   - Define binding site (auto-detected or user-defined)
   - Generate initial poses using smina/Vina
   - Optional: Generate poses using alternative methods

2. **Physics-Based Scoring**
   - Metal coordination scoring
   - Polarization effects
   - Solvation/desolvation
   - Torsion strain
   - Van der Waals and electrostatics

3. **Pose Refinement**
   - Energy minimization
   - Short MD relaxation
   - Rescoring with full physics model

### Outputs
- `docking_results.csv`: Scored poses
- `top_poses/`: Top-scoring poses (PDB format)
- `scoring_breakdown.json`: Detailed energy terms

---

## Stage 4: Binding Affinity Prediction
**Purpose**: Calculate binding free energies and compare with experimental data

### Inputs
- Top-scoring poses (from Stage 3)
- Experimental binding data (from Stage 1)
- Simulation parameters

### Weight Scans & ΔG Calibration
1. **Physics Weight Scan**
   - Configure per-preset weights and damping overrides in `ovmpk_results/single_pose_config.yml` under `scoring.weight_scan`.
   - Run the scan via `python -m ovmpk.cli.run_analysis scan-weights --config <config.yml> --output <results.json>` to generate physics totals for each preset.

2. **ChEMBL Benchmarking**
   - Fetch activities with `python -m ovmpk.cli.run_analysis chembl CHEMBL262 --out ovmpk_results/benchmarks/chembl_CHEMBL262.json`.
   - Annotate ΔG values using `python -m ovmpk.cli.run_analysis convert-dg <activities.json> --out <dg.json>` (defaults to 298.15 K).

3. **Linear Calibration**
   - Fit slope/intercept mapping physics totals to experimental ΔG values (e.g., notebook or script using the scan output and ChEMBL ΔG set).
   - Apply the calibration with `python -m ovmpk.cli.run_analysis map-dg --input <results.json> --out <mapped.json> --slope <m> --intercept <b>`.

4. **Analysis**
   - Review `mapped.json` for each preset’s `mapped_delta_g` alongside energy breakdowns.
   - Compare mapped values to the ΔG distribution to assess scaling quality before propagating into downstream reporting.

### Process
1. **Free Energy Calculation**
   - MM/PBSA or MM/GBSA calculations
   - Alchemical transformation (if applicable)
   - Entropy estimation

2. **Experimental Comparison**
   - Convert between Ki/Kd/ΔG as needed
   - Temperature correction
   - Statistical analysis

3. **Validation**
   - Compare predicted vs experimental ΔG
   - Calculate error metrics (RMSE, MAE, R²)
   - Generate correlation plots

### Outputs
- `binding_energies.csv`: Predicted and experimental values
- `validation_metrics.json`: Statistical analysis
- `correlation_plots/`: Visualizations

---

## Stage 5: Analysis & Reporting
**Purpose**: Generate comprehensive reports and analysis of results

### Inputs
- All previous stage outputs
- Configuration for analysis

### Process
1. **Performance Analysis**
   - Docking success rates
   - Scoring function performance
   - Computational efficiency

2. **Visualization**
   - 3D visualization of binding poses
   - Energy component analysis
   - Structural features vs. binding

3. **Report Generation**
   - Summary of inputs and results
   - Statistical analysis
   - Quality metrics and validation

### Outputs
- `final_report.html`: Interactive report
- `results_summary.json`: Key metrics
- `analysis_plots/`: Visualizations

---

## Integration with OVM-PK PRD v2

### Key Features Implemented
- **Pluggable Scoring System**: Modular architecture for different scoring functions
- **Physics-Aware Docking**: Integration of metal coordination, polarization, and solvation
- **Experimental Validation**: Direct comparison with experimental binding data
- **Reproducible Workflow**: Containerized environment and version control

### Configuration Example
```yaml
# Example configuration for physics-aware docking
docking:
  engine: smina
  scorer: new_physics
  box:
    auto_center: heme_fe
    size: [20, 20, 20]  # Ångströms
  physics:
    metal_coordination: true
    polarization: true
    solvation: true
    weights:
      vdw: 1.0
      elec: 0.6
      pol: 0.2
      coord: 1.5
      desolv: 0.5
      tors: 0.1
```

### Running the Pipeline
```bash
# Run the complete pipeline
python -m ovmpk.pipeline_v2 --config configs/example.yaml

# Run specific stages
python -m ovmpk.pipeline_v2 --stages prepare,dock,score,analyze
```

---

## Stage 7: Analysis
**Purpose**: Extract meaningful data from simulations

### Inputs
- `trajectory.dcd`
- `system_metal.prmtop`
- Analysis parameters

### Process
1. **Stability Analysis**
   - RMSD/RMSF
   - Secondary structure
   - Hydrogen bonds

2. **Binding Analysis**
   - Ligand interactions
   - Binding free energy
   - Residence times

### MM/GBSA quickstart

The `run_analysis` CLI now exposes an MM/GBSA helper for single-snapshot binding energies that reuses the force-field assets bundled with the project (RESP-heme + GAFF ligand parameters).

```bash
# Evaluate the weight-scan poses in place (no additional minimization)
python -m ovmpk.cli.run_analysis mmgbsa \
  ovmpk_results/benchmarks/weight_scan_results.json \
  --out ovmpk_results/benchmarks/mmgbsa_results.json \
  --gbsa-model GBn2

# Optionally minimize complex/receptor/ligand before scoring
python -m ovmpk.cli.run_analysis mmgbsa \
  ovmpk_results/benchmarks/weight_scan_results.json \
  --out ovmpk_results/benchmarks/mmgbsa_minimized_results.json \
  --gbsa-model GBn2 \
  --minimize

# Relax snapshots prior to MM/GBSA and archive intermediates
python -m ovmpk.cli.run_analysis mmgbsa \
  ovmpk_results/benchmarks/weight_scan_results.json \
  --out ovmpk_results/benchmarks/mmgbsa_relaxed_results.json \
  --gbsa-model GBn2 \
  --minimize \
  --relax \
  --relax-output ovmpk_results/benchmarks/relaxed_snapshots \
  --relax-md-steps 1000 \
  --relax-random-seed 314159
```

Outputs (`mmgbsa_results.json`, `mmgbsa_minimized_results.json`, `mmgbsa_relaxed_results.json`) provide complex/receptor/ligand energies and ΔG_bind in kcal/mol; relaxed runs also embed per-state energies and paths to archived PDBs under `relaxation`. With the current CYP3A4 pose the raw snapshot yields ΔG_bind ≈ −73.8 kcal/mol across all presets, while minimization flips the sign (≈ +30.7 kcal/mol), underscoring the sensitivity of single-frame MM/GBSA to geometry preparation. Relax-and-archive workflows ensure that follow-up ensemble calculations reuse consistent structures and enable reproducibility of preprocessing choices.

For validation against experiment, compare these ΔG_bind values to the ChEMBL activity panel converted to free energies (`chembl_CHEMBL262_dg.json`). The present snapshot differs from the ChEMBL mean (~−8.8 kcal/mol) by ≈ −65 kcal/mol, motivating further refinement of the receptor/ligand preparation.

**Recommended next-tier validation roadmap**

1. **Replica MM/GBSA ensembles** – Launch ≥3 Langevin replicas with distinct seeds (pass `--md-steps`, `--sample-interval`, and `--relax-random-seed`) to quantify sampling uncertainty and report per-replica ΔG summaries.
2. **Entropy corrections** – Apply normal-mode or quasi-harmonic corrections on relaxed snapshots to estimate ΔS contributions, delivering ΔG = ΔH − TΔS for direct comparison with experiment.
3. **Absolute binding free energy (ABFE)** – Transition to alchemical double-decoupling using the existing `FreeEnergyCalculator` skeleton once MM/GBSA trends stabilize; use relaxed snapshots as initial states and reuse the archived PDBs to guarantee consistent topology preparation.
4. **Constrained minimization/heating** – When large structural drifts are observed post-relaxation, adopt restrained minimization followed by short heating to 300 K before production to preserve key contacts.

### Outputs
- `analysis/`: Plots and data
- `binding_energy.txt`: ΔG binding
- `interactions.csv`: Contact frequencies

---

## Stage 8: Validation & Reporting
**Purpose**: Ensure simulation quality and document results

### Inputs
- All previous outputs
- Validation criteria

### Process
1. **Quality Checks**
   - Energy convergence
   - Structural stability
   - Artifact detection

2. **Report Generation**
   - Simulation parameters
   - Validation metrics
   - Analysis results

### Outputs
- `validation_report.pdf`
- `simulation_metadata.json`
- `archive/`: Complete workflow snapshot

---

## Configuration Reference

### Example YAML Configuration
```yaml
pipeline:
  structure_prep:
    ph: 7.4
    keep_water: false
    
  system_setup:
    water_model: "tip3p"
    padding: 1.0  # nm
    ion_concentration: 0.15  # M
    
  simulation:
    minimization_steps: 5000
    equilibration_time: 0.3  # ns
    production_time: 100  # ns
    
  analysis:
    skip_frames: 100  # frames to skip for equilibration
    stride: 10  # frame stride for analysis
```

### Environment Variables
```bash
# Required
EXPERIMENT_NAME="CYP3A4_Ketoconazole"
FORCE_FIELD="amber14"
WATER_MODEL="tip3p"

# Optional
NUM_CPUS=8
GPU_ENABLED=true
MEMORY_GB=32
```

This document provides complete context for any AI assistant to understand and help implement any part of the pipeline. Each stage is modular and can be developed/tested independently while maintaining compatibility with the overall workflow.
