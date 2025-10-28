# ovm-pk

A pragmatic protein–ligand pipeline that goes **Fetch → Prep → Dock (smina) → Pose select → MD (OpenMM) → Resume/Analyze** with opinionated defaults, reproducible configs, and an end-to-end pytest harness.

---

## ✨ Highlights

- **End-to-end tests** (pytest) for every stage with rich progress + live logs
- **Multi-seed docking** via `smina` (Vina-compatible) with automatic box centering from heme Fe
- **Pose selection** (Fe–N sanity gate + best-energy pick)
- **MD prep/equil/prod** on GPU (OpenMM, TIP3P, PME) with restart checkpoints
- **Deterministic configs** (`configs/*.yaml`) + environment flags for quick smoke vs. longer runs
- **Repro hygiene**: run artifacts and structure files are `.gitignore`d by default

---

## 🧱 Requirements

- Linux (WSL2 Ubuntu 22.04 OK) with an NVIDIA GPU + CUDA drivers (for OpenMM CUDA)
- Mamba/Conda (recommended)  
- External tools on `PATH`:
  - **smina** (Koes lab, Vina-compatible)
  - **OpenBabel** (`obabel`) for PDB/SDF↔PDBQT conversions
- Python libs (installed via `environment.yml`):
  - OpenMM, MDTraj, RDKit, OpenFF Toolkit
  - rcsbsearchapi / httpx (structure fetch)
  - rich / pytest-rich (pretty test output)

> Tip: `mamba env create -f environment.yml && mamba activate ovmpk`

---

## 📦 Install (editable)

```bash
mamba activate ovmpk
pip install -e .
```

---

## 🚀 Quickstart

Run a compact, GPU-light pipeline using the production-like test config:

```bash
# E2E with pretty output, live logs, and fail-fast
mkdir -p runlogs && OVMPK_EQUIL_STEPS_NVT=1000 OVMPK_EQUIL_STEPS_NPT=1500 OVMPK_PROD_STEPS=5000 pytest --rich --rich-capture=runlogs/pytest-rich.txt   -vv --capture=tee-sys -rA --durations=0 --maxfail=1 --color=yes   -o log_cli=true -o log_cli_level=INFO   --log-file=runlogs/pytest-e2e.log --log-file-level=INFO   tests/test_fetchers.py   tests/test_prep.py   tests/test_docking.py   tests/test_pose_selection.py   tests/test_md_prepare.py   tests/test_md_equil.py   tests/test_md_prod.py   tests/test_md_resume.py   tests/test_md_analysis.py
```

**Reset workspace (safe):**
```bash
rm -rf data/work/* data/output/* runlogs/*
```

---

## 🗺️ Pipeline Overview

### 1) Fetch
- **Proteins** from RCSB (filters: method, resolution, required non-polymers e.g., HEM).
- **Ligands** by name/CID/SMILES (PubChem/ChEMBL resolution and sanity).

### 2) Prep
- **Protein**: lightweight chain/HETATM clean → PDBFixer (missing atoms + H @ pH) → PDB
- **Ligand**: dimorphite-dl pH handling (via prep script) → SDF; later to PDBQT

### 3) Dock (smina)
- Receptor/ligand **PDBQT** generated via OpenBabel
- **Boxing**:
  - Configured bounding box **or**
  - `autocenter_box_from_pdb(...)` that detects **heme Fe** and centers box there
- **Multi-seed** runs (configurable seeds, poses, exhaustiveness)
- Outputs: SDF poses + logs per seed

### 4) Pose selection
- Picks the **best-scoring pose** and asserts **Fe–N** proximity when relevant

### 5) MD
- **md_gpu_prepare.py** stages: `ligand → receptor → merge → forcefield → solvate → system → minimize → all`
- **OpenMM** CUDA platform, TIP3P water, 0.15 M ions, XML system cache
- **Equilibration**: short NVT+NPT
- **Production**: short run; DCD + log + checkpoint
- **Resume**: loads `.chk` and continues
- **Analysis**: checks trajectory frames, thermo windows, and artifact integrity

---

## 🔧 Configs

See `configs/prod_test.yaml` (used by tests). Key knobs:

```yaml
fetch:
  protein:
    target_species_id: 9606
    exp_method: X-RAY DIFFRACTION
    max_resolution: 3.0
    required_ligands: ["HEM"]
  ligand:
    identifier_type: name     # name | cid | smiles | inchikey

prep:
  protein:
    ph: 7.4
    run_pdbfixer: true
    output_suffix: "_fixed_ph7.4"
  ligand:
    ph: 7.4
    charge_method: gasteiger
    force_3d: true
    output_suffix: "_prepared_ph7.4"
  strip_heme: true
  strip_resnames: [GOL, SO4, EDO, PEG, TRS, ACT, PG4, PO4, CL, NA]

docking:
  structure_policy: prefer_apo
  blind_mode: false
  num_seeds: 5
  poses: 20
  exhaustiveness: 16
  seed: 42
  box:
    # Either provide center/size …
    center: [-22.578, -26.236, -11.422]
    size: [12.0, 12.0, 12.0]
    # … or enable auto-centering (heme Fe); explicit center wins if both are present
    auto_center: heme_fe          # heme_fe | none
    fallback_center: [-15.846, -23.032, -11.293]

md:
  temperature: 310
  pressure: 1.0
  length_ns: 10
  water_model: tip3p
  padding: 10.0
  ionic: 0.15

runtime:
  serialize_gpu: true
  cpu_threads: 16

reporting:
  outdir: results/reports_prod
  run_root: results/runs_prod
```

You can also pass `OVMPK_LIGAND_SMILES` to force a trusted template SMILES (useful when a pose SDF carries odd radicals from RDKit perception).

---

## 🧪 Tests

- Rich, fail-fast, and timed:
  - **Fetchers**: structure + ligand retrieval
  - **Prep**: protein/ligand prep
  - **Docking**: multi-seed smina with **auto heme-Fe centering**
  - **Pose selection**: Fe–N distance sanity
  - **MD**: stepwise prep, short NVT/NPT, production, resume, and analysis (DCD frames/thermo)
- Minimal `conftest.py` wires **run root**, **output dirs**, **platform discovery**, and live logging.

---

## 📁 Repo Layout

```
src/ovmpk/
  docking/smina_wrapper.py     # obabel↔pdbqt, multi-seed docking, auto-center
  utils/logging.py             # unified logging
  utils/smiles_fetch.py        # trusted SMILES resolver (PubChem/ChEMBL)
scripts/
  md_gpu_prepare.py            # staged MD prep & system build (OpenMM)
configs/
  prod_test.yaml               # production-like test config
tests/
  test_fetchers.py
  test_prep.py
  test_docking.py
  test_pose_selection.py
  test_md_prepare.py
  test_md_equil.py
  test_md_prod.py
  test_md_resume.py
  test_md_analysis.py
data/
  input/    # kept empty via .gitkeep
  work/     # ephemeral (gitignored)
  output/   # ephemeral (gitignored)
```

---

## 🧹 Housekeeping

- **Ignore rules**: large/binary artifacts, structures, logs are excluded by `.gitignore`
- **Clean before re-runs**:
  ```bash
  rm -rf data/work/* data/output/* runlogs/*
  ```

---

## 🧭 Design Notes

- **Why smina (Vina-compatible)?** Portable, fast CPU search; multi-seed for coverage.
- **Why OpenBabel for PDBQT?** Stable PDB/SDF↔PDBQT conversions; gasteiger charges.
- **Why OpenMM?** Modern GPU MD with Pythonic APIs, easy XML system I/O and checkpointing.
- **Why OpenFF?** Parameterization robustness for drug-like ligands; integrates with RDKit pipelines.
- **Why pytest?** Treat the whole pipeline like code: repeatable, observable, testable.

---

## 🙏 Acknowledgments

- **smina** (Vina-compatible) — Koes lab  
- **OpenMM** — Eastman/Pande et al.  
- **OpenBabel** — O’Boyle et al.  
- **RDKit** — Landrum et al.  
- **Open Force Field Toolkit** — OpenFF consortium  
- **PDBFixer** / **MDTraj**  
- **RCSB PDB APIs** for search and structure retrieval

---

## 📄 License

TBD. Until then, assume research/educational use only; check licenses of third-party tools listed above.

---

## 🔗 Changelog (high-level)

- Added end-to-end pytest suite with rich UI and fail-fast
- Multi-seed docking; **auto-center box from heme Fe**
- Stepwise MD prep + short NVT/NPT/prod + resume + analysis
- `.gitignore` hardened for binary/trajectory outputs

### GNINA (CNN) integration

This repo can optionally use **GNINA** for docking or rescoring with CNN models.

- **Install**: place `gnina` on your `PATH` (preferred: official binary/Docker).  
- **Select engine**: set `docking.engine: gnina` in a config, or keep `smina` and enable `docking.rescore_with_gnina: true`.  
- **Device**: `docking.gnina.device: auto|gpu|cpu`. Auto prefers GPU if available.  
- **Tests**: `pytest -q tests/test_gnina_presence.py tests/test_gnina_rescore.py` (skips if GNINA not found).  

Artifacts land under `data/work/docking/*_gnina*.{sdf,log}` and per-run metrics JSON in `results/...` when enabled.
