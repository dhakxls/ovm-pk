# ovm-pk

A pragmatic protein–ligand pipeline that goes **Fetch → Prep → Dock (smina/GNINA opt-in) → Pose select → MD (OpenMM) → Resume/Analyze** with opinionated defaults, reproducible configs, and an end-to-end pytest harness.

---

## ✨ Highlights

* **End-to-end tests** (pytest) for every stage with rich progress + live logs
* **Multi-seed docking** via `smina` (Vina-compatible) with automatic box centering (heme Fe)
* **Pose selection** (Fe–N sanity gate + best-energy pick)
* **MD prep/equil/prod** on GPU (OpenMM, TIP3P, PME) with restart checkpoints
* **Deterministic configs** (`configs/*.yaml`) + env flags for smoke vs. longer runs
* **Repro hygiene**: run artifacts and structure files are `.gitignore`’d by default

---

## 🧱 Requirements (split by layer)

**APT (system)**

* Ubuntu 22.04+ (WSL2 OK), basic build & GUI bits: `build-essential`, `cmake`, `git`, `curl`, `wget`, `mesa-utils`, `xauth`, `x11-apps` (for X), optional `pymol` (viewer)
* NVIDIA driver (Windows host) + **CUDA toolkit (optional)** on WSL for GPU OpenMM

**Mamba/Conda (project env)**

* Python libs & CLI tools pinned via `environment.yml`:

  * `openmm`, `mdtraj`, `mdanalysis`, `rdkit`, `openbabel`, `pdbfixer`, `propka`, `pdb2pqr`, `smina`, `pytest`, etc.
* Pip extras: `dimorphite-dl`, `typer`, `pyyaml`, `requests`, `rcsb-api`

> TL;DR policy: **OS things with APT** (drivers, GL/X, system libs), **science stack with mamba**.

---

## 🧰 One-shot bootstrap (fresh WSL or new machine)

From repo root:

```bash
# 0) Clone and enter
git clone <your-fork-url>.git ovm-pk && cd ovm-pk

# 1) System deps (APT)
./scripts/setup_apt.sh       # installs build tools, X basics, optional PyMOL

# 2) WSL GPU (optional helper)
./scripts/setup_wsl_gpu.sh   # verifies NVIDIA & CUDA visibility inside WSL

# 3) Project env (mamba)
./scripts/setup_conda.sh     # creates/updates 'ovmpk' env from environment.yml

# 4) Install package (editable) & sanity checks
mamba activate ovmpk
pip install -e .
./scripts/verify_tools.sh    # asserts smina, obabel, openmm import, etc.
```

If you ever need to **rebuild**:

```bash
mamba deactivate || true
mamba env remove -n ovmpk || true
./scripts/setup_conda.sh
mamba activate ovmpk
pip install -e .
```

---

## 🚀 Quickstart (compact e2e run)

Use the production-like test config with shorter MD steps for speed:

```bash
mkdir -p runlogs
OVMPK_EQUIL_STEPS_NVT=1000 OVMPK_EQUIL_STEPS_NPT=1500 OVMPK_PROD_STEPS=5000 \
pytest --rich --rich-capture=runlogs/pytest-rich.txt \
  -vv --capture=tee-sys -rA --durations=0 --maxfail=1 --color=yes \
  -o log_cli=true -o log_cli_level=INFO \
  --log-file=runlogs/pytest-e2e.log --log-file-level=INFO \
  tests/test_fetchers.py \
  tests/test_prep.py \
  tests/test_docking.py \
  tests/test_pose_selection.py \
  tests/test_md_prepare.py \
  tests/test_md_equil.py \
  tests/test_md_prod.py \
  tests/test_md_resume.py \
  tests/test_md_analysis.py
```

**Reset workspace (safe):**

```bash
rm -rf data/work/* data/output/* runlogs/*
```

---

## 🗺️ Pipeline Overview

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

* `scripts/md_gpu_prepare.py` staged build → XML system cache.
* **Equilibration**: short NVT + NPT. **Production**: short run; DCD + `.chk`.
* **Resume**: continue from checkpoint. **Analysis**: basic integrity checks.

---

## 🔧 Configs

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

## 🧪 Tests

Covers the full stack:

* **Fetchers** (RCSB + ligand resolver)
* **Prep** (protein & ligand)
* **Docking** (multi-seed smina; auto heme-Fe centering)
* **Pose selection** (Fe–N sanity)
* **MD** (prepare → NVT → NPT → prod → resume → analysis)

All tests are under `tests/` and are runnable individually or in suites.

---

## 🖼️ Visualization (PyMOL)

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

## 🤖 GNINA (CNN) integration (optional)

* **Binary**: put `gnina` on `PATH` (official release or container).
* **Enable**: set `docking.engine: gnina` in a config, or keep `smina` and set `docking.rescore_with_gnina: true`.
* **Device**: `docking.gnina.device: auto|gpu|cpu` (auto prefers GPU).
* **Tests** (will skip if GNINA not found):

  ```bash
  pytest -q tests/test_gnina_presence.py tests/test_gnina_rescore.py
  ```

Artifacts will appear as `*_gnina*.{sdf,log}` under `data/work/docking/` with metrics alongside.

---

## 📁 Repo Layout

```
src/ovmpk/
  docking/…                 # engines, wrappers, selectors, scorers
  prep/…                    # protein & ligand
  fetchers/…                # PDB & ligand resolvers
  utils/…                   # logging, SMILES helpers, toolchain
scripts/
  setup_apt.sh              # system deps (APT)
  setup_conda.sh            # mamba env create/update from environment.yml
  setup_wsl_gpu.sh          # optional CUDA/GPU sanity for WSL
  verify_tools.sh           # asserts smina/obabel/openmm availability
  md_gpu_prepare.py         # staged MD build
configs/
  default.yaml | fast_test.yaml | prod_test.yaml
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

