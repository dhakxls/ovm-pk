# OVM-PK MVP: Product Specification — **Physics-Aware Docking & Validation**

> **Goal:** Extend the current OVM-PK pipeline so we can (1) **replicate experimental affinity (Ki) and ΔG** for benchmark protein–ligand pairs end-to-end, and (2) **invent and plug in a new docking physics** module, then **compare it head-to-head** against established baselines (e.g., smina/Vina scoring + short OpenMM relaxations) using the same pipeline, data, and metrics.

---

## 1) Why this matters (Vision & Outcomes)

* **Scientific validity:** Docking scores are often only *rank-correlated* with binding affinity. We want a pipeline that **closes the loop to experiment** by mapping predicted energies to **ΔG°** and **Ki** (or pKi) so we can quantify accuracy on an absolute scale, not just enrichment.
* **Modularity for innovation:** By **factoring scoring/physics into a pluggable module**, we can iterate on new force/energy ideas (e.g., metal coordination, polarization, explicit waters, SE(3) learned terms) **without rewriting the pipeline**.
* **Decision support:** For drug discovery and translational users, better physics = **better ranking and triage**, fewer false positives, and more reliable prioritization for synthesis/assays.

**Primary success criteria**

1. **Replication:** For a curated benchmark set, pipeline delivers **|ΔG_pred − ΔG_exp| ≤ 1.5–2.0 kcal/mol median** *or* **Pearson r ≥ 0.5** between predicted ΔG and experimental ΔG° (baseline target).
2. **Novel physics improvement:** The **new docking physics** achieves **statistically significant better correlation / lower error** vs the baseline on the same set, with **reproducible** configs and seeds.
3. **Operational excellence:** One-command **E2E run + rich logs + artifacts**, green test suite (pytest), reproducible environments, and pluggable scorers/engines.

---

## 2) Scope (MVP)

### Included

* **A. Experimental replication loop**

  * Curate N protein–ligand pairs with **Ki / Kd / ΔG°** (converted to a common reference; temperature-aware).
  * Extend pipeline to compute **ΔG_pred** from docking/relaxation outputs and **compare** to ΔG_exp with standardized metrics.
* **B. Physics plug-in system**

  * Define a stable interface for **Pose Generation** and **Energy/Score Evaluation**.
  * Implement baselines (current **smina/Vina** + optional short **OpenMM** relaxation/scoring).
  * Add **NewPhysicsScorer** (first candidate below) and A/B harness to evaluate it.
* **C. CI-style benchmarks**

  * Pytest-driven **calibration tests**, **regression tests**, and **performance tests** to keep results stable while we iterate.

### Deferred / Out-of-scope (for MVP)

* Full FEP/TI pipelines; long timescale sampling; ensemble receptors; automated multi-tautomer/microstate enumeration; full solvent free-energy methods. (We’ll architect hooks to add these later.)

---

## 3) Fundamentals (Physics, Chemistry, Biology) we’ll apply

* **Thermodynamics:** Standard binding free energy
  [
  \Delta G^\circ = RT \ln!\left(\frac{K_d}{C^\circ}\right) = -RT \ln!\left( K_a,C^\circ \right)
  ]
  with ( C^\circ = 1,\text{M} ), ( R ) gas constant, ( T ) in Kelvin.
  **At 298 K**: ( \Delta G^\circ,[\text{kcal/mol}] \approx 1.364,\log_{10}(K_d/\text{M}) ).
* **Ki vs Kd:** For competitive inhibitors under standard conditions, **Ki ≈ Kd** (careful with assay conditions; we’ll store temperature, ionic strength).
* **Scoring physics (baseline):** Steric/VdW (LJ), electrostatics (Coulomb), simplistic desolvation terms, empirical hydrogen bonds, rotatable bond penalties.
* **Where new physics can help:**
  Metal coordination (e.g., **heme Fe–N/O**), **polarization/induction**, **anisotropic VdW**, **pH-dependent protonation**, **explicit/bridging waters**, **side-chain flexibility**, and **microstate-aware scoring**.

---

## 4) Architecture & Interfaces

### 4.1 High-level flow (existing, proven)

Fetch (RCSB/ligand) → Prep (PDBFixer, Dimorphite/OpenBabel) → **Dock (smina)** → Pose selection (heme-aware) → MD prep (OpenFF/SMIRNOFF) → **Short OpenMM relax/eval** → Analysis (DCD/logs) → Reporting

### 4.2 New plugin seam

We introduce **two extension points**:

```text
ovmpk.docking.engines.BaseEngine
  - generate_poses(receptor: Receptor, ligand: Ligand, box: Box, n_poses:int) -> List[Pose]

ovmpk.docking.scorers.BaseScorer
  - score_pose(receptor: Receptor, pose: Pose, context: ScoreContext) -> EnergyBreakdown
  - score_batch(receptor, poses, context) -> List[EnergyBreakdown]
```

* **Engine**: how poses are proposed (default = smina/Vina).
* **Scorer**: how a pose is **physically evaluated** (default = baseline Vina term; optional OpenMM re-scoring).

**Config (YAML) additions**

```yaml
docking:
  engine: smina             # or newton, equivdock, etc. (future)
  scorer: baseline_vina     # or new_physics, openmm_rescore, ml_rescore
  box:
    auto_center: heme_fe
    fallback_center: [-22.58, -26.24, -11.42]
    size: [12.0, 12.0, 12.0]
scoring:
  new_physics:
    metal:
      enable: true
      metals: [FE]         # atom names/elements to recognize
      preferred_bonds: [{donor: N, ideal: 2.1, k: 50.0}, {donor: O, ideal: 2.0, k: 40.0}]
      angular_penalty: {ideal: 180, k: 1.0}    # linear Fe–N–X, tweak per chemistry
    polarization:
      enable: true
      model: point_dipoles # MVP: simple inducible dipoles / screened
      alpha_scale: 1.0
    desolvation:
      model: gb_like       # fast GB-style term for MVP
    waters:
      enable: optional     # stub; reserve hook for “water-place” term
    weights:               # tunable in calibration
      w_vdw: 1.0
      w_elec: 0.6
      w_pol: 0.2
      w_coord: 1.5
      w_desolv: 0.5
      w_tors: 0.1
```

---

## 5) The **New Docking Physics** (MVP candidate)

**Concept:** A **metal-aware, polarization-augmented, fast solvent** scoring term that can be evaluated per pose in milliseconds, designed for **heme proteins** and similar metalloproteins.

**Energy form (illustrative, tunable):**
[
E_{\text{total}} = w_{\text{vdw}}E_{\text{LJ}} + w_{\text{elec}}E_{\text{Coulomb}}^{\text{screened}}

* w_{\text{pol}}E_{\text{pol}} + w_{\text{coord}}E_{\text{coord}}^{\text{metal}}
* w_{\text{desolv}}E_{\text{GB-like}} + w_{\text{tors}}E_{\text{lig_int}}
  ]

- **(E_{\text{coord}}^{\text{metal}})**: harmonic distance + angular preference for Fe–donor(s), optional multi-dentate penalties/rewards.
- **(E_{\text{pol}})**: simple induced-dipole response (MVP: analytic screened variant or pre-tabulated atom-type polarizabilities).
- **Screened electrostatics + GB-like desolvation**: fast approximation to capture salt/polar environment.
- **Ligand internal penalty**: torsion strain to avoid over-strained poses.

**Implementation strategy (fast & modular):**

* Start as a **pure Python/C++ vectorized scorer** operating on pose coordinates (NumPy/Numba or tiny C++ ext); later add **OpenMM Custom*Force** re-scoring for GPU speed if needed.
* Use **existing prepared coordinates/topology** from the pipeline (no re-param needed for MVP).
* Wire into `BaseScorer` + config flag `docking.scorer=new_physics`.

---

## 6) Experimental replication & benchmarking plan

### 6.1 Data curation (lightweight MVP)

* Choose a **small panel (e.g., 20–50)** protein–ligand complexes with **PDB structures + Ki/Kd** for the exact ligand/construct where possible. Store:

  * PDB ID, chain(s), metal(s), assay **T**, buffer (if known), reported **Ki/Kd** (units), computed **ΔG°**.
  * Canonical **SMILES** and **protonation** at pH ~7.4 (record assumptions).

### 6.2 Mapping to ΔG°

* Store **T** (default 298.15 K if unknown).
  Convert **Ki/Kd → ΔG°** with ( \Delta G^\circ = RT \ln(K_d/C^\circ) ).
  For pKi: ( K_i = 10^{-pK_i},\text{M} ).

### 6.3 Baselines & variants

* **Baseline A:** smina/Vina best pose score.
* **Baseline B:** smina/Vina → short **OpenMM** relax (your current NVT/NPT mini-MD) → **OpenMM energy** (consistent force field).
* **Candidate:** **new_physics** scorer on smina poses (optionally after short relax).

### 6.4 Metrics

* **Primary:** Pearson r, Spearman ρ, RMSE of **ΔG_pred vs ΔG_exp**.
* **Secondary:** Top-k enrichment, AUC (if decoys), **success fraction** (|ΔG error| ≤ X).
* **Ops:** wall-time per complex; reproducibility across seeds.

---

## 7) Testing strategy (pytest)

Add to your existing green suite:

1. **`tests/test_scoring_interface.py`**

   * Asserts `BaseScorer` contract; smoke-tests `new_physics` with a toy receptor + ligand.
2. **`tests/test_metal_geometry.py`**

   * Unit tests for Fe–N/O distance/angle terms from synthetic coordinates.
3. **`tests/test_ki_mapping.py`**

   * Unit tests for **Ki/pKi → ΔG°** conversions at various T.
4. **`tests/test_calibration_dataset.py`**

   * Loads small calibration set; verifies pipeline runs end-to-end and emits JSON metrics.
5. **`tests/test_engine_integration.py`**

   * Runs `engine=smina`, `scorer=new_physics` vs `baseline_vina` on 1–2 targets; checks consistent artifacts.
6. **`tests/test_regression_new_physics.py`**

   * Locks expected **metric ranges** to catch accidental degradations.

**Artifacts**

* `results/bench/*.csv|json` per run (pose energies, ΔG_pred, ΔG_exp, residuals).
* Rich logs (already present), plus **seed**, **config snapshot**, **package versions**.

---

## 8) Implementation plan (milestones)

**M1 — Plug-in scaffolding (1–2 days)**

* Define `BaseEngine`, `BaseScorer`; wrap existing smina path as `SminaEngine` + `VinaScorer`.
* Config switch + CLI flag passthrough; update docs.

**M2 — NewPhysicsScorer v0 (3–5 days)**

* Implement **metal coordination** term + **screened Coulomb + LJ**.
* Add **torsion strain** from prepared ligand.
* Unit tests for geometry/energy pieces.

**M3 — GB-like & polarization terms (3–5 days)**

* Add approximate desolvation + simple inducible dipole (param tables).
* Calibrate **weights** on small calibration set (grid search / Optuna optional).

**M4 — Replication harness (2–4 days)**

* Add **benchmark runner** (`ovmpk bench run …`) to compute ΔG_pred vs ΔG_exp with CSV outputs.
* Pytest calibration + regression tests; HTML/Markdown summary report.

**M5 — Comparative report (2–3 days)**

* Fixed benchmark subset; run Baseline A/B vs **new_physics**; summarize metrics, violin plots of residuals.

---

## 9) Risks & mitigations

* **Ki context mismatch:** Assay conditions may differ (T, pH, ions).
  → Store metadata; where unknown, **flag**; expect variance; compare relatively (rank) as well as absolute (RMSE).
* **Overfitting weights:**
  → Train/val split; lock hyperparams; run blind test set; keep seeds/configs versioned.
* **Metal chemistry edge cases:**
  → Start with **heme Fe** (your immediate use case); generalize later.
* **Runtime:**
  → MVP scorer should be **milliseconds per pose**; if slow, batch and vectorize; consider OpenMM CustomForce re-score.

---

## 10) Configuration & CLI (proposed)

**YAML snippet**

```yaml
benchmark:
  dataset_csv: data/bench/kiset.csv   # PDBID, ligand_id, Ki, T(K), notes...
  temperature_default: 298.15
  metrics_out: results/bench/metrics.csv

docking:
  engine: smina
  scorer: new_physics
  num_seeds: 5
  poses: 20
  box:
    auto_center: heme_fe
    size: [12,12,12]

scoring:
  new_physics: { ... as above ... }
```

**CLI**

```
ovmpk bench run --config configs/prod_test.yaml \
  --dataset data/bench/kiset.csv \
  --out results/bench/run-<date> \
  --engine smina --scorer baseline_vina
ovmpk bench run --config configs/prod_test.yaml \
  --dataset data/bench/kiset.csv \
  --out results/bench/run-<date> \
  --engine smina --scorer new_physics
```

---

## 11) What we already have (leveraged components)

* **Fetching** (RCSB APIs), **prep** (PDBFixer, OpenBabel/Dimorphite), **docking** (smina), **pose selection** (heme-aware Fe–donor), **MD** (OpenMM + OpenFF/SMIRNOFF), **analysis** (mdtraj/OpenMM), **rich logging & pytest**.
* Passing, reproducible end-to-end tests; environment seeds; detailed logs; artifacts.

The new work *slots cleanly* into this: the **engine/scorer** seam sits between docking and pose selection/MD; the **benchmark harness** parallels your existing test orchestration.

---

## 12) Multi-audience framing

### 12.1 Bioinformatics

* You get a **reproducible, scripted workflow** from sequence/structure fetch to **quantitative ΔG comparison**, with versioned configs and tests. The **scorer plug-in** model lets you test **physics variants** like any other algorithmic module.

### 12.2 Pharmaceutical (drug developer)

* The pipeline turns docking into a **calibrated decision tool** tied to **experimental Ki/ΔG**. You can **A/B** new physics against established baselines on your targets (e.g., heme P450s), quantify deltas, and **automate reports** for project meetings and go/no-go triage.

### 12.3 Computer science (full-stack / infra)

* Clear **interfaces**, **unit/integration tests**, **artifacted runs**, and **config-driven** behavior. The new scorer is a **pluggable strategy** with well-defined IO, enabling micro-benchmarking and profiling. Easy to add **ML-based scorers** later (e.g., SE(3) GNN) without breaking the pipeline.

### 12.4 Physician (pediatric focus)

* This is a **pre-clinical modeling tool** aiming to **improve how we shortlist candidate molecules** before wet-lab work, potentially accelerating **safer, targeted therapies** for pediatric disease. It doesn’t replace experiments; it **prioritizes** them more intelligently.

**Common ground:** rigorous, reproducible, and comparable results (ΔG/Ki), with transparent logs and tests.

---

## 13) Deliverables (MVP)

* **Code**

  * `ovmpk/docking/engines/base.py`, `smina_engine.py`
  * `ovmpk/docking/scorers/base.py`, `vina_scorer.py`, `new_physics.py`
  * `ovmpk/bench/runner.py`
* **Configs & Data**

  * `configs/prod_test.yaml` (extended), `data/bench/kiset.csv` (small curated set)
* **Tests**

  * As listed in §7 (new test files)
* **Docs**

  * README update with **how-to run replication and A/B compare**
  * A **Method Card** for NewPhysicsScorer: equations, assumptions, params

---

## 14) Example pytest additions

```
pytest -q tests/test_scoring_interface.py::test_contract
pytest -q tests/test_metal_geometry.py
pytest -q tests/test_calibration_dataset.py -s --rich
pytest -q tests/test_regression_new_physics.py --maxfail=1
```

---

## 15) Roadmap after MVP

* Add **water-placement** (grand canonical or knowledge-based waters) term.
* Optionally re-score with **OpenMM Custom*Force** for GPU speed.
* Add **ML re-ranker** (equivariant GNN) as another `BaseScorer` implementation.
* Expand to **non-heme metals** and **non-metallo** targets; microstate ensembles.

---

### TL;DR (what changes now)

* You’ll add a **scorer plug-in** seam, implement **NewPhysicsScorer v0** (metal-aware + polarization + fast desolvation), and build a **benchmark harness** that converts experimental **Ki → ΔG°** and compares directly to **ΔG_pred**—all inside your existing, already-tested OVM-PK pipeline with rich logs and CI-style tests.
