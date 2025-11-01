#!/usr/bin/env python3
import argparse, os, sys, json, math, random, time, shutil, uuid, subprocess
from pathlib import Path
from typing import List, Tuple, Optional
import yaml

from ovmpk.utils.run_dirs import stage_dir

# ---------- Small types ----------

class Candidate:
    def __init__(self,
                 center: Tuple[float,float,float],
                 size: float,
                 min_fe_n: float = float("inf"),
                 best_pose_idx: Optional[int] = None,
                 best_sdf: Optional[str] = None,
                 best_alias: Optional[str] = None,
                 suggested_center: Optional[Tuple[float,float,float]] = None,
                 exhaustiveness: Optional[int] = None,
                 poses: Optional[int] = None,
                 seeds: Optional[int] = None):
        self.center = center
        self.size = size
        self.min_fe_n = min_fe_n
        self.best_pose_idx = best_pose_idx
        self.best_sdf = best_sdf
        self.best_alias = best_alias
        self.suggested_center = suggested_center
        self.exhaustiveness = exhaustiveness
        self.poses = poses
        self.seeds = seeds

# ---------- Utils ----------

def ensure_dir(p: Path) -> None:
    p.mkdir(parents=True, exist_ok=True)

def nowstamp() -> str:
    return time.strftime("%Y-%m-%d_%H-%M-%S")

def load_yaml(p: Path) -> dict:
    with p.open() as fh:
        return yaml.safe_load(fh)

def dump_yaml(d: dict, p: Path) -> None:
    with p.open("w") as fh:
        yaml.safe_dump(d, fh, sort_keys=False)

def write_tsv(rows: List[dict], path: Path) -> None:
    if not rows:
        Path(path).write_text("")
        return
    keys = list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")

def clamp(v: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, v))

def clamp_box_size(sz: float, lo: float = 8.0, hi: float = 30.0) -> float:
    return clamp(float(sz), lo, hi)

def jitter_center(c: Tuple[float,float,float], sigma: float) -> Tuple[float,float,float]:
    return (c[0] + random.gauss(0, sigma),
            c[1] + random.gauss(0, sigma),
            c[2] + random.gauss(0, sigma))

# ---------- RDKit scoring helpers ----------

def _rdkit_imports():
    from rdkit import Chem
    import numpy as np
    return Chem, np

def fe_n_distance_for_sdf(sdf_path: Path, fe_xyz: Tuple[float,float,float]) -> Tuple[float, Optional[int], Optional[Tuple[float,float,float]]]:
    """
    Returns (best_min_distance, best_pose_index, suggested_center_midpoint)
    suggested_center = midpoint(heme Fe, closest ligand N in that pose)
    """
    Chem, np = _rdkit_imports()
    if not sdf_path.exists() or sdf_path.stat().st_size == 0:
        return float("inf"), None, None

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    best_d = float("inf")
    best_i = None
    best_mid = None
    FE = np.array(fe_xyz, dtype=float)

    for i, m in enumerate(suppl, start=1):
        if m is None:
            continue
        conf = m.GetConformer()
        dmin = float("inf")
        best_n = None
        for a in m.GetAtoms():
            if a.GetSymbol() == "N":
                p = conf.GetAtomPosition(a.GetIdx())
                v = np.array([p.x, p.y, p.z])
                d = float(np.linalg.norm(v - FE))
                if d < dmin:
                    dmin = d
                    best_n = v
        if dmin < best_d:
            best_d = dmin
            best_i = i
            if best_n is not None:
                mid = (FE + best_n) / 2.0
                best_mid = (float(mid[0]), float(mid[1]), float(mid[2]))

    return best_d, best_i, best_mid

# ---------- Candidate proposal ----------

def propose_population(elites: List[Candidate], pop: int, sigma: float, size: float,
                       fe_xyz: Tuple[float,float,float]) -> List[Candidate]:
    out: List[Candidate] = []
    if not elites:
        elites = [Candidate(center=fe_xyz, size=size)]
    out.extend(Candidate(center=e.center, size=size) for e in elites)
    while len(out) < pop:
        base = random.choice(elites)
        out.append(Candidate(center=jitter_center(base.center, sigma), size=size))
    return out[:pop]

# ---------- smina/ovmpk plumbing ----------

def build_config_from_base(base_cfg: dict,
                           center: Tuple[float,float,float],
                           edge: float,
                           exhaustiveness: int,
                           poses: int,
                           seed: int,
                           cpu_threads_per_job: int) -> dict:
    cfg = yaml.safe_load(yaml.safe_dump(base_cfg))  # deep copy-ish
    cfg.setdefault("docking", {}).setdefault("box", {})
    cfg["docking"]["box"]["center"] = [float(center[0]), float(center[1]), float(center[2])]
    cfg["docking"]["box"]["size"] = [float(edge), float(edge), float(edge)]
    cfg["docking"]["exhaustiveness"] = int(exhaustiveness)
    cfg["docking"]["poses"] = int(poses)
    cfg["docking"]["seed"] = int(seed)
    cfg.setdefault("runtime", {})["cpu_threads"] = int(max(1, cpu_threads_per_job))
    return cfg

def run_one_candidate(gene: str,
                      ligand: str,
                      base_cfg: dict,
                      fe_xyz: Tuple[float,float,float],
                      center: Tuple[float,float,float],
                      edge: float,
                      exhaustiveness: int,
                      poses: int,
                      seeds: int,
                      alias_prefix: str,
                      cap_env: Optional[dict] = None) -> Candidate:
    work_glob_dir = stage_dir("docking")
    tmp_cfg_dir = stage_dir("tmp_cfgs")
    best_cand = Candidate(center=center, size=edge,
                          exhaustiveness=exhaustiveness, poses=poses, seeds=seeds)

    for _ in range(seeds):
        seed_val = random.randint(1, 10_000_000)
        alias = f"{alias_prefix}_s{seed_val:x}"
        cfg = build_config_from_base(
            base_cfg, center, edge, exhaustiveness, poses, seed_val,
            int(base_cfg.get("runtime", {}).get("cpu_threads", 1))
        )
        tmp_yaml = tmp_cfg_dir / f"{alias}.yaml"
        ensure_dir(tmp_yaml.parent)
        dump_yaml(cfg, tmp_yaml)

        env = os.environ.copy()
        env["OVM_ALIAS"] = alias
        if cap_env:
            env.update(cap_env)

        cmd = ["ovmpk", "--config", str(tmp_yaml), "--gene", gene, "--ligand", ligand]
        p = subprocess.run(cmd, text=True, capture_output=True, env=env)
        if p.returncode != 0:
            sys.stderr.write(f"[warn] ovmpk failed for {alias}\n{p.stderr}\n")
            continue

        sdf_list = sorted(work_glob_dir.glob(f"*{alias}_poses.sdf"))
        if not sdf_list:
            sdf_list = sorted(work_glob_dir.glob(f"*__{alias}_poses.sdf"))
        if not sdf_list:
            sys.stderr.write(f"[warn] poses SDF missing for {alias}\n")
            continue

        sdf = sdf_list[-1]
        d, idx, mid = fe_n_distance_for_sdf(sdf, fe_xyz)
        if d < best_cand.min_fe_n:
            best_cand.min_fe_n = d
            best_cand.best_pose_idx = idx
            best_cand.best_sdf = sdf.as_posix()
            best_cand.best_alias = alias
            best_cand.suggested_center = mid

    return best_cand

def evaluate_generation(gene: str,
                        ligand: str,
                        base_cfg: dict,
                        cand_list: List[Candidate],
                        seeds: int,
                        jobs: int,
                        fe_xyz: Tuple[float,float,float],
                        worktag: str,
                        cpu_threads_per_job: int,
                        cap_env: Optional[dict]) -> List[Candidate]:
    base_cfg = yaml.safe_load(yaml.safe_dump(base_cfg))
    base_cfg.setdefault("runtime", {})["cpu_threads"] = int(max(1, cpu_threads_per_job))

    results: List[Candidate] = []
    from concurrent.futures import ThreadPoolExecutor, as_completed
    with ThreadPoolExecutor(max_workers=jobs) as ex:
        futs = []
        for c in cand_list:
            alias_prefix = f"{ligand}__{worktag}_{uuid.uuid4().hex[:6]}"
            futs.append(ex.submit(
                run_one_candidate,
                gene, ligand, base_cfg, fe_xyz,
                c.center, c.size,
                int(base_cfg["docking"]["exhaustiveness"]),
                int(base_cfg["docking"]["poses"]),
                seeds,
                alias_prefix,
                cap_env
            ))
        for f in as_completed(futs):
            results.append(f.result())

    return results

# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser(description="Parallel evolutionary pocket search for ovmpk+smina")
    ap.add_argument("--ligand", required=True, help="Ligand base name (data/input/ligands/<name>.sdf)")
    ap.add_argument("--gene", required=True, help="Protein gene (e.g., CYP3A4)")
    ap.add_argument("--fe", nargs=3, type=float, required=True, metavar=("X", "Y", "Z"),
                    help="Heme Fe coordinates (Å) for scoring / initial center")
    ap.add_argument("--center", nargs=3, type=float, metavar=("X", "Y", "Z"),
                    help="Optional start center; defaults to --fe")
    ap.add_argument("--base-config", default="configs/heme_focus.yaml",
                    help="Base YAML config to modify per job")

    ap.add_argument("--gens", type=int, default=3, help="Generations")
    ap.add_argument("--pop", type=int, default=10, help="Population per generation")
    ap.add_argument("--seeds", type=int, default=3, help="Docking seeds per candidate")
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 8)//2),
                    help="Parallel jobs (candidates evaluated concurrently)")
    ap.add_argument("--cpu-threads-per-job", type=int, default=0,
                    help="Override per-job smina CPU threads. 0 = auto-scale.")
    ap.add_argument("--threads-cap", type=int, default=4,
                    help="Max threads per job when auto-scaling (ignored if --cpu-threads-per-job>0).")
    ap.add_argument("--cap-omp", action="store_true",
                    help="Set OMP/BLAS threading env vars (OMP, MKL, OPENBLAS, etc.) to 1 for child runs.")

    ap.add_argument("--start-size", type=float, default=18.0, help="Start box edge (Å)")
    ap.add_argument("--shrink", type=float, default=0.8, help="Per-generation shrink factor for size/sigma")
    ap.add_argument("--sigma-frac", type=float, default=0.35,
                    help="Exploration sigma = sigma_frac * box_size")
    ap.add_argument("--topk", type=int, default=None,
                    help="How many elites seed the next generation (default pop//3, min 2)")
    ap.add_argument("--outdir", default=None,
                    help="Output folder for summaries (default results/autopocket/<ligand>_<timestamp>)")

    ap.add_argument("--exh-sched", default="16,16,32,64",
                    help="Comma list of exhaustiveness per generation")
    ap.add_argument("--poses-sched", default="10,10,20,40",
                    help="Comma list of poses per generation")
    ap.add_argument("--stop-fen", type=float, default=None,
                    help="Stop early if best Fe–N <= this Å")

    # ---- Adaptive mode ----
    ap.add_argument("--adaptive", action="store_true",
                    help="Adapt jobs (and optionally seeds) as box size shrinks.")
    ap.add_argument("--adaptive-alpha", type=float, default=1.0,
                    help="Jobs scale factor exponent; jobs ~= base_jobs*(start_size/current_size)^alpha.")
    ap.add_argument("--jobs-min", type=int, default=1,
                    help="Lower bound for adaptive jobs.")
    ap.add_argument("--jobs-max", type=int, default=0,
                    help="Upper bound for adaptive jobs (0 = auto to CPU threads or CPU//threads_per_job).")
    ap.add_argument("--adaptive-seeds", action="store_true",
                    help="Increase seeds as size shrinks.")
    ap.add_argument("--adaptive-seeds-beta", type=float, default=0.5,
                    help="Seeds ~= base_seeds*(start_size/current_size)^beta.")
    ap.add_argument("--seeds-max", type=int, default=6,
                    help="Upper bound for seeds when adaptive seeds enabled.")

    args = ap.parse_args()

    fe_xyz = (args.fe[0], args.fe[1], args.fe[2])
    start_center = tuple(args.center) if args.center else fe_xyz

    base_cfg_path = Path(args.base_config)
    if not base_cfg_path.exists():
        print(f"[FATAL] Base config not found: {base_cfg_path}", file=sys.stderr)
        sys.exit(2)
    base_cfg = load_yaml(base_cfg_path)

    # Hardware + defaults
    hw_threads = int(os.cpu_count() or 8)

    # Derive per-job threads if not provided
    fixed_threads = args.cpu_threads_per_job if args.cpu_threads_per_job > 0 else 0
    if fixed_threads:
        cpu_threads_per_job0 = int(max(1, fixed_threads))
        auto_note = ""
    else:
        # initial guess for gen1 when not adaptive
        cpu_threads_per_job0 = max(1, min(args.threads-cap if hasattr(args, "threads-cap") else args.threads_cap,
                                          hw_threads // max(1, args.jobs)))
        auto_note = " (auto)"

    # Derive jobs_max default
    if args.jobs_max > 0:
        jobs_max_default = args.jobs_max
    else:
        if fixed_threads:
            jobs_max_default = max(1, hw_threads // cpu_threads_per_job0)
        else:
            jobs_max_default = hw_threads

    # Optional library thread caps for child runs
    cap_env = None
    if args.cap_omp:
        cap_env = {
            "OMP_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "RDKIT_NUM_THREADS": "1",
            "NUMEXPR_NUM_THREADS": "1",
        }

    # Output root
    out_root = Path(args.outdir) if args.outdir else Path("results") / "autopocket" / f"{args.ligand}_{nowstamp()}"
    ensure_dir(out_root)

    # Schedules
    exh_sched = [int(x) for x in args.exh_sched.split(",") if x.strip()]
    poses_sched = [int(x) for x in args.poses_sched.split(",") if x.strip()]

    # Evolution parameters
    curr_size = float(args.start_size)
    sigma = args.sigma_frac * curr_size
    topk = args.topk if args.topk is not None else max(2, args.pop // 3)
    elites: List[Candidate] = [Candidate(center=start_center, size=curr_size)]

    # Planner helpers
    def plan_jobs_threads(size_now: float) -> Tuple[int,int]:
        if not args.adaptive:
            return int(args.jobs), int(cpu_threads_per_job0)
        ratio = max(1e-6, args.start_size / max(1e-6, size_now))
        jobs_target = int(round(args.jobs * (ratio ** args.adaptive_alpha)))
        if fixed_threads:
            # Respect fixed per-job threads; keep total <= hw if possible
            jobs_hi = jobs_max_default
            jobs_cur = int(clamp(jobs_target, args.jobs_min, jobs_hi))
            threads_cur = int(cpu_threads_per_job0)
        else:
            jobs_hi = jobs_max_default
            jobs_cur = int(clamp(jobs_target, args.jobs_min, jobs_hi))
            threads_cur = int(max(1, min(args.threads_cap, hw_threads // max(1, jobs_cur))))
        return jobs_cur, threads_cur

    def plan_seeds(size_now: float) -> int:
        if not args.adaptive_seeds:
            return int(args.seeds)
        ratio = max(1e-6, args.start_size / max(1e-6, size_now))
        seeds_target = int(round(args.seeds * (ratio ** args.adaptive_seeds_beta)))
        return int(clamp(seeds_target, 1, args.seeds_max))

    print(f"[info] Start center = {start_center}, size = {curr_size:.2f} Å  (sigma={sigma:.2f})")
    print(f"[info] Generations={args.gens}, pop={args.pop}, topk={topk}")
    print(f"[info] Schedules: exhaustiveness={exh_sched}, poses={poses_sched}")
    if args.adaptive:
        print(f"[info] ADAPTIVE ON — alpha={args.adaptive_alpha}, jobs_min={args.jobs_min}, jobs_max={jobs_max_default}")
        if args.adaptive_seeds:
            print(f"[info] ADAPTIVE SEEDS — beta={args.adaptive_seeds_beta}, seeds_max={args.seeds_max}")
    else:
        print(f"[info] jobs={args.jobs}, per-job threads={cpu_threads_per_job0}{auto_note} (hw={hw_threads})")
    if args.cap_omp:
        print(f"[info] Capping OMP/BLAS env threads for child runs: {cap_env}")
    print(f"[info] Results -> {out_root}\n")

    all_rows: List[dict] = []

    for gen in range(1, args.gens + 1):
        worktag = f"gen{gen}"
        exh = exh_sched[min(gen-1, len(exh_sched)-1)]
        poses = poses_sched[min(gen-1, len(poses_sched)-1)]
        seeds_cur = plan_seeds(curr_size)

        jobs_cur, threads_cur = plan_jobs_threads(curr_size)

        print(f"\n=== Generation {gen}/{args.gens} ===")
        print(f"[info] box size={curr_size:.2f} Å, sigma={sigma:.2f}")
        print(f"[info] proposing {args.pop} candidates around {len(elites)} elite(s)")
        print(f"[info] per-gen: exhaustiveness={exh}, poses={poses}, seeds={seeds_cur}")
        print(f"[info] resources: jobs={jobs_cur}, per-job threads={threads_cur} (hw={hw_threads})")

        # per-gen base config
        base_cfg_gen = yaml.safe_load(yaml.safe_dump(base_cfg))
        base_cfg_gen.setdefault("docking", {})
        base_cfg_gen["docking"]["exhaustiveness"] = exh
        base_cfg_gen["docking"]["poses"] = poses
        base_cfg_gen.setdefault("runtime", {})["cpu_threads"] = int(threads_cur)

        # Propose + evaluate
        cands = propose_population(elites, args.pop, sigma, curr_size, fe_xyz)
        cands = evaluate_generation(
            gene=args.gene,
            ligand=args.ligand,
            base_cfg=base_cfg_gen,
            cand_list=cands,
            seeds=seeds_cur,
            jobs=jobs_cur,
            fe_xyz=fe_xyz,
            worktag=worktag,
            cpu_threads_per_job=threads_cur,
            cap_env=cap_env,
        )

        # Summarize
        rows = []
        for i, c in enumerate(sorted(cands, key=lambda x: x.min_fe_n)):
            rows.append({
                "gen": gen,
                "rank": i + 1,
                "center_x": round(c.center[0], 3),
                "center_y": round(c.center[1], 3),
                "center_z": round(c.center[2], 3),
                "size": round(c.size, 2),
                "min_FeN_A": round(c.min_fe_n, 3) if math.isfinite(c.min_fe_n) else "inf",
                "best_pose_idx": c.best_pose_idx,
                "best_sdf": c.best_sdf or "",
                "alias": c.best_alias or "",
                "suggested_center": ",".join(f"{v:.3f}" for v in c.suggested_center) if c.suggested_center else "",
                "exhaustiveness": exh,
                "poses": poses,
                "seeds": seeds_cur,
                "jobs": jobs_cur,
                "threads": threads_cur,
            })
        all_rows.extend(rows)

        tsv_path = out_root / f"gen{gen}_summary.tsv"
        write_tsv(rows, tsv_path)
        print(f"[info] wrote {tsv_path}")

        # Elites -> next gen
        elites = sorted(cands, key=lambda x: x.min_fe_n)[:max(2, topk)]
        best_now = elites[0].min_fe_n if elites else float("inf")
        print(f"[info] best min Fe–N this gen = {best_now:.3f} Å")

        for e in elites:
            if e.suggested_center:
                e.center = e.suggested_center

        if args.stop_fen is not None and best_now <= float(args.stop_fen):
            print(f"[info] early stop: best Fe–N ≤ {args.stop_fen} Å")
            break

        # Shrink search
        curr_size = clamp_box_size(curr_size * args.shrink, lo=8.0)
        sigma = max(0.5, args.sigma_frac * curr_size)

    # Overall summary
    overall_tsv = out_root / "all_generations.tsv"
    write_tsv(all_rows, overall_tsv)
    with (out_root / "summary.json").open("w") as fh:
        json.dump({"params": vars(args), "outdir": str(out_root)}, fh, indent=2)

    # Print best overall
    best_overall = None
    for r in all_rows:
        if isinstance(r["min_FeN_A"], str):
            continue
        if (best_overall is None) or (r["min_FeN_A"] < best_overall["min_FeN_A"]):
            best_overall = r
    if best_overall:
        print("\n=== BEST OVERALL ===")
        print(json.dumps(best_overall, indent=2))
    else:
        print("\n[warn] No successful candidates to report.")

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("\n[abort] Interrupted by user.", file=sys.stderr)
        sys.exit(130)

