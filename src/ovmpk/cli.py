from pathlib import Path
import os, json, time, platform, yaml, typer

app = typer.Typer(help="OVM-PK pipeline")

def _write_manifest(outdir: Path, params: dict):
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "host": platform.node(),
        "params": params,
        "timings": {},
        "dry_run": bool(os.getenv("OVM_DRY_RUN", ""))
    }
    (outdir / "manifest.json").write_text(json.dumps(manifest, indent=2))

def pipeline_mvp(config: str, gene: str, ligand: str) -> str:
    cfg = yaml.safe_load(Path(config).read_text())
    # Phase 1: fetch
    from ovmpk.fetchers import protein_fetcher, ligand_fetcher
    p = protein_fetcher.fetch(gene, cfg)
    l = ligand_fetcher.fetch(ligand, cfg)
    # Phase 2: prep
    from ovmpk.prep import protein_prep, ligand_prep
    prot = protein_prep.prepare(p, cfg)
    ligp = ligand_prep.prepare(l, cfg)
    # Phase 3: docking
    from ovmpk.docking import smina_wrapper, pose_selector
    poses = smina_wrapper.run(prot, ligp, cfg)
    best = pose_selector.pick(poses, cfg)
    outdir = Path(cfg["reporting"]["outdir"]) / f"{Path(prot).stem}_{ligand}"
    _write_manifest(outdir, cfg)
    (outdir / "summary.json").write_text(json.dumps({"best_pose": best}, indent=2))
    return outdir.as_posix()

@app.command("run")
def run(config: str = "configs/default.yaml", gene: str = "CYP3A4", ligand: str = "ketoconazole"):
    out = pipeline_mvp(config, gene, ligand)
    typer.echo(out)

# Test helper (imported in tests)
def _smoke_run():
    return run

if __name__ == "__main__":
    app()
