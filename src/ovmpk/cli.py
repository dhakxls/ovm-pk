from pathlib import Path
import os
import json
import time
import platform
import yaml
import typer
import traceback # Import traceback

# Initialize the Typer application
app = typer.Typer(help="OVM-PK pipeline: Docking -> MD -> Free Energy")

def _write_manifest(outdir: Path, params: dict):
    """Helper function to write a manifest file with run parameters."""
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "host": platform.node(),
        "params": params,
        "timings": {},
        "dry_run": bool(os.getenv("OVM_DRY_RUN", ""))
    }
    manifest_path = outdir / "manifest.json"
    try:
        manifest_path.write_text(json.dumps(manifest, indent=2))
        print(f"[info] Wrote manifest: {manifest_path}")
    except Exception as e:
        print(f"[error] Failed to write manifest {manifest_path}: {e}")


def pipeline_mvp(config: str, gene: str, ligand: str) -> str:
    """Core logic for the Minimum Viable Pipeline."""
    print(f"[info] Starting pipeline with config='{config}', gene='{gene}', ligand='{ligand}'")
    config_path = Path(config)
    if not config_path.exists():
        raise FileNotFoundError(f"Configuration file not found: {config}")

    try:
        cfg = yaml.safe_load(config_path.read_text())
        print("[info] Loaded configuration.")
    except Exception as e:
        raise ValueError(f"Error loading configuration file {config}: {e}")

    # Ensure reporting keys exist
    if "reporting" not in cfg or "outdir" not in cfg["reporting"]:
         raise ValueError("Configuration missing 'reporting.outdir' key.")

    # Phase 1: fetch
    print("[info] Phase 1: Fetching inputs...")
    # It's generally better practice to import at the top,
    # but keep imports inside function if needed for specific reasons (like avoiding circular imports)
    from ovmpk.fetchers import protein_fetcher, ligand_fetcher
    protein_paths = protein_fetcher.fetch(gene, cfg) # Returns a dict like {"apo": Path(...)}
    ligand_path = ligand_fetcher.fetch(ligand, cfg)
    print(f"[info] Fetched protein: {protein_paths}, Ligand: {ligand_path}")

    # Phase 2: prep
    print("[info] Phase 2: Preparing inputs...")
    from ovmpk.prep import protein_prep, ligand_prep
    # Assuming prepare expects the path, not the dict for now
    protein_input_path = protein_paths.get("apo") # Adjust if policy changes
    if not protein_input_path:
         raise ValueError(f"Could not determine protein input path from fetched paths: {protein_paths}")
    # Ensure protein_prep.prepare receives what it expects. If it needs a dict:
    prepared_protein_path = protein_prep.prepare({"apo": protein_input_path}, cfg)
    # If it just needs the path string/Path object:
    # prepared_protein_path = protein_prep.prepare(protein_input_path, cfg) # Check protein_prep.prepare definition
    prepared_ligand_path = ligand_prep.prepare(ligand_path, cfg)
    print(f"[info] Prepared Protein: {prepared_protein_path}, Ligand: {prepared_ligand_path}")


    # Phase 3: docking
    print("[info] Phase 3: Running docking...")
    from ovmpk.docking import smina_wrapper, pose_selector
    # Ensure paths are strings for subprocess calls if needed
    pose_sdf_paths = smina_wrapper.run(str(prepared_protein_path), str(prepared_ligand_path), cfg)
    print(f"[info] Docking generated poses: {pose_sdf_paths}")
    best_pose_path = pose_selector.pick(pose_sdf_paths, cfg) # Should return a path string or None
    print(f"[info] Selected best pose: {best_pose_path}")

    # --- Reporting ---
    # Use the *prepared* protein path stem for naming consistency if prep modifies it
    # If prep just passes through, Path(prepared_protein_path).stem is fine.
    protein_stem = Path(prepared_protein_path).stem
    outdir = Path(cfg["reporting"]["outdir"]) / f"{protein_stem}_{ligand}"
    print(f"[info] Writing results to: {outdir}")
    _write_manifest(outdir, cfg)

    summary_data = {
        "config_file": config,
        "gene": gene,
        "ligand": ligand,
        "fetched_protein_paths": {k: str(v) for k, v in protein_paths.items()},
        "fetched_ligand_path": str(ligand_path),
        "prepared_protein_path": str(prepared_protein_path),
        "prepared_ligand_path": str(prepared_ligand_path),
        "docking_pose_sdf_paths": pose_sdf_paths,
        "selected_best_pose_path": best_pose_path,
        "output_directory": str(outdir)
        }
    summary_path = outdir / "summary.json"
    try:
        summary_path.write_text(json.dumps(summary_data, indent=2))
        print(f"[info] Wrote summary: {summary_path}")
    except Exception as e:
        print(f"[error] Failed to write summary {summary_path}: {e}")


    return str(outdir) # Return path as string

# Define the 'run' command
# Using @app.command() implicitly uses the function name "run" as the command name
@app.command() # MODIFIED: Removed explicit "run" string
def run(
    config: str = typer.Option("configs/default.yaml", help="Path to the YAML configuration file."),
    gene: str = typer.Option("CYP3A4", help="Gene name for protein fetching."),
    ligand: str = typer.Option("ketoconazole", help="Ligand name for fetching/docking.")
):
    """
    Runs the OVM-PK pipeline: Fetch -> Prep -> Dock.
    """
    try:
        output_directory = pipeline_mvp(config, gene, ligand)
        typer.echo(f"Pipeline finished successfully. Results in: {output_directory}")
    except FileNotFoundError as e:
        typer.echo(f"Error: Input file not found - {e}", err=True)
        raise typer.Exit(code=1)
    except ValueError as e:
        typer.echo(f"Error: Configuration or data issue - {e}", err=True)
        raise typer.Exit(code=1)
    except RuntimeError as e:
        # Catch errors from external tools like smina/obabel
        typer.echo(f"Error: Runtime issue (external tool failed?) - {e}", err=True)
        raise typer.Exit(code=1)
    except Exception as e:
        # Catch any other unexpected errors
        typer.echo(f"An unexpected error occurred: {e}", err=True)
        traceback.print_exc() # Print traceback for debugging
        raise typer.Exit(code=1)


# Test helper (can be imported in tests)
def _smoke_run():
    # This is slightly awkward for testing, maybe adjust test later
    # For now, just return the function reference
    return run

# Main entry point for running the script directly
if __name__ == "__main__":
    app()

