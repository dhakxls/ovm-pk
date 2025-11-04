"""Interactive and scripted analysis workflow for OVM-PK.

Provides interactive CLI and programmatic commands for fetching benchmarks,
running physics scoring, and scanning energy weights.
"""

import argparse
import json
import os
import re
from copy import deepcopy
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import questionary
from rich.console import Console
from rich.panel import Panel
from rich.progress import Progress, SpinnerColumn, TextColumn

from ovmpk import Pipeline
from ovmpk.config import load_config, update_config
from ovmpk.physics.mmgbsa import MMGBSACalculator
from ovmpk.cli.integration import (
    fetch_protein_structure,
    fetch_ligand_structure,
    create_config,
    run_analysis as run_analysis_integrated,
)
from ovmpk.fetchers.chembl_client import ChEMBLClient


def build_parser() -> argparse.ArgumentParser:
    """Create the argument parser for command-line execution."""
    parser = argparse.ArgumentParser(description="OVM-PK analysis utilities")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run the physics scoring pipeline")
    run_parser.add_argument("--config", "-c", type=Path, required=True, help="Pipeline config YAML")

    chembl_parser = subparsers.add_parser(
        "chembl", help="Fetch bioactivity benchmarks from ChEMBL"
    )
    chembl_parser.add_argument("target", help="Target ChEMBL ID, e.g. CHEMBL262")
    chembl_parser.add_argument(
        "--out", "-o", type=Path, default=Path("chembl_activities.json"), help="Output JSON file"
    )
    chembl_parser.add_argument(
        "--max-records", type=int, default=1000, help="Maximum number of activities to retrieve"
    )
    chembl_parser.add_argument(
        "--standard-types",
        nargs="*",
        default=["Ki", "IC50", "Kd"],
        help="Filter to these standard types",
    )
    chembl_parser.add_argument(
        "--assay-type",
        default="B",
        help="Assay type filter (default binding 'B')",
    )

    scan_parser = subparsers.add_parser(
        "scan-weights", help="Run weight-scan defined in pipeline config"
    )
    scan_parser.add_argument("--config", "-c", type=Path, required=True)
    scan_parser.add_argument(
        "--output", "-o", type=Path, default=Path("weight_scan_results.json")
    )

    delta_g_parser = subparsers.add_parser(
        "convert-dg", help="Compute ΔG from ChEMBL activities JSON"
    )
    delta_g_parser.add_argument("input", type=Path, help="ChEMBL activities JSON")
    delta_g_parser.add_argument(
        "--out", "-o", type=Path, default=Path("chembl_with_delta_g.json")
    )
    delta_g_parser.add_argument(
        "--temperature", type=float, default=298.15, help="Temperature in K"
    )

    map_dg_parser = subparsers.add_parser(
        "map-dg", help="Map physics totals to ΔG using a linear calibration"
    )
    map_dg_parser.add_argument("input", type=Path, help="Physics results JSON")
    map_dg_parser.add_argument(
        "--out", "-o", type=Path, default=Path("mapped_results.json")
    )
    map_dg_parser.add_argument(
        "--slope",
        type=float,
        required=True,
        help="Linear slope applied to physics total scores",
    )
    map_dg_parser.add_argument(
        "--intercept",
        type=float,
        required=True,
        help="Linear intercept added after scaling",
    )

    mmgbsa_parser = subparsers.add_parser(
        "mmgbsa", help="Compute MM/GBSA binding energies for weight scan results"
    )
    mmgbsa_parser.add_argument("input", type=Path, help="Weight scan results JSON")
    mmgbsa_parser.add_argument(
        "--out", "-o", type=Path, default=Path("mmgbsa_results.json")
    )
    mmgbsa_parser.add_argument(
        "--gbsa-model",
        choices=["GBn2", "GBn", "OBC1", "none"],
        default="GBn2",
        help="Implicit solvent model for MM/GBSA",
    )
    mmgbsa_parser.add_argument(
        "--labels",
        nargs="*",
        help="Optional subset of scan labels to evaluate",
    )
    mmgbsa_parser.add_argument(
        "--minimize",
        action="store_true",
        help="Perform a short energy minimization of each state before evaluating energies",
    )
    mmgbsa_parser.add_argument(
        "--md-steps",
        type=int,
        default=0,
        help="If >0, run Langevin MD for this many steps per state and accumulate an ensemble",
    )
    mmgbsa_parser.add_argument(
        "--sample-interval",
        type=int,
        default=50,
        help="Stride (in MD steps) between ensemble samples",
    )
    mmgbsa_parser.add_argument(
        "--relax",
        action="store_true",
        help="Relax complex/receptor/ligand snapshots before energy evaluation",
    )
    mmgbsa_parser.add_argument(
        "--relax-output",
        type=Path,
        help="Directory to write relaxed complex/receptor/ligand PDB snapshots",
    )
    mmgbsa_parser.add_argument(
        "--relax-md-steps",
        type=int,
        default=0,
        help="Number of MD steps to run during relaxation (0 = skip MD)",
    )
    mmgbsa_parser.add_argument(
        "--relax-random-seed",
        type=int,
        help="Random seed for relaxation MD sampling",
    )

    return parser


def _run_pipeline_with_config(config_path: Path) -> Dict[str, Any]:
    pipeline = Pipeline(config_path)
    return pipeline.run()


def run_cli(args: argparse.Namespace) -> None:
    if args.command == "chembl":
        client = ChEMBLClient()
        activities = client.fetch_activities(
            args.target,
            standard_types=args.standard_types,
            assay_type=args.assay_type,
            max_records=args.max_records,
        )
        payload = [record.__dict__ for record in activities]
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(payload, indent=2))
        console.print(
            Panel.fit(
                f"Saved {len(payload)} activities for {args.target} to {args.out}",
                title="ChEMBL",
            )
        )
        return

    if args.command == "scan-weights":
        # Lazy import to avoid circular dependency during CLI usage
        import yaml

        with args.config.open() as fh:
            base_config = yaml.safe_load(fh)

        scoring_section = base_config.get("scoring", {})
        weight_scan = scoring_section.get("weight_scan", [])
        if not weight_scan:
            console.print("[red]No weight_scan entries found in config[/red]")
            return

        results: List[Dict[str, Any]] = []
        for idx, overrides in enumerate(weight_scan, 1):
            console.print(f"[bold]Running scan {idx}/{len(weight_scan)}[/bold]: {overrides}")
            scan_config = deepcopy(base_config)
            scan_scoring = deepcopy(scan_config.get("scoring", {}))
            scan_scoring.pop("weight_scan", None)

            weight_overrides = {k: v for k, v in overrides.items() if k.startswith("w_")}
            param_overrides = overrides.get("parameters") or {}

            scan_parameters = dict(scan_scoring.get("parameters", {}))
            scan_parameters.update(param_overrides)

            if scan_parameters:
                scan_scoring["parameters"] = scan_parameters

            scan_scoring.update(weight_overrides)
            scan_config["scoring"] = scan_scoring

            tmp_path = args.config.parent / f"_tmp_weight_scan_{idx}.yml"
            with tmp_path.open("w") as tmp_file:
                yaml.safe_dump(scan_config, tmp_file)

            try:
                run_result = _run_pipeline_with_config(tmp_path)
                results.append({
                    "weights": weight_overrides,
                    "label": overrides.get("name", f"scan_{idx}"),
                    "result": run_result,
                })
            finally:
                tmp_path.unlink(missing_ok=True)

        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, indent=2, default=str))
        console.print(
            Panel.fit(
                f"Saved weight scan results to {args.output}", title="Weight Scan"
            )
        )
        return

    if args.command == "convert-dg":
        from ovmpk.fetchers.chembl_client import ChEMBLClient

        data = json.loads(args.input.read_text())
        client = ChEMBLClient()
        results = []
        for entry in data:
            record = client._to_activity_record(entry, args.temperature)  # type: ignore
            results.append(record.__dict__)

        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(results, indent=2))
        console.print(
            Panel.fit(
                f"Annotated ΔG for {len(results)} activities to {args.out}",
                title="ChEMBL ΔG",
            )
        )
        return

    if args.command == "map-dg":
        raw_data = json.loads(args.input.read_text())
        mapping_info = {"slope": args.slope, "intercept": args.intercept}

        def _apply_mapping(scoring: Dict[str, Any]) -> None:
            if not isinstance(scoring, dict):
                return
            for score in scoring.get("scores", []):
                total = score.get("total_score")
                if total is not None:
                    score["mapped_delta_g"] = args.slope * total + args.intercept
            best = scoring.get("best_pose")
            if isinstance(best, dict):
                total = best.get("total_score")
                if total is not None:
                    best["mapped_delta_g"] = args.slope * total + args.intercept
                best.setdefault("delta_g_mapping", mapping_info)
            scoring.setdefault("delta_g_mapping", mapping_info)

        def _annotate_container(container: Any) -> None:
            if isinstance(container, dict):
                if "scoring" in container:
                    _apply_mapping(container.get("scoring"))
                if "result" in container and isinstance(container["result"], dict):
                    _annotate_container(container["result"])
                container.setdefault("delta_g_mapping", mapping_info)

        if isinstance(raw_data, list):
            for item in raw_data:
                _annotate_container(item)
        else:
            _annotate_container(raw_data)

        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(raw_data, indent=2, default=str))
        console.print(
            Panel.fit(
                f"Annotated ΔG mapping (slope={args.slope:.6g}, intercept={args.intercept:.6g})",
                title="ΔG Mapping",
            )
        )
        return

    if args.command == "mmgbsa":
        weight_data = json.loads(args.input.read_text())
        labels_filter = set(args.labels) if args.labels else None

        project_root = Path(__file__).resolve().parents[3]
        ff_files = ["amber14-all.xml", "amber14/tip3pfb.xml"]
        heme_ff = project_root / "forcefields" / "shahrokh_heme_ic6.ffxml"
        ligand_ff = project_root / "forcefields" / "UNL_gaff211.ffxml"
        if heme_ff.exists():
            ff_files.append(str(heme_ff))
        if ligand_ff.exists():
            ff_files.append(str(ligand_ff))

        gbsa_model = None if args.gbsa_model == "none" else args.gbsa_model
        calculator = MMGBSACalculator(ff_files, gbsa_model=gbsa_model)

        if args.relax_output:
            args.relax_output.mkdir(parents=True, exist_ok=True)

        results: List[Dict[str, Any]] = []
        for entry in weight_data:
            label = entry.get("label", "unknown")
            if labels_filter and label not in labels_filter:
                continue
            entry_data = entry.get("result", entry)
            poses = entry_data.get("poses", [])
            if not poses:
                console.print(f"[yellow]No poses found for entry {label}; skipping[/yellow]")
                continue
            pose_info = poses[0]
            pdb_path = pose_info.get("pdb_file")
            ligand_indices = pose_info.get("ligand_indices")
            if not pdb_path or not ligand_indices:
                console.print(f"[yellow]Missing pose PDB or ligand indices for entry {label}; skipping[/yellow]")
                continue

            try:
                modellers = None
                relaxation_records = None
                if args.relax:
                    entry_output_dir = None
                    if args.relax_output:
                        entry_output_dir = args.relax_output / label
                        entry_output_dir.mkdir(parents=True, exist_ok=True)

                    relaxation_records, modellers = calculator.relax_states(
                        pdb_path,
                        ligand_indices,
                        minimize=args.minimize,
                        md_steps=args.relax_md_steps,
                        random_seed=args.relax_random_seed,
                        output_dir=entry_output_dir,
                    )

                if args.md_steps > 0:
                    ensemble = calculator.calculate_ensemble(
                        pdb_path,
                        ligand_indices,
                        minimize=args.minimize,
                        md_steps=args.md_steps,
                        sample_interval=args.sample_interval,
                        modellers=modellers,
                    )
                    energies = {
                        "ensemble": ensemble.to_kcal_summary(),
                        "samples": len(ensemble.delta_g_series),
                    }
                else:
                    mmgbsa_result = calculator.calculate(
                        pdb_path,
                        ligand_indices,
                        minimize=args.minimize,
                        modellers=modellers,
                    )
                    energies = {
                        "single_frame": mmgbsa_result.to_kcal_per_mol(),
                    }
            except Exception as exc:
                console.print(f"[red]MM/GBSA evaluation failed for {label}: {exc}[/red]")
                continue

            results.append({
                "label": label,
                "pdb": pdb_path,
                "energies_kcal_per_mol": energies,
                "gbsa_model": gbsa_model or "none",
                "minimized": args.minimize,
                "md_steps": args.md_steps,
                "sample_interval": args.sample_interval,
            })
            if relaxation_records is not None:
                results[-1]["relaxation"] = relaxation_records

        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(results, indent=2))
        console.print(
            Panel.fit(
                f"Computed MM/GBSA energies for {len(results)} entries", title="MM/GBSA"
            )
        )
        return

    if args.command == "run":
        run_cli_simple(args.config)


def run_cli_simple(config_path: Path) -> None:
    pipeline = Pipeline(config_path)
    pipeline.run()


# Initialize console
console = Console()

def _looks_like_uniprot(s: str) -> bool:
    """Check if a string looks like a UniProt ID."""
    # Basic check for UniProt format (e.g., P12345, A0A0A0)
    uniprot_pattern = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$')
    return bool(uniprot_pattern.match(s.upper()))

def _looks_like_pdb(s: str) -> bool:
    """Check if a string looks like a PDB ID."""
    # PDB IDs are 4 alphanumeric characters
    return len(s) == 4 and s.isalnum()

def _select_from_list(items: List[Any], prompt: str, formatter=None) -> Any:
    """Present a list of items for selection."""
    if not items:
        return None
        
    if len(items) == 1:
        return items[0]
    
    # Format items for display
    if formatter:
        choices = [formatter(i, idx) for idx, i in enumerate(items, 1)]
    else:
        choices = [f"{i+1}. {str(item)}" for i, item in enumerate(items)]
    
    # Add a cancel option
    choices.append(questionary.Choice("Cancel", value=None))
    
    selected = questionary.select(
        prompt,
        choices=choices,
        use_shortcuts=True,
    ).ask()
    
    if selected is None:
        return None
    
    # Extract the index from the formatted string if needed
    if formatter:
        return selected
    
    try:
        idx = int(selected.split('.')[0]) - 1
        return items[idx]
    except (ValueError, IndexError):
        return None

def _format_protein_choice(protein: Dict, idx: int) -> str:
    """Format a protein choice for display."""
    name = protein.get('name', 'Unknown')
    uniprot_id = protein.get('uniprot_id', 'N/A')
    organism = protein.get('organism', 'N/A')
    return f"{idx}. {name} ({uniprot_id}) - {organism}"

def _format_structure_choice(structure: Dict, idx: int) -> str:
    """Format a structure choice for display."""
    pdb_id = structure.get('id', 'UNKNOWN')
    resolution = structure.get('resolution', 'N/A')
    method = structure.get('experimental_method', 'N/A')
    return f"{idx}. {pdb_id} - {method} - {resolution}Å"

def _is_valid_uniprot_id(uniprot_id: str) -> bool:
    """Check if the given string is a valid UniProt ID."""
    import re
    uniprot_pattern = re.compile(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$')
    return bool(uniprot_pattern.match(uniprot_id.upper()))

def run_interactive_analysis():
    """Run the interactive analysis workflow."""
    console.print(Panel.fit(
        "[bold blue]OVM-PK: Interactive Protein-Ligand Analysis[/bold blue]\n"
        "This tool will guide you through analyzing protein-ligand interactions.\n"
        "You'll need to provide a protein and optionally a ligand of interest.",
        title="Welcome"
    ))
    
    try:
        # Get protein information
        protein_query = questionary.text(
            "Enter protein name, UniProt ID, or PDB ID:",
            validate=lambda x: len(x.strip()) > 0
        ).ask()
        
        if not protein_query:
            console.print("[yellow]Protein query cannot be empty. Exiting.[/yellow]")
            return
        
        # Create output directory
        output_dir = Path("ovmpk_results") / protein_query.lower()
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Check if the input is a PDB ID
        if _looks_like_pdb(protein_query):
            pdb_id = protein_query.upper()
            console.print(f"[info] Treating input as PDB ID: {pdb_id}")
            pdb_path = output_dir / f"{pdb_id}.pdb"
            
            # Download the PDB file
            with console.status(f"[bold green]Downloading {pdb_id}...[/bold green]"):
                pdb_path = fetch_protein_structure(pdb_id, output_dir, {"fetch": {"protein": {"results_limit": 1}}})
                
            if not pdb_path or not pdb_path.exists():
                console.print(f"[red]Failed to download PDB {pdb_id}[/red]")
                return
                
            protein_name = f"PDB_{pdb_id}"
            
        # Check if the input is a UniProt ID
        elif _looks_like_uniprot(protein_query):
            uniprot_id = protein_query.upper()
            console.print(f"[info] Treating input as UniProt ID: {uniprot_id}")
            
            # Fetch the protein structure
            with console.status(f"[bold green]Fetching structure for {uniprot_id}...[/bold green]"):
                pdb_path = fetch_protein_structure(uniprot_id, output_dir)
                
            if not pdb_path or not pdb_path.exists():
                console.print(f"[red]Failed to fetch structure for {uniprot_id}[/red]")
                return
                
            protein_name = f"UniProt_{uniprot_id}"
            
        # Otherwise, treat as a protein name
        else:
            # Search UniProt for the protein name
            console.print(f"[info]Searching UniProt for: {protein_query}[/info]")
            try:
                from Bio import Entrez, SwissProt
                from io import StringIO
                
                # Set up Entrez email (required by NCBI)
                Entrez.email = "ovmpk@example.com"  # Replace with a valid email
                
                # Try different search strategies, with more specific terms for CYP3A4
                search_terms = [
                    f"{protein_query} AND reviewed:yes AND organism_id:9606",  # Human, reviewed
                    f"cyp3a4 AND cytochrome p450 3a4 AND human",  # Common name and full name
                    f"cyp3a4 AND p450 AND human",  # More general
                    f"{protein_query}* AND organism_id:9606",  # Wildcard with human filter
                    f"{protein_query} AND p450"  # General P450 search
                ]
                
                record = None
                for term in search_terms:
                    try:
                        handle = Entrez.esearch(db="protein", term=term, retmax=10)
                        record = Entrez.read(handle)
                        handle.close()
                        if record["IdList"]:
                            break
                    except Exception as e:
                        console.print(f"[yellow]Search with term '{term}' failed: {e}[/yellow]")
                        continue
                
                if not record or not record["IdList"]:
                    console.print(f"[red]No results found for {protein_query} in UniProt[/red]")
                    console.print("[yellow]Trying to find similar proteins...[/yellow]")
                    
                    # Try a more general search
                    try:
                        handle = Entrez.espell(term=protein_query, db="protein")
                        spell_result = Entrez.read(handle)
                        handle.close()
                        
                        if spell_result["Query"] != spell_result["ReplacedQuery"]:
                            suggestion = spell_result["ReplacedQuery"]
                            console.print(f"[info]Did you mean: {suggestion}?[/info]")
                            if questionary.confirm("Try with this suggestion?").ask():
                                handle = Entrez.esearch(db="protein", term=f"{suggestion} AND reviewed:yes", retmax=10)
                                record = Entrez.read(handle)
                                handle.close()
                    except Exception as e:
                        console.print(f"[yellow]Could not get spelling suggestions: {e}[/yellow]")
                
                if not record or not record["IdList"]:
                    console.print("[red]No matching proteins found. Please try a different search term.[/red]")
                    return
                
                # Get the best matching UniProt ID
                uniprot_id = None
                
                # First try to find a reviewed human protein
                for prot_id in record["IdList"]:
                    try:
                        # Fetch the record
                        handle = Entrez.efetch(db="protein", id=prot_id, rettype="gb", retmode="text")
                        record_text = handle.read()
                        handle.close()
                        
                        # Extract info
                        lines = record_text.split('\n')
                        
                        # Check if it's a reviewed (Swiss-Prot) entry
                        is_reviewed = any('dbxref="UniProt/Swiss-Prot"' in line for line in lines)
                        is_human = any('/organism="Homo sapiens"' in line for line in lines)
                        
                        # Extract UniProt ID
                        for line in lines:
                            if line.startswith('VERSION'):
                                uniprot_id = line.split()[1].split('.')[0]
                                break
                        
                        if is_reviewed and is_human and uniprot_id:
                            console.print(f"[info]Found reviewed human UniProt entry: {uniprot_id}[/info]")
                            break
                            
                    except Exception as e:
                        console.print(f"[yellow]Error processing protein record: {e}[/yellow]")
                        continue
                
                if not uniprot_id and record["IdList"]:
                    # Fallback to first result if no reviewed human protein found
                    uniprot_id = record["IdList"][0]
                    console.print(f"[yellow]Using first available UniProt ID: {uniprot_id}[/yellow]")
                elif not uniprot_id:
                    console.print("[yellow]No UniProt ID found. Will use protein name for PDB search.[/yellow]")
                
                # Configure search parameters
                search_config = {
                    "fetch": {
                        "protein": {
                            "results_limit": 5,
                            "max_resolution": 3.5,
                            "exp_method": "X-RAY DIFFRACTION",
                            "sort_by": "resolution_asc",
                            "target_species_id": 9606,  # Human
                            "required_ligands": []
                        }
                    }
                }

                # Try to find PDB structures
                with console.status("[bold green]Searching for PDB structures...[/bold green]"):
                    try:
                        # First try with protein name
                        console.print(f"[info]Searching RCSB for structures matching: {protein_query}[/info]")
                        pdb_path = fetch_protein_structure(protein_query, output_dir, search_config)
                        
                        # If no structure found, try with UniProt ID if available
                        if (not pdb_path or not pdb_path.exists()) and uniprot_id and _is_valid_uniprot_id(uniprot_id):
                            console.print(f"[info]No structures found with protein name. Trying with UniProt ID: {uniprot_id}[/info]")
                            pdb_path = fetch_protein_structure(uniprot_id, output_dir, search_config)
                    except Exception as e:
                        console.print(f"[red]Error fetching structure: {e}[/red]")
                        return
                
                if not pdb_path or not pdb_path.exists():
                    console.print("[red]Failed to fetch protein structure. Please try a different search term or PDB ID.[/red]")
                    console.print("[yellow]You can try searching with a PDB ID directly, e.g. '1TQN' for CYP3A4 with ketoconazole[/yellow]")
                    return
                
                protein_name = f"{protein_query}_{uniprot_id}"
                
            except ImportError:
                console.print("[yellow]Biopython is not installed. Using placeholder PDB ID.[/yellow]")
                console.print("Please install Biopython for better protein name resolution: pip install biopython")
                pdb_id = "1ABC"  # Fallback to a placeholder
                console.print(f"Using placeholder PDB ID: {pdb_id}")
                pdb_path = output_dir / f"{pdb_id}.pdb"
                protein_name = f"{protein_query}_1ABC"
                
                # Try to fetch the structure with the placeholder ID
                with console.status(f"[bold green]Fetching structure for PDB ID: {pdb_id}...[/bold green]"):
                    pdb_path = fetch_protein_structure(pdb_id, output_dir, {"fetch": {"protein": {"results_limit": 1}}})
                
                if not pdb_path or not pdb_path.exists():
                    console.print(f"[red]Failed to fetch structure for PDB ID: {pdb_id}[/red]")
                    return
        
        console.print(f"\n[bold]Selected Protein:[/bold] {protein_name}")
        console.print(f"[bold]Structure:[/bold] {pdb_path.name}")
        
        # Get ligand information
        ligand_name = questionary.text(
            "Enter ligand name or SMILES (leave empty for apo protein):"
        ).ask().strip()
        
        ligand_path = None
        if ligand_name:
            # Check if the input looks like a SMILES string
            if "=" in ligand_name or "(" in ligand_name:
                # Looks like a SMILES string
                smiles = ligand_name
                ligand_name = "LIG"  # Default name for SMILES
            else:
                # Treat as a name
                smiles = None
                
            with console.status("[bold green]Fetching ligand structure...[/bold green]"):
                try:
                    if smiles:
                        # Use SMILES directly
                        ligand_path = fetch_ligand_structure(
                            smiles, "smiles", output_dir
                        )
                    else:
                        # Try to find by name
                        ligand_path = fetch_ligand_structure(
                            ligand_name, "name", output_dir
                        )
                    
                    if not ligand_path or not ligand_path.exists():
                        console.print(f"[yellow]Could not find structure for {ligand_name}.[/yellow]")
                        if not questionary.confirm("Continue with apo protein?").ask():
                            return
                        ligand_name = None
                except Exception as e:
                    console.print(f"[yellow]Error fetching ligand: {e}[/yellow]")
                    if not questionary.confirm("Continue with apo protein?").ask():
                        return
                    ligand_name = None
        
        # Create configuration
        config = create_config(pdb_path, ligand_path, output_dir)
        
        # Save configuration
        config_path = output_dir / "config.yaml"
        with open(config_path, "w") as f:
            import yaml
            yaml.dump(config, f, default_flow_style=False)
            
        console.print(f"\n[bold]Configuration saved to:[/bold] {config_path}")
        
        # Confirm before running
        if not questionary.confirm("\nStart the analysis? This may take a while.").ask():
            console.print("[yellow]Analysis cancelled by user.[/yellow]")
            return
        
        # Run the analysis
        with Progress(
            SpinnerColumn(),
            "•",
            "[progress.description]{task.description}",
            "•",
            "[progress.percentage]{task.percentage:>3.0f}%"
        ) as progress:
            task = progress.add_task("[green]Running analysis...", total=100)
            
            # Run the pipeline
            try:
                results = run_analysis_integrated(config)
                progress.update(task, completed=100)
                
                # Generate reports
                if ligand_name and results.get("success", False):
                    console.print("\n[bold green]Analysis complete![/bold green]")
                    console.print(f"[bold]Results saved to:[/bold] {output_dir.absolute()}")
                    
                    # Display key results
                    if "binding_energy" in results:
                        dg = results["binding_energy"]
                        console.print(f"[bold]Binding Free Energy (ΔG):[/bold] {dg:.2f} kcal/mol")
                    
                    if "interactions" in results:
                        console.print("\n[bold]Key Interactions:[/bold]")
                        for i, interaction in enumerate(results["interactions"][:5], 1):
                            console.print(f"{i}. {interaction}")
                    
            except Exception as e:
                progress.stop()
                console.print(f"\n[red]Error during analysis: {str(e)}[/red]")
                if hasattr(e, '__traceback__'):
                    import traceback
                    console.print("\n[red]Traceback:[/red]")
                    console.print(traceback.format_exc())
                return

    except KeyboardInterrupt:
        console.print("\n[red]Analysis was cancelled by the user.[/red]")
    except Exception as e:
        console.print(f"\n[red]An error occurred: {str(e)}[/red]")
        if 'output_dir' in locals():
            console.print(f"Check the logs in {output_dir / 'ovmpk.log'} for more details.")

def main():
    """Entry point for running the analysis workflow."""
    import sys

    if len(sys.argv) > 1 and sys.argv[1] in {"chembl", "scan-weights", "run", "convert-dg", "map-dg", "mmgbsa"}:
        parser = build_parser()
        args = parser.parse_args()
        run_cli(args)
        return 0

    parser = argparse.ArgumentParser(description="OVM-PK analysis runner")
    parser.add_argument(
        "-c",
        "--config",
        type=str,
        help="Path to a YAML configuration file for non-interactive execution",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=str,
        help="Optional output directory override for the pipeline",
    )
    args = parser.parse_args()

    if args.config:
        config_path = Path(args.config).expanduser().resolve()
        if not config_path.exists():
            console.print(f"[red]Config file not found: {config_path}[/red]")
            return 1

        try:
            pipeline = Pipeline(config_path=config_path, output_dir=args.output_dir)
            pipeline.run()
            console.print("[green]Pipeline completed successfully.[/green]")
            return 0
        except Exception as e:
            console.print(f"[red]Fatal error during pipeline execution: {e}[/red]")
            import traceback
            console.print("\nTraceback:")
            console.print(traceback.format_exc())
            return 1

    try:
        run_interactive_analysis()
    except Exception as e:
        console.print(f"[red]Fatal error: {str(e)}[/red]")
        import traceback
        console.print("\nTraceback:")
        console.print(traceback.format_exc())
        return 1
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
