"""OVM-PK Command Line Interface.

This module serves as the main entry point for the OVM-PK command line tools.
"""
import sys
import click
from rich.console import Console

console = Console()

@click.group()
def cli():
    """OVM-PK: Physics-Aware Docking & Validation"""
    pass

@cli.command()
@click.option('--config', '-c', type=click.Path(exists=True), help='Path to configuration file')
@click.option('--interactive', '-i', is_flag=True, help='Run in interactive mode')
def run(config, interactive):
    """Run the OVM-PK pipeline."""
    if interactive or not config:
        from .run_analysis import run_interactive_analysis
        run_interactive_analysis()
    else:
        # Non-interactive mode with config file
        from ..pipeline_v2 import run_pipeline
        import json
        
        try:
            with open(config) as f:
                config_data = json.load(f)
            
            results = run_pipeline(config_data)
            console.print("\n[green]Analysis completed successfully![/green]")
            
        except Exception as e:
            console.print(f"[red]Error running pipeline: {str(e)}[/red]")
            sys.exit(1)

if __name__ == "__main__":
    cli()
