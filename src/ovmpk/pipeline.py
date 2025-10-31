"""Pipeline with ligand info in outputs."""
from pathlib import Path
from typing import Dict, Optional
import json
import time
import platform
from datetime import datetime
from .physics.registry import physics_registry
from .logging_config import configure_logging
from .reporting.html_report import generate_html_report
import logging

from .config.env_defaults import resolve_default_inputs

# Initialize logging
configure_logging()
logger = logging.getLogger(__name__)

class Pipeline:
    """Pipeline with comprehensive output reporting."""
    
    def __init__(self, config: Dict, work_dir: Optional[str] = None, output_dir: Optional[str] = None):
        logger.info("Initializing pipeline")
        self.config = config
        self.work_dir = Path(work_dir).resolve() if work_dir else Path("data/work").resolve()
        self.start_time = time.time()

        try:
            if 'metal_model' not in config.get('physics', {}):
                raise ValueError("Config missing required 'physics.metal_model' section")
            physics_registry.register_metal_handler(config['physics']['metal_model'])

            # Resolve protein/ligand defaults from configuration file
            resolved_inputs = resolve_default_inputs()
            system_cfg = self.config.setdefault("system", {})
            system_cfg.setdefault("protein_token", resolved_inputs.get("protein_token"))
            system_cfg.setdefault("pdb_id", resolved_inputs.get("pdb_id"))
            system_cfg.setdefault("resolver_preferences", resolved_inputs.get("resolver_preferences", {}))

            ligand_cfg = self.config.setdefault("ligand", {})
            ligand_cfg.setdefault("name", resolved_inputs.get("ligand"))

            # Handle output directory configuration
            reporting = config.get('reporting', {})
            base_output = Path(output_dir).resolve() if output_dir else self.work_dir.parent / "results"
            
            self.report_dir = Path(reporting.get('outdir', str(base_output / "reports"))).resolve()
            self.run_dir = Path(reporting.get('run_root', str(base_output / "runs"))).resolve() / datetime.now().strftime("%Y%m%d_%H%M%S")
            
            # Ensure parent directories exist
            self.report_dir.parent.mkdir(parents=True, exist_ok=True)
            self.run_dir.parent.mkdir(parents=True, exist_ok=True)
            
            logger.info(f"Output will be saved to {self.run_dir}")
            
        except Exception as e:
            logger.error(f"Pipeline initialization failed: {str(e)}")
            raise

    def run(self) -> Dict:
        """Execute pipeline with enhanced reporting."""
        try:
            logger.info("Starting pipeline execution")
            
            # Create directories
            self.work_dir.mkdir(parents=True, exist_ok=True)
            self.report_dir.mkdir(parents=True, exist_ok=True)
            self.run_dir.mkdir(parents=True, exist_ok=True)
            
            # Run pipeline
            pdb_id = self.config.get("system", {}).get("pdb_id")
            if not pdb_id:
                raise ValueError("System configuration missing resolved 'pdb_id'.")

            protein = self.work_dir / f"protein/{pdb_id}.pdb"
            protein.parent.mkdir(parents=True, exist_ok=True)
            
            exec_time = time.time() - self.start_time
            logger.info(f"Pipeline executed in {exec_time:.2f} seconds")
            
            results = {
                "system_info": self._get_system_info(),
                "performance": {
                    "execution_time_sec": exec_time,
                    "platform": platform.platform()
                },
                "config": {
                    "physics": self.config["physics"]["metal_model"],
                    "system": self.config.get("system", {}),
                    "ligand": self.config.get("ligand", {"name": "unknown"}),
                    "docking": self.config.get("docking", {})
                },
                "timestamp": datetime.now().isoformat(),
                "status": "success"
            }
            
            # Run free energy calculation
            system = self._get_system()
            free_energy = self._run_free_energy(system)
            results["free_energy"] = free_energy
            
            # Write outputs
            self._write_outputs(results)
            self.generate_report(results)
            logger.info(f"Results saved to {self.run_dir}")
            
            return results
            
        except Exception as e:
            logger.error(f"Pipeline execution failed: {str(e)}", exc_info=True)
            return {
                "status": "failed",
                "error": str(e),
                "timestamp": datetime.now().isoformat()
            }
    
    def _get_system_info(self) -> Dict:
        """Get system hardware/software info."""
        return {
            "python_version": platform.python_version(),
            "os": platform.system(),
            "cpu": platform.processor() or "Unknown",
            "openmm_version": self._get_openmm_version()
        }
    
    def _get_openmm_version(self) -> str:
        """Get OpenMM version if available."""
        try:
            import openmm
            return openmm.version.short_version
        except ImportError:
            return "Not available"
    
    def _write_outputs(self, results: Dict):
        """Serialize pipeline outputs including HTML report."""
        try:
            # Main results JSON
            results_json = self.run_dir / "results.json"
            with open(results_json, "w") as f:
                json.dump(results, f, indent=2)
            
            # Summary report
            with open(self.report_dir / "latest_run.json", "w") as f:
                json.dump({
                    "run_path": str(self.run_dir),
                    "status": results["status"],
                    "execution_time": results.get("performance", {}).get("execution_time_sec", 0),
                    "timestamp": results["timestamp"]
                }, f, indent=2)
            
            # HTML report
            html_report = self.run_dir / "report.html"
            try:
                generate_html_report(results_json, html_report)
                logger.info(f"Generated HTML report at {html_report}")
            except Exception as e:
                logger.warning(f"Failed to generate HTML report: {str(e)}")
                
        except Exception as e:
            logger.error(f"Failed to write outputs: {str(e)}")
            raise
    
    def _run_free_energy(self, system):
        """Perform alchemical free energy calculations"""
        from openmmtools.alchemy import AbsoluteAlchemicalFactory
        alchemical_system = AbsoluteAlchemicalFactory(system)
        # MBAR analysis implementation
        return deltaG
    
    def generate_report(self, results):
        """Create JSON report with key metrics"""
        import json
        report = {
            "timestamp": datetime.now().isoformat(),
            "deltaG": results["free_energy"],
            "predicted_Ki": self._calculate_ki(results["free_energy"]),
            "system": self.config.get("system", {}).get("pdb_id")
        }
        with open(self.run_dir/"report.json", "w") as f:
            json.dump(report, f)
    
    def _get_system(self):
        # TO DO: implement system retrieval
        pass
    
    def _calculate_ki(self, free_energy):
        # TO DO: implement Ki calculation
        pass
