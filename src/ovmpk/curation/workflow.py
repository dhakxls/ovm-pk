"""Core workflow for automated ligand/enzyme asset curation."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

from ..fetchers.core import FlexibleFetcher
from ..prep import ligand_prep, protein_prep
from ..analysis.metal_parameterization import apply_metal_parameterization, CoordinationReport
from .input_manager import AutoCurateConfig
from .ligand_parameters import build_forcefield_bundle, LigandParameterBundle

logger = logging.getLogger(__name__)


@dataclass
class AutoCurateResult:
    """Bundle of curated asset paths."""

    run_dir: Path
    protein_path: Path
    ligand_path: Path
    protein_prepped: Path
    ligand_prepped: Path
    system_xml: Optional[Path]
    solvated_pdb: Optional[Path]
    metal_report: Optional[CoordinationReport]
    ffxml_bundle: Optional[LigandParameterBundle]


class AutoCurator:
    """Coordinate automated data prep, parameterization, and packaging."""

    def __init__(self, config: AutoCurateConfig):
        self.config = config
        self.run_dir = config.run_dir
        self.stage_dirs = {
            "stage1": self.run_dir / "stage1",
            "protein_prep": self.run_dir / "protein_prep",
            "ligand_prep": self.run_dir / "ligand_prep",
            "system_setup": self.run_dir / "system_setup",
            "metal": self.run_dir / "metal_parameters",
            "assets": self.run_dir / "curated_assets",
        }
        self.ffxml_bundle: Optional[LigandParameterBundle] = None

    def run(self) -> AutoCurateResult:
        logger.info("Starting auto-curation for protein=%s ligand=%s", self.config.protein_input, self.config.ligand_input)
        self._prepare_directories()

        protein_raw, ligand_raw = self._stage_inputs()
        protein_prepped = self._prepare_protein(protein_raw)
        ligand_prepped = self._prepare_ligand(ligand_raw)

        self.ffxml_bundle = self._generate_forcefield_bundle(ligand_prepped)
        system_xml, solvated_pdb = self._setup_system(protein_prepped, ligand_prepped)
        metal_report = self._parameterize_metals(protein_prepped)

        logger.info("Auto-curation completed. Assets stored under %s", self.run_dir)
        return AutoCurateResult(
            run_dir=self.run_dir,
            protein_path=protein_raw,
            ligand_path=ligand_raw,
            protein_prepped=protein_prepped,
            ligand_prepped=ligand_prepped,
            system_xml=system_xml,
            solvated_pdb=solvated_pdb,
            metal_report=metal_report,
            ffxml_bundle=self.ffxml_bundle,
        )

    # ------------------------------------------------------------------
    # Stage helpers
    # ------------------------------------------------------------------
    def _prepare_directories(self) -> None:
        self.run_dir.mkdir(parents=True, exist_ok=True)
        for path in self.stage_dirs.values():
            path.mkdir(parents=True, exist_ok=True)

    def _stage_inputs(self) -> tuple[Path, Path]:
        fetch_config = self._build_fetch_config()
        fetcher = FlexibleFetcher(fetch_config)

        protein_path = fetcher.fetch_protein(self.config.protein_input)
        ligand_path = fetcher.fetch_ligand(self.config.ligand_input)

        logger.info("Stage 1 complete: protein=%s ligand=%s", protein_path, ligand_path)
        return protein_path, ligand_path

    def _build_fetch_config(self) -> Dict[str, Any]:
        stage1_dir = self.stage_dirs["stage1"]
        protein_cfg = {
            "save_dir": str(stage1_dir / "protein"),
            "resolver": self.config.fetch.get("protein_resolver", {}),
        }
        ligand_cfg: Dict[str, Any] = {
            "save_dir": str(stage1_dir / "ligand"),
            "identifier_type": self.config.fetch.get("ligand_identifier", "name"),
        }
        if self.config.fetch.get("ligand_smiles"):
            ligand_cfg["smiles"] = self.config.fetch["ligand_smiles"]
        return {"protein": protein_cfg, "ligand": ligand_cfg}

    def _prepare_protein(self, protein_path: Path) -> Path:
        base_cfg = {
            "ph": 7.4,
            "output_suffix": "_stage2_prepped",
            "run_pdbfixer": True,
            "output_dir": str(self.stage_dirs["protein_prep"]),
            "heterogens": self.config.prep.get("heterogens"),
        }
        user_overrides = self.config.prep.get("protein", {})
        base_cfg.update(user_overrides)
        prep_cfg = {"prep": {"protein": base_cfg}}

        result = protein_prep.prepare({"apo": protein_path}, prep_cfg)
        logger.info("Stage 2 protein prep complete -> %s", result)
        return Path(result)

    def _prepare_ligand(self, ligand_path: Path) -> Path:
        base_cfg: Dict[str, Any] = {
            "ph": 7.4,
            "charge_method": "gasteiger",
            "output_suffix": "_stage2_prepped",
            "force_3d": True,
            "use_dimorphite": False,
            "output_dir": str(self.stage_dirs["ligand_prep"]),
        }
        user_overrides = self.config.prep.get("ligand", {})
        base_cfg.update(user_overrides)
        prep_cfg = {"prep": {"ligand": base_cfg}}

        result = ligand_prep.prepare(ligand_path, prep_cfg)
        logger.info("Stage 2 ligand prep complete -> %s", result)
        return Path(result)

    def _generate_forcefield_bundle(self, ligand_prepped: Path) -> Optional[LigandParameterBundle]:
        try:
            bundle = build_forcefield_bundle(
                ligand_prepped,
                self.stage_dirs["assets"],
                self.config.ligand,
            )
        except Exception as exc:  # pragma: no cover - depends on filesystem/tooling
            logger.warning("Force-field bundling failed: %s", exc)
            return None

        logger.info("Generated ligand force-field bundle -> %s", bundle.ffxml)
        return bundle

    def _setup_system(self, protein_prepped: Path, ligand_prepped: Path) -> Tuple[Optional[Path], Optional[Path]]:
        try:
            from openmm import XmlSerializer, app, unit
            from rdkit import Chem
        except ImportError as exc:  # pragma: no cover - environment dependent
            logger.warning("Stage 3 skipped (missing dependency): %s", exc)
            return None, None

        ligand_supplier = Chem.SDMolSupplier(str(ligand_prepped), removeHs=False)
        ligand_mol = ligand_supplier[0] if ligand_supplier else None
        if ligand_mol is None:
            raise RuntimeError(f"Ligand could not be loaded from {ligand_prepped}")

        ligand_pdb_path = self.stage_dirs["ligand_prep"] / f"{ligand_prepped.stem}.pdb"
        Chem.MolToPDBFile(ligand_mol, str(ligand_pdb_path))

        protein_pdb = app.PDBFile(str(protein_prepped))
        ligand_pdb = app.PDBFile(str(ligand_pdb_path))

        modeller = app.Modeller(protein_pdb.topology, protein_pdb.positions)
        modeller.add(ligand_pdb.topology, ligand_pdb.positions)

        forcefield_files = self._resolve_forcefield_files()
        try:
            forcefield = app.ForceField(*forcefield_files)
        except Exception as exc:  # pragma: no cover - file availability
            raise RuntimeError(f"Failed to load force fields: {exc}") from exc

        modeller.addSolvent(
            forcefield,
            model="tip3p",
            padding=1.0 * unit.nanometer,
            ionicStrength=0.15 * unit.molar,
            neutralize=True,
        )

        system = forcefield.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            constraints=app.HBonds,
        )

        output_dir = self.stage_dirs["system_setup"]
        output_dir.mkdir(parents=True, exist_ok=True)

        solvated_pdb = output_dir / f"{protein_prepped.stem}_stage3_solvated.pdb"
        with open(solvated_pdb, "w", encoding="utf-8") as fh:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)

        system_xml = output_dir / f"{protein_prepped.stem}_stage3_system.xml"
        system_xml.write_text(XmlSerializer.serialize(system), encoding="utf-8")

        logger.info("Stage 3 system setup complete -> %s", system_xml)
        return system_xml, solvated_pdb

    def _resolve_forcefield_files(self) -> list[str]:
        custom = self.config.ligand.get("forcefield_files")
        if custom:
            return [str(Path(path)) for path in custom]
        if self.ffxml_bundle is not None:
            return [
                "amber14/protein.ff14SB.xml",
                "amber14/tip3p_standard.xml",
                "amber14/lipid17.xml",
                str(self.ffxml_bundle.ffxml),
            ]
        return [
            "amber14/protein.ff14SB.xml",
            "amber14/tip3p_standard.xml",
            "amber14/lipid17.xml",
            str(Path("forcefields/shahrokh_heme_ic6_unl.ffxml")),
        ]

    def _parameterize_metals(self, protein_prepped: Path) -> Optional[CoordinationReport]:
        if not self.config.metal.get("enable", True):
            logger.info("Stage 4 metal parameterization disabled by config")
            return None

        output_dir = self.stage_dirs["metal"]
        profile = self.config.metal.get("profile")
        report = apply_metal_parameterization(protein_prepped, output_dir, profile=profile)
        if report is None:
            logger.warning("Metal parameterization produced no report for %s", protein_prepped)
            return None

        logger.info("Stage 4 metal parameterization complete -> %s", report.bonded_xml_path)
        return report
