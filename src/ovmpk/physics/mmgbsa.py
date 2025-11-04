"""Simplified MM/GBSA calculator for single-snapshot binding free energies."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Sequence, Tuple

import numpy as np

from openmm import unit
from openmm import app
from openmm.app import ForceField, Modeller, PDBFile, Simulation
from openmm import LangevinIntegrator, LocalEnergyMinimizer


@dataclass
class MMGBSAResult:
    """Energies for complex, receptor, and ligand."""

    complex_energy: unit.Quantity
    receptor_energy: unit.Quantity
    ligand_energy: unit.Quantity

    @property
    def delta_g_bind(self) -> unit.Quantity:
        """Return ΔG_bind = E_complex − E_receptor − E_ligand."""
        return self.complex_energy - self.receptor_energy - self.ligand_energy

    def to_kcal_per_mol(self) -> dict:
        """Convenience conversion to kcal/mol scalars."""
        factor = unit.kilojoule_per_mole.conversion_factor_to(unit.kilocalorie_per_mole)
        return {
            "E_complex": self.complex_energy.value_in_unit(unit.kilojoule_per_mole) * factor,
            "E_receptor": self.receptor_energy.value_in_unit(unit.kilojoule_per_mole) * factor,
            "E_ligand": self.ligand_energy.value_in_unit(unit.kilojoule_per_mole) * factor,
            "deltaG_bind": self.delta_g_bind.value_in_unit(unit.kilojoule_per_mole) * factor,
        }


@dataclass
class MMGBSAEnsembleResult:
    """Ensemble of MM/GBSA energies with summary statistics."""

    complex_energies: np.ndarray  # kJ/mol
    receptor_energies: np.ndarray  # kJ/mol
    ligand_energies: np.ndarray  # kJ/mol

    @property
    def delta_g_series(self) -> np.ndarray:
        return self.complex_energies - self.receptor_energies - self.ligand_energies

    @staticmethod
    def _stats(values: np.ndarray) -> Dict[str, float]:
        return {
            "mean": float(values.mean()),
            "std": float(values.std(ddof=0)),
            "min": float(values.min()),
            "max": float(values.max()),
        }

    def to_kcal_summary(self) -> Dict[str, Dict[str, float]]:
        factor = unit.kilojoule_per_mole.conversion_factor_to(unit.kilocalorie_per_mole)

        def convert(stats: Dict[str, float]) -> Dict[str, float]:
            return {key: value * factor for key, value in stats.items()}

        return {
            "E_complex": convert(self._stats(self.complex_energies)),
            "E_receptor": convert(self._stats(self.receptor_energies)),
            "E_ligand": convert(self._stats(self.ligand_energies)),
            "deltaG_bind": convert(self._stats(self.delta_g_series)),
        }


class MMGBSACalculator:
    """Single-frame MM/GBSA style energy estimator."""

    def __init__(
        self,
        forcefield_files: Sequence[str],
        *,
        gbsa_model: Optional[str] = None,
        temperature: unit.Quantity = 300 * unit.kelvin,
        friction: unit.Quantity = 1 / unit.picosecond,
        minimization_tolerance_kj: float = 1.0,
        minimization_max_iterations: int = 200,
    ) -> None:
        self.forcefield = ForceField(*forcefield_files)
        self.gbsa_model = gbsa_model
        self.temperature = temperature
        self.friction = friction
        self.minimization_tolerance_kj = minimization_tolerance_kj
        self.minimization_max_iterations = minimization_max_iterations

    def _create_system(self, topology: app.Topology) -> app.System:
        kwargs = {"nonbondedMethod": app.NoCutoff}
        if self.gbsa_model:
            gb_models = {
                "GBn2": app.GBn2,
                "GBn": app.GBn,
                "OBC1": app.HCT,
            }
            kwargs["implicitSolvent"] = gb_models.get(self.gbsa_model, app.GBn2)
        try:
            return self.forcefield.createSystem(topology, **kwargs)
        except ValueError as exc:
            if "implicitSolvent" in str(exc) and "never used" in str(exc) and "implicitSolvent" in kwargs:
                # Custom force-fields (e.g., GAFF/porphyrin) may not supply GB parameters; retry in vacuum.
                fallback_kwargs = {k: v for k, v in kwargs.items() if k != "implicitSolvent"}
                return self.forcefield.createSystem(topology, **fallback_kwargs)
            raise

    def _potential_energy(
        self,
        topology: app.Topology,
        positions,
        *,
        minimize: bool = False,
    ) -> unit.Quantity:
        system = self._create_system(topology)
        integrator = LangevinIntegrator(
            self.temperature,
            self.friction,
            0.002 * unit.picoseconds,
        )
        simulation = Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        if minimize:
            LocalEnergyMinimizer.minimize(
                simulation.context,
                self.minimization_tolerance_kj,
                self.minimization_max_iterations,
            )
        state = simulation.context.getState(getEnergy=True)
        potential = state.getPotentialEnergy()
        del simulation, integrator, system
        return potential

    def _build_modellers(
        self,
        complex_pdb: str,
        ligand_indices: Sequence[int],
    ) -> Tuple[Modeller, Modeller, Modeller]:
        pdb = PDBFile(complex_pdb)
        complex_modeller = Modeller(pdb.topology, pdb.positions)

        ligand_modeller = self._subset_modeller(complex_modeller, ligand_indices)
        protein_indices = [atom.index for atom in complex_modeller.topology.atoms() if atom.index not in ligand_indices]
        receptor_modeller = self._subset_modeller(complex_modeller, protein_indices)

        return complex_modeller, receptor_modeller, ligand_modeller

    def _single_frame_result(
        self,
        modellers: Tuple[Modeller, Modeller, Modeller],
        *,
        minimize: bool,
    ) -> MMGBSAResult:
        complex_modeller, receptor_modeller, ligand_modeller = modellers
        complex_energy = self._potential_energy(
            complex_modeller.topology,
            complex_modeller.positions,
            minimize=minimize,
        )
        receptor_energy = self._potential_energy(
            receptor_modeller.topology,
            receptor_modeller.positions,
            minimize=minimize,
        )
        ligand_energy = self._potential_energy(
            ligand_modeller.topology,
            ligand_modeller.positions,
            minimize=minimize,
        )

        return MMGBSAResult(
            complex_energy=complex_energy,
            receptor_energy=receptor_energy,
            ligand_energy=ligand_energy,
        )

    def _ensemble_energies(
        self,
        topology: app.Topology,
        positions,
        *,
        minimize: bool,
        md_steps: int,
        sample_interval: int,
        random_seed: Optional[int],
    ) -> np.ndarray:
        system = self._create_system(topology)
        integrator = LangevinIntegrator(
            self.temperature,
            self.friction,
            0.002 * unit.picoseconds,
        )
        if random_seed is not None:
            integrator.setRandomNumberSeed(random_seed)
        simulation = Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        if minimize:
            LocalEnergyMinimizer.minimize(
                simulation.context,
                self.minimization_tolerance_kj,
                self.minimization_max_iterations,
            )

        energies: list[float] = []
        if md_steps <= 0:
            state = simulation.context.getState(getEnergy=True)
            energies.append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))
        else:
            simulation.context.setVelocitiesToTemperature(
                self.temperature,
                random_seed if random_seed is not None else 0,
            )
            total_steps = 0
            while total_steps < md_steps:
                steps = min(sample_interval, md_steps - total_steps)
                simulation.step(steps)
                total_steps += steps
                state = simulation.context.getState(getEnergy=True)
                energies.append(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole))

        del simulation, integrator, system
        return np.asarray(energies, dtype=float)

    def _relax_structure(
        self,
        topology: app.Topology,
        positions,
        *,
        minimize: bool,
        md_steps: int,
        random_seed: Optional[int],
    ) -> Tuple[unit.Quantity, unit.Quantity]:
        system = self._create_system(topology)
        integrator = LangevinIntegrator(
            self.temperature,
            self.friction,
            0.002 * unit.picoseconds,
        )
        if random_seed is not None:
            integrator.setRandomNumberSeed(random_seed)
        simulation = Simulation(topology, system, integrator)
        simulation.context.setPositions(positions)
        if minimize:
            LocalEnergyMinimizer.minimize(
                simulation.context,
                self.minimization_tolerance_kj,
                self.minimization_max_iterations,
            )
        if md_steps > 0:
            simulation.context.setVelocitiesToTemperature(
                self.temperature,
                random_seed if random_seed is not None else 0,
            )
            simulation.step(md_steps)

        state = simulation.context.getState(getPositions=True, getEnergy=True)
        relaxed_positions = state.getPositions()
        potential = state.getPotentialEnergy()

        del simulation, integrator, system
        return relaxed_positions, potential

    @staticmethod
    def _subset_modeller(original: Modeller, atom_indices: Iterable[int]) -> Modeller:
        atoms_to_keep = {idx for idx in atom_indices}
        modeller = Modeller(original.topology, original.positions)
        modeller.delete([atom for atom in modeller.topology.atoms() if atom.index not in atoms_to_keep])
        return modeller

    def calculate(
        self,
        complex_pdb: str,
        ligand_indices: Sequence[int],
        *,
        minimize: bool = False,
        modellers: Optional[Tuple[Modeller, Modeller, Modeller]] = None,
    ) -> MMGBSAResult:
        if not ligand_indices:
            raise ValueError("Ligand indices must be provided for MM/GBSA calculation")

        if modellers is None:
            modellers = self._build_modellers(complex_pdb, ligand_indices)

        return self._single_frame_result(modellers, minimize=minimize)

    def calculate_ensemble(
        self,
        complex_pdb: str,
        ligand_indices: Sequence[int],
        *,
        minimize: bool = True,
        md_steps: int = 5000,
        sample_interval: int = 50,
        modellers: Optional[Tuple[Modeller, Modeller, Modeller]] = None,
        random_seed: Optional[int] = None,
    ) -> MMGBSAEnsembleResult:
        if not ligand_indices:
            raise ValueError("Ligand indices must be provided for MM/GBSA calculation")

        if md_steps < 0 or sample_interval <= 0:
            raise ValueError("md_steps must be >= 0 and sample_interval must be positive")

        if modellers is None:
            modellers = self._build_modellers(complex_pdb, ligand_indices)

        complex_modeller, receptor_modeller, ligand_modeller = modellers

        complex_series = self._ensemble_energies(
            complex_modeller.topology,
            complex_modeller.positions,
            minimize=minimize,
            md_steps=md_steps,
            sample_interval=sample_interval,
            random_seed=random_seed,
        )
        receptor_series = self._ensemble_energies(
            receptor_modeller.topology,
            receptor_modeller.positions,
            minimize=minimize,
            md_steps=md_steps,
            sample_interval=sample_interval,
            random_seed=None if random_seed is None else random_seed + 1,
        )
        ligand_series = self._ensemble_energies(
            ligand_modeller.topology,
            ligand_modeller.positions,
            minimize=minimize,
            md_steps=md_steps,
            sample_interval=sample_interval,
            random_seed=None if random_seed is None else random_seed + 2,
        )

        return MMGBSAEnsembleResult(
            complex_energies=complex_series,
            receptor_energies=receptor_series,
            ligand_energies=ligand_series,
        )

    def relax_states(
        self,
        complex_pdb: str,
        ligand_indices: Sequence[int],
        *,
        minimize: bool = True,
        md_steps: int = 0,
        random_seed: Optional[int] = None,
        output_dir: Optional[Path] = None,
    ) -> Tuple[Dict[str, Any], Tuple[Modeller, Modeller, Modeller]]:
        if not ligand_indices:
            raise ValueError("Ligand indices must be provided for MM/GBSA relaxation")

        modellers = self._build_modellers(complex_pdb, ligand_indices)
        complex_modeller, receptor_modeller, ligand_modeller = modellers

        if output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

        relaxation_records: Dict[str, Any] = {}
        for offset, (label, modeller) in enumerate(
            (
                ("complex", complex_modeller),
                ("receptor", receptor_modeller),
                ("ligand", ligand_modeller),
            )
        ):
            relaxed_positions, potential = self._relax_structure(
                modeller.topology,
                modeller.positions,
                minimize=minimize,
                md_steps=md_steps,
                random_seed=None if random_seed is None else random_seed + offset,
            )
            modeller.positions = relaxed_positions
            record = {
                "energy_kj_per_mol": float(potential.value_in_unit(unit.kilojoule_per_mole)),
                "energy_kcal_per_mol": float(
                    potential.value_in_unit(unit.kilocalorie_per_mole)
                ),
            }
            if output_dir is not None:
                file_path = output_dir / f"{label}_relaxed.pdb"
                with open(file_path, "w") as handle:
                    PDBFile.writeFile(modeller.topology, relaxed_positions, handle)
                record["pdb"] = str(file_path)
            relaxation_records[label] = record

        return relaxation_records, modellers
