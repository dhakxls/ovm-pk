"""Stage 4 metal parameterization utilities."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Tuple

import json

import numpy as np
from openmm import unit
from openmm.app import PDBFile


@dataclass
class CoordinationBond:
    """Metal–ligand bond definition with validation metadata."""

    metal_label: str
    partner_label: str
    expected_length_nm: float
    tolerance_nm: float
    type1: str
    type2: str
    k_kj_per_mol_nm2: float
    observed_length_nm: float | None = None

    def within_tolerance(self) -> bool:
        if self.observed_length_nm is None:
            return False
        delta = abs(self.observed_length_nm - self.expected_length_nm)
        return delta <= self.tolerance_nm


@dataclass
class CoordinationReport:
    """Summary of a metal coordination environment."""

    metal_site: str
    bonds: List[CoordinationBond] = field(default_factory=list)
    bonded_xml_path: Path | None = None
    summary_path: Path | None = None

    @property
    def all_within_tolerance(self) -> bool:
        return all(bond.within_tolerance() for bond in self.bonds)


DEFAULT_PROFILE: Dict[str, Dict] = {
    "metal_residue": "HEM",
    "metal_atom": "FE",
    "bonds": [
        {
            "partner_residue": "HEM",
            "partner_atom": label,
            "expected_length_nm": 0.200,
            "tolerance_nm": 0.030,
            "type1": "fe",
            "type2": "nc" if label in {"NA", "NC"} else "nd",
            "k_kj_per_mol_nm2": 83680.0,
        }
        for label in ("NA", "NB", "NC", "ND")
    ]
    + [
        {
            "partner_residue": "CYP",
            "partner_residues": ["CYP", "CYS"],
            "partner_atom": "SG",
            "expected_length_nm": 0.266,
            "tolerance_nm": 0.040,
            "type1": "fe",
            "type2": "SH",
            "k_kj_per_mol_nm2": 66944.0,
        }
    ],
}


def _atom_lookup(pdb: PDBFile) -> Dict[Tuple[str, int, str], Tuple[int, unit.Quantity]]:
    mapping: Dict[Tuple[str, int, str], Tuple[int, unit.Quantity]] = {}
    for atom, pos in zip(pdb.topology.atoms(), pdb.positions):
        res = atom.residue
        mapping[(res.name.strip(), res.index, atom.name.strip())] = (atom.index, pos)
    return mapping


def _format_label(residue_name: str, atom_name: str, residue_index: int) -> str:
    return f"{residue_name}:{atom_name}@{residue_index}"


def _distance_nm(pos_a: unit.Quantity, pos_b: unit.Quantity) -> float:
    vec = pos_a - pos_b
    return float(np.linalg.norm(vec.value_in_unit(unit.nanometer)))


def identify_coordination(pdb_path: Path, profile: Dict[str, Dict] | None = None) -> CoordinationReport:
    if profile is None:
        profile = DEFAULT_PROFILE

    pdb = PDBFile(str(pdb_path))
    lookup = _atom_lookup(pdb)

    metal_candidates = [
        key for key in lookup if key[0] == profile["metal_residue"] and key[2] == profile["metal_atom"]
    ]
    if not metal_candidates:
        raise ValueError(
            f"No {profile['metal_residue']}:{profile['metal_atom']} metal centre found in {pdb_path}"
        )
    if len(metal_candidates) > 1:
        raise ValueError(
            f"Expected single metal centre, found {len(metal_candidates)} copies of"
            f" {profile['metal_residue']}:{profile['metal_atom']}"
        )

    metal_resname, metal_resid, metal_atom = metal_candidates[0]
    _, metal_pos = lookup[metal_candidates[0]]
    metal_label = _format_label(metal_resname, metal_atom, metal_resid)

    bonds: List[CoordinationBond] = []

    for entry in profile["bonds"]:
        partner_names = entry.get("partner_residues") or [entry["partner_residue"]]
        candidates = [
            key
            for key in lookup
            if key[0] in partner_names and key[2] == entry["partner_atom"]
        ]
        if not candidates:
            raise ValueError(
                f"Missing coordinating atom {entry['partner_residue']}:{entry['partner_atom']} in {pdb_path}"
            )

        distances = []
        for key in candidates:
            _, pos = lookup[key]
            distances.append((_distance_nm(metal_pos, pos), key))
        distance_nm, partner_key = min(distances, key=lambda x: x[0])
        partner_resname, partner_resid, partner_atom = partner_key

        bonds.append(
            CoordinationBond(
                metal_label=metal_label,
                partner_label=_format_label(partner_resname, partner_atom, partner_resid),
                expected_length_nm=entry["expected_length_nm"],
                tolerance_nm=entry["tolerance_nm"],
                type1=entry["type1"],
                type2=entry["type2"],
                k_kj_per_mol_nm2=entry["k_kj_per_mol_nm2"],
                observed_length_nm=distance_nm,
            )
        )

    return CoordinationReport(metal_site=metal_label, bonds=bonds)


def validate_geometry(report: CoordinationReport) -> None:
    failures = [bond for bond in report.bonds if not bond.within_tolerance()]
    if failures:
        details = ", ".join(
            f"{bond.metal_label}-{bond.partner_label}: observed {bond.observed_length_nm:.3f} nm"
            f" (expected {bond.expected_length_nm:.3f} ± {bond.tolerance_nm:.3f})"
            for bond in failures
        )
        raise ValueError(f"Metal coordination geometry outside tolerance: {details}")


def _serialize_bonded_xml(report: CoordinationReport) -> str:
    lines = ["<ForceField>", "  <HarmonicBondForce>"]
    for bond in report.bonds:
        lines.append(
            "    <Bond type1=\"{bond.type1}\" type2=\"{bond.type2}\" length=\"{length:.6f}\" k=\"{k:.6f}\"/>".format(
                bond=bond, length=bond.expected_length_nm, k=bond.k_kj_per_mol_nm2
            )
        )
    lines.extend(["  </HarmonicBondForce>", "</ForceField>"])
    return "\n".join(lines) + "\n"


def apply_metal_parameterization(
    pdb_path: Path,
    output_dir: Path,
    profile: Dict[str, Dict] | None = None,
) -> CoordinationReport:
    report = identify_coordination(pdb_path, profile)
    validate_geometry(report)

    output_dir.mkdir(parents=True, exist_ok=True)
    bonded_xml = output_dir / "metal_bonds.xml"
    bonded_xml.write_text(_serialize_bonded_xml(report), encoding="utf-8")
    report.bonded_xml_path = bonded_xml

    summary = {
        "metal_site": report.metal_site,
        "bonds": [
            {
                "partner": bond.partner_label,
                "expected_length_nm": bond.expected_length_nm,
                "observed_length_nm": bond.observed_length_nm,
                "tolerance_nm": bond.tolerance_nm,
            }
            for bond in report.bonds
        ],
    }
    summary_path = output_dir / "metal_coordination.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    report.summary_path = summary_path

    return report
