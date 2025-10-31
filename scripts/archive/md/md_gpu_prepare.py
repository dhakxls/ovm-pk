#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Prepare a solvated, GPU-ready OpenMM system from a receptor PDB and a docked ligand pose (SDF),
with stage-by-stage timers.

Key features
------------
- Reads a specific pose (1-based) from an SDF.
- Optionally reassigns correct bond orders to the docked pose using a trusted SMILES
  via RDKit's AssignBondOrdersFromTemplate (--ligand-smiles) or auto-resolve it (--auto-ligand-smiles).
- Robust fallback that de-radicalizes, re-sanitizes, and adds explicit hydrogens.
- Builds an OpenFF Molecule for ligand parameterization (SMIRNOFF 2.0.0).
- Optional heme ffxml loading for HEM/HEC systems when present, or strip heme by flag.
- Optional stripping of common crystallographic residues (e.g., GOL/EDO/SO4).
- Stage-wise workflow with --stop-after markers for debugging.
- Timers: per-stage durations printed, saved to durations.json, and embedded in prepare.json.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import sys
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import requests

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem, rdmolops
from rdkit.Chem import rdFMCS
from rdkit.Chem.rdmolops import GetMolFrags

# Optional standardization (MolVS-style); present in many RDKit builds
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as Std
    HAS_STD = True
except Exception:
    try:
        from rdkit.Chem import rdMolStandardize as Std
        HAS_STD = True
    except Exception:
        HAS_STD = False

# RDKit InChI availability (optional but helpful)
try:
    from rdkit.Chem.inchi import MolToInchiKey
    HAS_INCHI = True
except Exception:
    HAS_INCHI = False

# OpenMM / OpenFF
import openmm as mm
from openmm import Platform, XmlSerializer, app, unit
from openff.toolkit.topology import Molecule
from openff.toolkit.utils.exceptions import RadicalsNotSupportedError
from openmmforcefields.generators import SMIRNOFFTemplateGenerator

# PDBFixer is optional
try:
    import pdbfixer
    HAS_FIXER = True
except Exception:
    HAS_FIXER = False

# Membrane MD settings
from ..prep.membrane import is_membrane_system
from openmm import MonteCarloBarostat, Lipid14Force

STAGES = ["ligand", "receptor", "merge", "forcefield", "solvate", "system", "minimize", "all"]


# ---------- helpers ----------

def pick_platform(preferred: str | None = None) -> Platform:
    names = [Platform.getPlatform(i).getName() for i in range(Platform.getNumPlatforms())]
    order = ([preferred] if preferred else []) + ["CUDA", "OpenCL", "CPU"]
    for name in order:
        if name and name in names:
            return Platform.getPlatformByName(name)
    return Platform.getPlatformByName("CPU")


def mark_ok(outdir: Path, stage: str):
    (outdir / f"{stage}.ok").write_text("ok\n")
    print(f"STEP_OK: {stage}", flush=True)


def sha1sum(path: Path) -> str:
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write_rdkit_pdb(mol: Chem.Mol, out_pdb: Path):
    w = Chem.rdmolfiles.PDBWriter(str(out_pdb))
    try:
        w.write(mol)
    finally:
        w.close()


def try_build_forcefield(off_lig: Molecule, water_model_for_xml: str):
    """
    Amber14SB (protein), water (user-chosen model -> ffxml name), SMIRNOFF 2.0.0 for ligand.
    Tries to load a heme ffxml if available.
    """
    # Map modeller model name to ffxml filename for water
    water_xml = {
        "tip3p": "amber14/tip3p.xml",
        "tip3pfb": "amber14/tip3pfb.xml",
        "spce": "amber14/spce.xml",
        "tip4pew": "amber14/tip4pew.xml",
    }.get(water_model_for_xml.lower(), "amber14/tip3p.xml")

    smirnoff = SMIRNOFFTemplateGenerator(molecules=[off_lig], forcefield="openff-2.0.0.offxml")
    heme_ffxml_candidates = [
        "amber/hemhemea.xml",
        "amber/heme.xml",
        "openmmforcefields/ffxml/amber/hemhemea.xml",
        "openmmforcefields/ffxml/amber/heme.xml",
        "openmmforcefields/ffxml/heme.xml",
    ]
    base = ["amber14/protein.ff14SB.xml", water_xml]
    ff, heme_loaded = None, None
    for h in [None] + heme_ffxml_candidates:
        try:
            files = base if h is None else base + [h]
            ff = app.ForceField(*files)
            ff.registerTemplateGenerator(smirnoff.generator)
            heme_loaded = h
            break
        except Exception:
            continue
    if heme_loaded:
        print(f"[info] Heme parameters loaded from: {heme_loaded}", flush=True)
    else:
        print("[warn] No heme ffxml found; HEM/HEC may fail if present.", flush=True)
    return ff, heme_loaded


def _std_cleanup(mol: Chem.Mol) -> Chem.Mol:
    """Standardize/uncharged cleanup to help resolve odd valence/charge encodings."""
    if not HAS_STD:
        return mol
    try:
        mol = Std.Cleanup(mol)
        mol = Std.Normalizer().normalize(mol)
        mol = Std.Reionizer().reionize(mol)
        mol = Std.Uncharger().uncharge(mol)
    except Exception:
        pass
    return mol


# ---------- SMILES auto-resolution helpers ----------

def _largest_organic_fragment(m: Chem.Mol) -> Chem.Mol:
    try:
        frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
        if frags and len(frags) > 1:
            frags = sorted(frags, key=lambda x: sum(a.GetAtomicNum() > 1 for a in x.GetAtoms()), reverse=True)
            return Chem.Mol(frags[0])
    except Exception:
        pass
    return m


def _rdkit_canonical_smiles_for_openff(m: Chem.Mol) -> str | None:
    try:
        mH = Chem.AddHs(m, addCoords=True)
        Chem.SanitizeMol(mH, catchErrors=True)
        mH = _std_cleanup(mH)
        for a in mH.GetAtoms():
            if a.GetNumRadicalElectrons():
                a.SetNumRadicalElectrons(0)
        Chem.SanitizeMol(mH, catchErrors=True)
        s = Chem.MolToSmiles(Chem.RemoveHs(mH), isomericSmiles=True)
        Molecule.from_smiles(s, allow_undefined_stereo=True)
        return s
    except Exception:
        return None


def _sdf_embedded_smiles(sdf_path: Path, pose_index: int) -> str | None:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    i = max(0, pose_index - 1)
    mol = suppl[i] if 0 <= i < len(suppl) else None
    if mol is None:
        return None
    for k in mol.GetPropNames():
        if "SMILES" in k.upper():
            val = mol.GetProp(k).strip()
            try:
                Molecule.from_smiles(val, allow_undefined_stereo=True)
                return val
            except Exception:
                continue
    return _rdkit_canonical_smiles_for_openff(Chem.Mol(mol))


def _pubchem_smiles_by_inchikey(inchikey: str) -> str | None:
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/IsomericSMILES/JSON"
    r = requests.get(url, timeout=10)
    if r.ok:
        js = r.json()
        props = js.get("PropertyTable", {}).get("Properties", [])
        if props:
            return props[0].get("IsomericSMILES")
    return None


def _pubchem_smiles_by_name(name: str) -> str | None:
    # name → properties (prefer IsomericSMILES, fallback CanonicalSMILES)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{requests.utils.quote(name)}/property/IsomericSMILES,CanonicalSMILES/JSON"
    r = requests.get(url, timeout=10)
    if r.ok:
        js = r.json()
        props = js.get("PropertyTable", {}).get("Properties", [])
        if props:
            for key in ("IsomericSMILES", "CanonicalSMILES"):
                s = props[0].get(key)
                if s:
                    try:
                        Molecule.from_smiles(s, allow_undefined_stereo=True)
                        return s
                    except Exception:
                        pass
    return None


def _chembl_smiles_by_inchikey(inchikey: str) -> str | None:
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
    r = requests.get(url, params={"molecule_structures__standard_inchi_key": inchikey}, timeout=10)
    if r.ok:
        js = r.json()
        recs = js.get("molecules", [])
        if recs:
            return recs[0].get("molecule_structures", {}).get("canonical_smiles")
    return None


def _chembl_smiles_by_name(name: str) -> str | None:
    # fuzzy synonym match; keep small limit
    url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
    r = requests.get(url, params={"molecule_synonyms__icontains": name, "limit": 5}, timeout=10)
    if r.ok:
        js = r.json()
        for rec in js.get("molecules", []):
            s = (rec.get("molecule_structures") or {}).get("canonical_smiles")
            if s:
                try:
                    Molecule.from_smiles(s, allow_undefined_stereo=True)
                    return s
                except Exception:
                    continue
    return None


def _cactus_smiles_by_name(name: str) -> str | None:
    url = f"https://cactus.nci.nih.gov/chemical/structure/{name}/smiles"
    r = requests.get(url, timeout=10)
    if r.ok and r.text.strip():
        return r.text.strip()
    return None


def resolve_trusted_smiles(sdf_path: Path, pose_index: int = 1, ligand_name_hint: str | None = None) -> str:
    """
    Returns a SMILES that OpenFF can parse (raises on failure).
    Strategy:
      1) embedded SDF property → validate with OpenFF
      2) PubChem by InChIKey (computed from pose) → validate
      3) PubChem by name hint (filename or provided) → validate
      4) ChEMBL by InChIKey → validate
      5) ChEMBL by name → validate
      6) NIH Cactus by name → validate
      7) RDKit canonicalized (std/uncharged) → validate
    """
    # 1) from SDF properties / RDKit canonical
    s = _sdf_embedded_smiles(sdf_path, pose_index)
    if s:
        try:
            Molecule.from_smiles(s, allow_undefined_stereo=True)
            return s
        except Exception:
            pass

    # compute inchikey if RDKit has InChI enabled
    inchikey = None
    if HAS_INCHI:
        try:
            suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
            mol = suppl[max(0, pose_index - 1)]
            mol = _largest_organic_fragment(Chem.Mol(mol))
            mol = _std_cleanup(mol)
            inchikey = MolToInchiKey(Chem.RemoveHs(mol))
        except Exception:
            inchikey = None

    # 2) PubChem by InChIKey
    if inchikey:
        s = _pubchem_smiles_by_inchikey(inchikey)
        if s:
            Molecule.from_smiles(s, allow_undefined_stereo=True)
            return s

    # build name hint (filename often contains the ligand)
    if not ligand_name_hint:
        ligand_name_hint = Path(sdf_path).stem
        ligand_name_hint = re.sub(r"^ligand_best_pose_\d+_from_|_poses?$", "", ligand_name_hint)
        m = re.search(r"[A-Za-z][A-Za-z0-9\-]{2,}", ligand_name_hint)
        ligand_name_hint = m.group(0) if m else ligand_name_hint

    # 3) PubChem by name
    if ligand_name_hint:
        s = _pubchem_smiles_by_name(ligand_name_hint)
        if s:
            Molecule.from_smiles(s, allow_undefined_stereo=True)
            return s

    # 4) ChEMBL by InChIKey
    if inchikey:
        s = _chembl_smiles_by_inchikey(inchikey)
        if s:
            Molecule.from_smiles(s, allow_undefined_stereo=True)
            return s

    # 5) ChEMBL by name
    if ligand_name_hint:
        s = _chembl_smiles_by_name(ligand_name_hint)
        if s:
            Molecule.from_smiles(s, allow_undefined_stereo=True)
            return s

    # 6) NIH Cactus by name
    if ligand_name_hint:
        s = _cactus_smiles_by_name(ligand_name_hint)
        if s:
            Molecule.from_smiles(s, allow_undefined_stereo=True)
            return s

    # 7) offline RDKit canonical fallback (already validated inside)
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = _largest_organic_fragment(Chem.Mol(suppl[max(0, pose_index - 1)]))
    s = _rdkit_canonical_smiles_for_openff(mol)
    if s:
        return s

    raise RuntimeError("Failed to resolve an OpenFF-compatible SMILES via all strategies.")


def _deradicalize_and_add_hs(mol: Chem.Mol) -> Chem.Mol:
    Chem.SanitizeMol(mol, catchErrors=True)
    mol.UpdatePropertyCache(strict=False)
    if hasattr(rdmolops, "RemoveRadicals"):
        try:
            rdmolops.RemoveRadicals(mol)
        except Exception:
            pass
    for a in mol.GetAtoms():
        if a.GetNumRadicalElectrons():
            a.SetNumRadicalElectrons(0)
    Chem.SanitizeMol(mol, catchErrors=True)
    mol.UpdatePropertyCache(strict=False)
    molH = Chem.AddHs(mol, addCoords=True)
    try:
        Chem.Kekulize(molH, clearAromaticFlags=True)
    except Exception:
        pass
    Chem.SanitizeMol(molH, catchErrors=True)
    molH.UpdatePropertyCache(strict=False)
    return molH


def _strip_residues_inplace(top: app.Topology, pos, names_to_strip: set[str]) -> tuple[app.Topology, Any]:
    """Return a new Topology/positions with selected residue names removed."""
    modeller = app.Modeller(top, pos)
    to_delete = []
    for res in modeller.topology.residues():
        if res.name in names_to_strip:
            to_delete.append(res)
    if to_delete:
        modeller.delete(to_delete)
    return modeller.topology, modeller.positions


# ---------- timers ----------

@contextmanager
def stage_timer(stage: str, durations: dict, outdir: Path):
    t0 = time.perf_counter()
    try:
        yield
    finally:
        dt = time.perf_counter() - t0
        durations[stage] = dt
        print(f"[time] {stage}: {dt:.3f}s", flush=True)
        # write rolling durations so we see progress even on early exit
        (outdir / "durations.json").write_text(json.dumps(durations, indent=2))


# ---------- ligand loading & repair ----------

def load_ligand_from_sdf(sdf_path: Path, pose_index: int, ref_smiles: str | None = None):
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    idx0 = max(0, pose_index - 1)
    mol = suppl[idx0] if 0 <= idx0 < len(suppl) else None
    if mol is None:
        raise RuntimeError(f"Could not read pose {pose_index} from {sdf_path}")

    mol = Chem.Mol(mol)
    try:
        frags = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        if frags and len(frags) > 1:
            frags = sorted(frags, key=lambda m: sum(a.GetAtomicNum() > 1 for a in m.GetAtoms()), reverse=True)
            mol = Chem.Mol(frags[0])
            print(f"[info] Kept largest fragment: {mol.GetNumAtoms()} atoms", flush=True)
    except Exception:
        pass

    molH_pose = Chem.AddHs(mol, addCoords=True)
    Chem.SanitizeMol(molH_pose, catchErrors=True)
    molH_pose.UpdatePropertyCache(strict=False)

    # A) Preferred: OFF from trusted SMILES + coordinates from pose via MCS
    if ref_smiles:
        try:
            offmol = Molecule.from_smiles(ref_smiles, allow_undefined_stereo=True)
            off_rd = offmol.to_rdkit()

            poseH = Chem.AddHs(Chem.Mol(mol), addCoords=True)
            Chem.SanitizeMol(poseH, catchErrors=True)
            poseH.UpdatePropertyCache(strict=False)

            q = Chem.RemoveHs(off_rd)
            t = Chem.RemoveHs(poseH)
            mcs = rdFMCS.FindMCS(
                [q, t],
                ringMatchesRingOnly=True,
                completeRingsOnly=True,
                bondCompare=rdFMCS.BondCompare.CompareAny,
                atomCompare=rdFMCS.AtomCompare.CompareElements,
            )
            patt = Chem.MolFromSmarts(mcs.smartsString)
            q_idx = q.GetSubstructMatch(patt)
            t_idx = t.GetSubstructMatch(patt)
            if not q_idx or not t_idx or len(q_idx) != len(t_idx):
                raise RuntimeError("Could not map OFF(SMILES) to pose by MCS elements.")

            pose_conf = poseH.GetConformer()
            if off_rd.GetNumConformers() == 0:
                off_rd.AddConformer(Chem.Conformer(off_rd.GetNumAtoms()), assignId=True)
            off_conf = off_rd.GetConformer()
            for off_hvy, pose_hvy in zip(q_idx, t_idx):
                p = pose_conf.GetAtomPosition(pose_hvy)
                off_conf.SetAtomPosition(off_hvy, p)
            for atom in off_rd.GetAtoms():
                if atom.GetAtomicNum() == 1:
                    nbr = atom.GetNeighbors()[0]
                    p = off_conf.GetAtomPosition(nbr.GetIdx())
                    off_conf.SetAtomPosition(atom.GetIdx(), (p.x + 0.1, p.y, p.z))

            offmol = Molecule.from_rdkit(off_rd, allow_undefined_stereo=True, hydrogens_are_explicit=True)
            rdkit_for_pdb = offmol.to_rdkit()
            smiles_log = Chem.MolToSmiles(Chem.RemoveHs(poseH), isomericSmiles=True) or "NA"
            print(f"[info] OFFMol graph (from SMILES+MCS coords) SMILES_REF: {ref_smiles}", flush=True)
            return rdkit_for_pdb, offmol, smiles_log

        except Exception as e:
            print(f"[warn] OFF-from-SMILES+MCS-coords failed: {e}; will try template-based bond orders.", flush=True)

    # B) Template-based bond orders with tautomer enumeration
    if ref_smiles:
        ref0 = Chem.MolFromSmiles(ref_smiles)
        if ref0 is None:
            raise RuntimeError(f"Could not parse --ligand-smiles: {ref_smiles}")
        ref0 = _std_cleanup(ref0)
        ref_templates = [ref0]
        if HAS_STD:
            try:
                te = Std.TautomerEnumerator()
                tauts = list(te.Enumerate(ref0))[:64]
                if tauts:
                    ref_templates = tauts
            except Exception:
                pass
        pose_noH = Chem.RemoveHs(Chem.Mol(molH_pose))
        for i, ref in enumerate(ref_templates, start=1):
            try:
                ref_noH = Chem.RemoveHs(ref)
                bo_fixed = AllChem.AssignBondOrdersFromTemplate(ref_noH, pose_noH)
            except Exception as e:
                print(f"[info] Template {i}: AssignBondOrdersFromTemplate failed ({e}); trying next template...", flush=True)
                continue
            cand = Chem.AddHs(bo_fixed, addCoords=True)
            Chem.SanitizeMol(cand, catchErrors=True)
            cand.UpdatePropertyCache(False)
            if any(a.GetNumRadicalElectrons() for a in cand.GetAtoms()):
                print(f"[info] Template {i}: radicals remain after template assignment; trying next template.", flush=True)
                continue
            smiles = Chem.MolToSmiles(Chem.RemoveHs(cand), isomericSmiles=True) or "NA"
            print(f"[info] OFFMol graph (from template #{i}) SMILES: {smiles}", flush=True)
            offmol = Molecule.from_rdkit(cand, allow_undefined_stereo=True, hydrogens_are_explicit=True)
            return cand, offmol, smiles

        print("[warn] No tautomer template yielded a radical-free assignment; falling back to de-radicalization workflow.", flush=True)

    # C) De-radicalize + last-resort kekulization
    molH = _deradicalize_and_add_hs(Chem.Mol(molH_pose))
    base_for_smiles = Chem.RemoveHs(molH)
    smiles = Chem.MolToSmiles(base_for_smiles, isomericSmiles=True) or "NA"
    print(f"[info] OFFMol graph (logging only) SMILES: {smiles}", flush=True)
    try:
        offmol = Molecule.from_rdkit(molH, allow_undefined_stereo=True, hydrogens_are_explicit=True)
    except RadicalsNotSupportedError as e:
        try:
            kek = Chem.Mol(molH)
            Chem.Kekulize(kek, clearAromaticFlags=True)
            offmol = Molecule.from_rdkit(kek, allow_undefined_stereo=True, hydrogens_are_explicit=True)
        except Exception:
            tmp = Path(sdf_path).with_suffix(".deradicalized.sdf")
            w = Chem.SDWriter(str(tmp))
            try:
                w.write(molH)
            finally:
                w.close()
            raise RuntimeError(
                f"Ligand still contains radicals after cleanup; wrote triage SDF to {tmp}. Original error: {e}"
            ) from e
    return molH, offmol, smiles


# ---------- main pipeline ----------

def main():
    ap = argparse.ArgumentParser(description="Prepare solvated MD system (GPU-ready) from receptor + docked SDF pose")
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--sdf", required=True)
    ap.add_argument("--pose-index", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ph", type=float, default=7.4)
    ap.add_argument("--padding", type=float, default=10.0)
    ap.add_argument("--ionic", type=float, default=0.15)
    ap.add_argument("--neutralize", action="store_true")
    ap.add_argument("--platform", choices=["CUDA", "OpenCL", "CPU"])
    ap.add_argument("--stop-after", choices=STAGES, default="all")
    ap.add_argument("--ligand-smiles", help="Trusted ligand SMILES used to assign bond orders to the docked pose")
    ap.add_argument(
        "--auto-ligand-smiles",
        action="store_true",
        help="Auto-resolve a trusted SMILES from the selected SDF (SDF props → PubChem/ChEMBL/Cactus → RDKit canonical), "
             "validated with OpenFF; used only if --ligand-smiles is not provided.",
    )
    ap.add_argument("--ligand-name-hint", help="Optional name hint for PubChem/ChEMBL/CACTUS (defaults from SDF filename)")
    ap.add_argument("--water-model", default="tip3p", help="Water model for modeller.addSolvent (tip3p, tip3pfb, spce, tip4pew)")
    ap.add_argument("--strip-heme", action="store_true", help="Strip HEM/HEC from the receptor before parameterization")
    ap.add_argument("--strip-resnames", default="", help="Comma-separated list of residue names to strip (e.g., GOL,SO4,EDO)")

    args = ap.parse_args()
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)

    # Provenance: script hash
    script_path = Path(__file__).resolve()
    print(f"[info] md_gpu_prepare.py: {script_path} sha1={sha1sum(script_path)}", flush=True)

    rec_pdb = Path(args.receptor)
    sdf = Path(args.sdf)
    if not rec_pdb.exists():
        sys.exit(f"[FATAL] receptor not found: {rec_pdb}")
    if not sdf.exists():
        sys.exit(f"[FATAL] sdf not found: {sdf}")

    durations: dict[str, float] = {}
    t_total0 = time.perf_counter()

    # 1) ligand (optionally auto-resolve SMILES)
    ref_smiles = args.ligand_smiles
    if not ref_smiles and args.auto_ligand_smiles:
        try:
            ref_smiles = resolve_trusted_smiles(
                Path(args.sdf),
                pose_index=args.pose_index,
                ligand_name_hint=args.ligand_name_hint,
            )
            print("[info] Auto-resolved SMILES via resolver.", flush=True)
        except Exception as e:
            print(f"[warn] auto SMILES resolution failed: {e}", flush=True)

    with stage_timer("ligand", durations, outdir):
        rdkit_lig, off_lig, smiles = load_ligand_from_sdf(sdf, args.pose_index, ref_smiles=ref_smiles)
        lig_pdb = outdir / "ligand_pose.pdb"
        write_rdkit_pdb(rdkit_lig, lig_pdb)
        (outdir / "ligand_smiles.txt").write_text(smiles + "\n")
        if ref_smiles:
            (outdir / "ligand_smiles_ref.txt").write_text(ref_smiles + "\n")
        mark_ok(outdir, "ligand")
    if args.stop_after == "ligand":
        return 0

    # 2) receptor (fix hydrogens/missing atoms; optionally strip)
    with stage_timer("receptor", durations, outdir):
        fixed_path = outdir / "receptor_fixed.pdb"
        if HAS_FIXER:
            print("[info] running PDBFixer (addH, missing atoms)", flush=True)
            fixer = pdbfixer.PDBFixer(filename=str(rec_pdb))
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms()
            fixer.addMissingHydrogens(pH=args.ph)
            with fixed_path.open("w") as fh:
                app.PDBFile.writeFile(fixer.topology, fixer.positions, fh)
        else:
            print("[info] PDBFixer not available; copying receptor to receptor_fixed.pdb", flush=True)
            shutil.copy2(rec_pdb, fixed_path)

        prot = app.PDBFile(str(fixed_path))

        # Optional stripping (heme/resnames)
        to_strip = set()
        if args.strip_heme:
            to_strip.update({"HEM", "HEC"})
            print("[info] stripping HEM/HEC from receptor", flush=True)
        if args.strip_resnames:
            more = {r.strip().upper() for r in args.strip_resnames.split(",") if r.strip()}
            if more:
                print(f"[info] stripping residues: {','.join(sorted(more))}", flush=True)
                to_strip.update(more)
        if to_strip:
            top, pos = _strip_residues_inplace(prot.topology, prot.positions, to_strip)
            # overwrite the fixed receptor with stripped version
            with fixed_path.open("w") as fh:
                app.PDBFile.writeFile(top, pos, fh)
            prot = app.PDBFile(str(fixed_path))

        mark_ok(outdir, "receptor")
    if args.stop_after == "receptor":
        return 0

    # 3) merge
    with stage_timer("merge", durations, outdir):
        print("[info] merging receptor + ligand", flush=True)
        lig = app.PDBFile(str(outdir / "ligand_pose.pdb"))
        modeller = app.Modeller(prot.topology, prot.positions)
        modeller.add(lig.topology, lig.positions)
        with (outdir / "complex_merged.pdb").open("w") as fh:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)
        mark_ok(outdir, "merge")
    if args.stop_after == "merge":
        return 0

    # 4) forcefield
    with stage_timer("forcefield", durations, outdir):
        ff, heme_loaded = try_build_forcefield(off_lig, args.water_model)
        (outdir / "ff.info.json").write_text(json.dumps({"heme_ffxml": heme_loaded}, indent=2))
        mark_ok(outdir, "forcefield")
    if args.stop_after == "forcefield":
        return 0

    # 5) solvate
    with stage_timer("solvate", durations, outdir):
        print("[info] adding solvent and ions", flush=True)
        water_model_used = args.water_model
        try:
            modeller.addSolvent(
                ff,
                model=water_model_used,
                padding=args.padding * unit.angstroms,
                ionicStrength=args.ionic * unit.molar,
                neutralize=args.neutralize,
            )
        except Exception as e:
            print(f"[warn] addSolvent(model='{water_model_used}') failed ({e}); retrying with 'tip3p'", flush=True)
            water_model_used = "tip3p"
            modeller.addSolvent(
                ff,
                model=water_model_used,
                padding=args.padding * unit.angstroms,
                ionicStrength=args.ionic * unit.molar,
                neutralize=args.neutralize,
            )
        with (outdir / "complex_solvated_premin.pdb").open("w") as fh:
            app.PDBFile.writeFile(modeller.topology, modeller.positions, fh)
        mark_ok(outdir, "solvate")
    if args.stop_after == "solvate":
        return 0

    # 6) system
    with stage_timer("system", durations, outdir):
        print("[info] creating System", flush=True)
        system = ff.createSystem(
            modeller.topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0 * unit.nanometer,
            constraints=app.HBonds,
            hydrogenMass=3.0 * unit.amu,  # HMR
        )
        configure_simulation(system, args)
        (outdir / "system.xml").write_text(XmlSerializer.serialize(system))
        mark_ok(outdir, "system")
    if args.stop_after == "system":
        return 0

    # 7) minimize
    with stage_timer("minimize", durations, outdir):
        print("[info] minimizing", flush=True)
        integrator = mm.LangevinMiddleIntegrator(310 * unit.kelvin, 1.0 / unit.picosecond, 0.002 * unit.picoseconds)
        platform = pick_platform(args.platform)
        sim = app.Simulation(modeller.topology, system, integrator, platform)
        sim.context.setPositions(modeller.positions)
        sim.minimizeEnergy(maxIterations=2000)
        (outdir / "integrator_min.xml").write_text(XmlSerializer.serialize(integrator))
        state = sim.context.getState(getPositions=True, getVelocities=True)
        (outdir / "state_min.xml").write_text(XmlSerializer.serialize(state))
        with (outdir / "complex_solvated.pdb").open("w") as fh:
            app.PDBFile.writeFile(modeller.topology, sim.context.getState(getPositions=True).getPositions(), fh)
        mark_ok(outdir, "minimize")
    if args.stop_after == "minimize":
        return 0

    # 8) all (metadata + total time)
    total_s = time.perf_counter() - t_total0
    durations["total"] = total_s
    (outdir / "durations.json").write_text(json.dumps(durations, indent=2))

    meta = {
        "receptor": str(rec_pdb),
        "receptor_fixed": str(outdir / "receptor_fixed.pdb"),
        "sdf": str(sdf),
        "pose_index": args.pose_index,
        "ph": args.ph,
        "padding_A": args.padding,
        "ionic_M": args.ionic,
        "neutralize": args.neutralize,
        "water_model": water_model_used,
        "platform": platform.getName(),
        "heme_ffxml": heme_loaded,
        "ligand_smiles_ref": ref_smiles,
        "durations_s": durations,
    }
    (outdir / "prepare.json").write_text(json.dumps(meta, indent=2))
    mark_ok(outdir, "all")
    print(f"[time] total: {total_s:.3f}s", flush=True)
    print(f"[done] Prepared MD folder: {outdir}", flush=True)
    return 0


def configure_simulation(system, config):
    """Set up simulation based on membrane presence."""
    if is_membrane_system(config):
        # Membrane-specific settings
        system.setDefaultPeriodicBoxVectors(*membrane_box_vectors)
        system.addForce(Lipid14Force())
    else:
        # Standard aqueous settings
        system.setDefaultPeriodicBoxVectors(*water_box_vectors)
    
    # Common settings
    system.addForce(MonteCarloBarostat(1.0, config['temperature_K']))


if __name__ == "__main__":
    sys.exit(main())
