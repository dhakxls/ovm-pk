# scripts/utils/smiles_fetch.py
from __future__ import annotations
import os, re, json, time
from pathlib import Path
from typing import Optional, Tuple, List

# stdlib HTTP (no extra deps)
import urllib.request, urllib.error

from rdkit import Chem
from rdkit.Chem import rdmolops
try:
    from rdkit.Chem import inchi as RDInchi  # optional RDKit InChI support
    HAS_INCHI = True
except Exception:
    HAS_INCHI = False

# Optional standardization
try:
    from rdkit.Chem.MolStandardize import rdMolStandardize as Std
    HAS_STD = True
except Exception:
    try:
        from rdkit.Chem import rdMolStandardize as Std
        HAS_STD = True
    except Exception:
        HAS_STD = False

CACHE_DIR = Path("data/cache/smiles")
CACHE_DIR.mkdir(parents=True, exist_ok=True)

def _http_get(url: str, timeout: float = 6.0) -> Optional[bytes]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as r:
            return r.read()
    except urllib.error.HTTPError:
        return None
    except urllib.error.URLError:
        return None

def _candidate_names_from_filename(path: Path) -> List[str]:
    # e.g., "..._ketoconazole_name_prepared_..." -> ["ketoconazole"]
    base = path.stem.lower()
    picks = set()
    for token in re.split(r"[^a-z0-9]+", base):
        if len(token) >= 5 and token.isalpha():
            picks.add(token)
    # Favor tokens near “name” or “ligand”
    ordered = sorted(picks, key=lambda t: (t not in base.split("_"), len(t)))
    return ordered

def _sdf_props_smiles_or_name(sdf_path: Path, pose_index: int = 1) -> Tuple[Optional[str], Optional[str]]:
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    i = max(0, pose_index - 1)
    mol = suppl[i] if (0 <= i < len(suppl)) else None
    if mol is None:
        return None, None
    # try direct SMILES tag first
    for k in ("SMILES", "IsomericSMILES", "CANONICAL_SMILES", "PUBCHEM_OPENEYE_CAN_SMILES"):
        if mol.HasProp(k):
            s = mol.GetProp(k).strip()
            if s:
                return s, None
    # else try a name-like tag
    for k in ("NAME", "MOLNAME", "PUBCHEM_IUPAC_TRADITIONAL_NAME", "DRUG_NAME", "TITLE"):
        if mol.HasProp(k):
            nm = mol.GetProp(k).strip()
            if nm:
                return None, nm
    return None, None

def _std_cleanup(mol: Chem.Mol) -> Chem.Mol:
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

def _is_openff_friendly(smiles: str) -> bool:
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return False
    m = _std_cleanup(m)
    # radicals?
    for a in m.GetAtoms():
        if a.GetNumRadicalElectrons():
            return False
    # minimal valence sanity check
    try:
        Chem.SanitizeMol(m, catchErrors=True)
        m.UpdatePropertyCache(strict=False)
    except Exception:
        pass
    return True

def _cache_key(term: str) -> Path:
    return CACHE_DIR / (re.sub(r"[^A-Za-z0-9._-]+", "_", term) + ".txt")

def _cache_get(term: str) -> Optional[str]:
    p = _cache_key(term)
    if p.exists():
        return p.read_text().strip() or None
    return None

def _cache_put(term: str, smiles: str) -> None:
    _cache_key(term).write_text(smiles.strip() + "\n")

# -------- PubChem --------

def _pubchem_smiles_by_inchikey(inchikey: str) -> Optional[str]:
    term = f"pubchem_inchikey_{inchikey}"
    c = _cache_get(term)
    if c:
        return c
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/property/IsomericSMILES,CanonicalSMILES/JSON"
    data = _http_get(url)
    if not data:
        return None
    try:
        js = json.loads(data)
        props = js["PropertyTable"]["Properties"][0]
        for k in ("IsomericSMILES", "CanonicalSMILES"):
            s = props.get(k)
            if s and _is_openff_friendly(s):
                _cache_put(term, s)
                return s
    except Exception:
        return None
    return None

def _pubchem_smiles_by_name(name: str) -> Optional[str]:
    term = f"pubchem_name_{name.lower()}"
    c = _cache_get(term)
    if c:
        return c
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{urllib.parse.quote(name)}/property/IsomericSMILES,CanonicalSMILES/JSON"
    data = _http_get(url)
    if not data:
        return None
    try:
        js = json.loads(data)
        props = js["PropertyTable"]["Properties"][0]
        for k in ("IsomericSMILES", "CanonicalSMILES"):
            s = props.get(k)
            if s and _is_openff_friendly(s):
                _cache_put(term, s)
                return s
    except Exception:
        return None
    return None

# -------- ChEMBL (fallback) --------

def _chembl_smiles_by_name(name: str) -> Optional[str]:
    term = f"chembl_name_{name.lower()}"
    c = _cache_get(term)
    if c:
        return c
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?molecule_synonyms__icontains={urllib.parse.quote(name)}&limit=5"
    data = _http_get(url)
    if not data:
        return None
    try:
        js = json.loads(data)
        for rec in js.get("molecules", []):
            structs = rec.get("molecule_structures") or {}
            s = structs.get("canonical_smiles")
            if s and _is_openff_friendly(s):
                _cache_put(term, s)
                return s
    except Exception:
        return None
    return None

def _inchikey_from_mol(mol: Chem.Mol) -> Optional[str]:
    if not HAS_INCHI:
        return None
    try:
        return RDInchi.MolToInchiKey(Chem.RemoveHs(mol))
    except Exception:
        return None

def get_trusted_smiles(
    sdf_path: Path,
    pose_index: int = 1,
    guess_from_filename: bool = True,
) -> Tuple[Optional[str], str]:
    """
    Returns (smiles, source) or (None, reason).
    Order:
      1) SDF property SMILES
      2) PubChem by InChIKey (if RDKit can compute it)
      3) PubChem by name (from SDF or filename tokens)
      4) ChEMBL by name (fallback)
    All candidates must pass _is_openff_friendly().
    """
    # 1) SDF property SMILES
    sdf_smiles, sdf_name = _sdf_props_smiles_or_name(sdf_path, pose_index)
    if sdf_smiles and _is_openff_friendly(sdf_smiles):
        return sdf_smiles, "sdf:SMILES"

    # Prepare RDKit mol for InChIKey
    mol = None
    try:
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
        idx0 = max(0, pose_index - 1)
        mol = suppl[idx0] if (0 <= idx0 < len(suppl)) else None
    except Exception:
        mol = None

    # 2) PubChem by InChIKey
    if mol is not None:
        ik = _inchikey_from_mol(mol)
        if ik:
            s = _pubchem_smiles_by_inchikey(ik)
            if s:
                return s, f"pubchem:inchikey:{ik}"

    # 3) PubChem by name (SDF prop or filename token)
    name_candidates = []
    if sdf_name:
        name_candidates.append(sdf_name)
    if guess_from_filename:
        name_candidates.extend(_candidate_names_from_filename(sdf_path))
    tried = set()
    for nm in name_candidates:
        nm = nm.strip()
        if not nm or nm in tried:
            continue
        tried.add(nm)
        s = _pubchem_smiles_by_name(nm)
        if s:
            return s, f"pubchem:name:{nm}"

    # 4) ChEMBL by name (fallback)
    for nm in name_candidates:
        s = _chembl_smiles_by_name(nm)
        if s:
            return s, f"chembl:name:{nm}"

    return None, "no-hit"
