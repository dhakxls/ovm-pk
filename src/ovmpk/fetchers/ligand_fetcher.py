# src/ovmpk/fetchers/ligand_fetcher.py
import urllib.parse
import requests
from pathlib import Path
from typing import Dict, Any

# PubChem PUG REST URLs for different identifier types
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_URLS = {
    "name":     f"{PUBCHEM_BASE}/compound/name/{{identifier}}/SDF?record_type=3d",
    "cid":      f"{PUBCHEM_BASE}/compound/cid/{{identifier}}/SDF?record_type=3d",
    "smiles":   f"{PUBCHEM_BASE}/compound/smiles/{{identifier}}/SDF?record_type=3d",
    "inchikey": f"{PUBCHEM_BASE}/compound/inchikey/{{identifier}}/SDF?record_type=3d",
    # add more types if needed (e.g., "inchi")
}

def _bytes_look_like_sdf(b: bytes) -> bool:
    """Basic SDF sanity check: presence of '$$$$' (record separator)."""
    if not b or len(b) < 5:
        return False
    tail = b[-4096:] if len(b) > 4096 else b
    return b"$$$$" in tail

def _file_looks_like_sdf(path: Path) -> bool:
    try:
        data = path.read_bytes()
    except OSError:
        return False
    return _bytes_look_like_sdf(data)

def fetch(identifier_value: str, cfg: Dict[str, Any]) -> Path:
    """
    Fetch a ligand 3D SDF from PubChem using various identifiers.

    Args:
        identifier_value: The identifier value (e.g., 'ketoconazole', '12345', 'C1=CC=...')
        cfg: Configuration dict. Expects 'fetch.ligand.identifier_type' (default: 'name').

    Returns:
        Path to the downloaded SDF file.
    """
    fetch_cfg = cfg.get("fetch", {}).get("ligand", {})

    # Determine identifier type (defaults to treating CLI --ligand as a name)
    identifier_type = str(fetch_cfg.get("identifier_type", "name")).lower()
    if identifier_type not in PUBCHEM_URLS:
        raise ValueError(
            f"Unsupported ligand identifier_type: '{identifier_type}'. "
            f"Supported: {list(PUBCHEM_URLS.keys())}"
        )

    # Sanitize filename from the raw identifier (avoid giant names from SMILES)
    safe = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in identifier_value)
    if len(safe) > 50:
        safe = safe[:50] + "_truncated"
    filename = f"{safe}_{identifier_type}.sdf"

    base = Path("data/input/ligands")
    base.mkdir(parents=True, exist_ok=True)
    outp = base / filename

    # If an SDF already exists and looks valid, reuse it
    if outp.exists() and _file_looks_like_sdf(outp):
        print(f"[info] Actual SDF file found at {outp}. Skipping download.")
        return outp

    # Otherwise (missing or invalid), download a fresh copy
    encoded_identifier = urllib.parse.quote(identifier_value)
    url = PUBCHEM_URLS[identifier_type].format(identifier=encoded_identifier)

    print(f"[info] Downloading {identifier_type} '{identifier_value}' SDF from PubChem to {outp}...")
    try:
        r = requests.get(url, timeout=60)
        r.raise_for_status()
        content = r.content

        # Validate: PubChem SDFs should contain '$$$$'
        if not _bytes_look_like_sdf(content):
            # It might be an HTML error page or an empty/invalid SDF
            try:
                text_preview = content[:2000].decode("utf-8", errors="ignore")
            except Exception:
                text_preview = "<binary>"
            raise RuntimeError(
                "Downloaded content does not look like a valid SDF "
                f"(no '$$$$' terminator). URL: {url}\nPreview:\n{text_preview}"
            )

        outp.write_bytes(content)
        print(f"[info] Successfully downloaded SDF to {outp}")
        return outp

    except requests.exceptions.HTTPError as e:
        raise RuntimeError(
            f"HTTP error fetching SDF for {identifier_type}='{identifier_value}' from {url} "
            f"(status={e.response.status_code})."
        ) from e
    except requests.exceptions.RequestException as e:
        raise RuntimeError(
            f"Network error fetching SDF for {identifier_type}='{identifier_value}' from {url}: {e}"
        ) from e
