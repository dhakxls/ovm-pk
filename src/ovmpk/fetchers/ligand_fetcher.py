import os
import urllib.parse
import requests
from pathlib import Path
from typing import Dict, Any, Optional

# PubChem PUG REST URLs for different identifier types
PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_URLS = {
    "name":   f"{PUBCHEM_BASE}/compound/name/{{identifier}}/SDF?record_type=3d",
    "cid":    f"{PUBCHEM_BASE}/compound/cid/{{identifier}}/SDF?record_type=3d",
    "smiles": f"{PUBCHEM_BASE}/compound/smiles/{{identifier}}/SDF?record_type=3d",
    "inchikey": f"{PUBCHEM_BASE}/compound/inchikey/{{identifier}}/SDF?record_type=3d",
    # Add other types like InChI if needed
}

# Placeholder content for dry runs
DRY_RUN_SDF_CONTENT = "$$$$ placeholder SDF for dry run" # Placeholder content for SDF

def _dry_sdf_file(path: Path, identifier: str, contents: str = DRY_RUN_SDF_CONTENT):
    """Creates a placeholder SDF file for dry runs."""
    path.parent.mkdir(parents=True, exist_ok=True)
    # Basic SDF structure requires at least molecule name and $$$$ delimiter
    placeholder_text = f"{identifier}\n  Placeholder\n\n  0  0  0  0  0  0  0  0  0  0 V2000\nM  END\n{contents}\n"
    path.write_text(placeholder_text)

def _is_ligand_placeholder(path: Path) -> bool:
    """Checks if a file contains only the dry run placeholder SDF content."""
    if not path.exists() or path.stat().st_size > 200: # Placeholder SDF is small
        return False
    try:
        # Check if the specific placeholder comment line exists
        content = path.read_text()
        return DRY_RUN_SDF_CONTENT in content
    except Exception:
        return False # Error reading, treat as not a placeholder

def fetch(identifier_value: str, cfg: Dict[str, Any]) -> Path:
    """
    Fetches a ligand 3D structure SDF from PubChem using various identifiers.

    Args:
        identifier_value: The value of the identifier (e.g., 'ketoconazole', '12345', 'C1=CC=...').
                          Passed directly from the --ligand CLI argument.
        cfg: Configuration dictionary, expects settings under 'fetch.ligand'.

    Returns:
        Path to the downloaded (or placeholder) SDF file.
    """
    dry = bool(os.getenv("OVM_DRY_RUN", ""))
    fetch_cfg = cfg.get("fetch", {}).get("ligand", {})

    # Determine identifier type and value
    # Default to interpreting the CLI --ligand argument as a 'name' unless specified otherwise
    identifier_type = fetch_cfg.get("identifier_type", "name").lower()
    # Use the value passed from the CLI directly
    identifier_for_fetch = identifier_value

    if identifier_type not in PUBCHEM_URLS:
        raise ValueError(f"Unsupported ligand identifier_type in config: '{identifier_type}'. Supported types: {list(PUBCHEM_URLS.keys())}")

    # Use a clean filename based on the identifier, replacing tricky characters for filesystem
    safe_filename = "".join(c if c.isalnum() or c in ('-', '_') else '_' for c in identifier_value)
    # Cap length to avoid excessively long filenames (e.g., from SMILES)
    max_len = 50
    if len(safe_filename) > max_len:
         safe_filename = safe_filename[:max_len] + "_truncated"

    filename = f"{safe_filename}_{identifier_type}.sdf"
    base = Path("data/input/ligands")
    base.mkdir(parents=True, exist_ok=True)
    outp = base / filename

    if dry:
        print(f"[info][dry-run] Creating placeholder SDF ({identifier_type}={identifier_value}) at {outp}")
        _dry_sdf_file(outp, identifier_value)
        return outp # Return placeholder path
    else:
        should_download = True
        if outp.exists():
            if _is_ligand_placeholder(outp):
                print(f"[info] Placeholder SDF file found at {outp}. Will overwrite.")
                try:
                    outp.unlink() # Delete placeholder
                except OSError as e:
                    print(f"[warn] Could not delete placeholder SDF file {outp}: {e}. Download may fail.")
            else:
                print(f"[info] Actual SDF file found at {outp}. Skipping download.")
                should_download = False # File exists and is not the placeholder

        if should_download:
            print(f"[info] Downloading {identifier_type} '{identifier_for_fetch}' SDF from PubChem to {outp}...")
            try:
                # URL encode the identifier for safety, especially for names/SMILES
                encoded_identifier = urllib.parse.quote(identifier_for_fetch)
                url_template = PUBCHEM_URLS[identifier_type]
                url = url_template.format(identifier=encoded_identifier)

                r = requests.get(url, timeout=60)
                r.raise_for_status() # Check for HTTP errors (like 404 Not Found)

                # Basic check for valid SDF content (presence of $$$$)
                content_bytes = r.content
                # Allow empty file if PubChem returns 200 OK but no structure (can happen for abstract compounds)
                # Downstream RDKit/obabel steps should handle empty SDF if necessary
                if content_bytes and b"$$$$" not in content_bytes[-100:] and len(content_bytes) > 10: # Check if content seems SDF-like
                     # Try decoding to check for HTML error messages
                     try:
                         text_content = content_bytes.decode('utf-8', errors='ignore')
                         if '<!DOCTYPE html>' in text_content or '<html>' in text_content:
                             raise ValueError(f"Downloaded content for {identifier_type} '{identifier_for_fetch}' appears to be an HTML error page, not SDF.")
                     except UnicodeDecodeError:
                         pass # Ignore if it wasn't utf-8 text

                     # If not obviously HTML, still raise a general warning/error
                     raise ValueError(f"Downloaded content for {identifier_type} '{identifier_for_fetch}' does not look like valid SDF (missing $$$$ marker?). Size: {len(content_bytes)}")


                outp.write_bytes(content_bytes) # Use write_bytes as SDF is often text but can have encoding issues
                if content_bytes:
                     print(f"[info] Successfully downloaded {identifier_type} '{identifier_for_fetch}' SDF.")
                else:
                     print(f"[warn] Download successful (200 OK) but received empty content for {identifier_type} '{identifier_for_fetch}'. Downstream tools might fail.")

            except requests.exceptions.HTTPError as e:
                 # Specifically catch HTTP errors like 404
                 raise RuntimeError(f"Failed to download SDF for {identifier_type} '{identifier_for_fetch}' from {url}. Status Code: {e.response.status_code}. Check identifier.") from e
            except requests.exceptions.RequestException as e:
                raise RuntimeError(f"Network or request error downloading SDF for {identifier_type} '{identifier_for_fetch}' from {url}: {e}") from e
            except ValueError as e:
                # Reraise content validation errors
                raise RuntimeError(str(e)) from e
            except Exception as e:
                # Catch any other unexpected errors
                raise RuntimeError(f"An unexpected error occurred during SDF download for {identifier_type} '{identifier_for_fetch}': {e}") from e

    return outp # Return the path to the SDF file

