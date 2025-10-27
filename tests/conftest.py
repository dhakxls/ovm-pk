# tests/conftest.py
from __future__ import annotations
import os
from pathlib import Path
from datetime import datetime
import pytest
import sys, time
from math import ceil
# --- repo & run folders -------------------------------------------------------
REPO = Path(__file__).resolve().parents[1]
RUN_ID = datetime.now().strftime("%Y%m%d-%H%M%S")
OUT_ROOT = REPO / "data" / "output"
RUN_ROOT = OUT_ROOT / f"run-{RUN_ID}"
LOG_ROOT = RUN_ROOT / "logs"
for p in (OUT_ROOT, RUN_ROOT, LOG_ROOT):
    p.mkdir(parents=True, exist_ok=True)

# --- optional environment banner in pytest header ------------------------------
def pytest_report_header(config):
    lines = [f"OVM-PK test run: {RUN_ID}", f"Output root: {RUN_ROOT}"]
    try:
        import openmm as mm
        plats = [mm.Platform.getPlatform(i).getName()
                 for i in range(mm.Platform.getNumPlatforms())]
        lines.append("OpenMM platforms: " + ", ".join(plats))
    except Exception as e:
        lines.append(f"OpenMM platforms: n/a ({e})")
    try:
        import rdkit
        lines.append(f"RDKit: {rdkit.__version__}")
    except Exception:
        pass
    try:
        import rich
        lines.append(f"Rich: {rich.__version__}")
    except Exception:
        pass
    return "\n".join(lines)

# --- session bootstrap: chdir + unbuffered stdout for live progress -----------
@pytest.hookimpl(tryfirst=True)
def pytest_sessionstart(session):
    os.chdir(str(REPO))
    os.environ.setdefault("PYTHONUNBUFFERED", "1")
    os.environ.setdefault("OVMPK_RUN_ID", RUN_ID)
    os.environ.setdefault("OVMPK_OUT_ROOT", str(RUN_ROOT))
    os.environ.setdefault("OVMPK_LOG_ROOT", str(LOG_ROOT))
    try:
        from rich.console import Console
        Console().rule(f"[bold cyan]OVM-PK tests — run {RUN_ID}")
    except Exception:
        print("=" * 60)
        print(f"OVM-PK tests — run {RUN_ID}")
        print("=" * 60)

# --- handy fixture for tests that want the paths ------------------------------
@pytest.fixture(scope="session")
def run_paths():
    return {
        "repo": REPO,
        "run_id": RUN_ID,
        "out_root": RUN_ROOT,
        "log_root": LOG_ROOT,
    }

_started = None
_total = None
_done = 0

def pytest_sessionstart(session):
    global _started
    _started = time.time()

def pytest_collection_finish(session):
    global _total
    _total = len(session.items)

def _print_progress():
    if not _total:
        return
    elapsed = time.time() - _started
    avg = elapsed / max(_done, 1)
    eta = avg * max(_total - _done, 0)
    bar_len = 24
    ratio = _done / _total
    filled = ceil(bar_len * ratio)
    bar = "█" * filled + "·" * (bar_len - filled)
    msg = f"\r[ {bar} ] {_done}/{_total} | {elapsed:6.1f}s elapsed | ETA {eta:6.1f}s"
    sys.stderr.write(msg)
    sys.stderr.flush()

def pytest_runtest_logreport(report):
    global _done
    if report.when == "call" and report.passed or report.failed or report.skipped:
        _done += 1
        _print_progress()

def pytest_sessionfinish(session, exitstatus):
    # move to a new line after final progress line
    sys.stderr.write("\n")
    sys.stderr.flush()
