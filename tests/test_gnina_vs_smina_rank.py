from __future__ import annotations
import pytest
pytestmark = pytest.mark.order(11)

import os
import pytest
from pathlib import Path

from ovmpk.docking.gnina_wrapper import has_gnina

REPO = Path(__file__).resolve().parents[1]

@pytest.mark.skipif(
    not has_gnina() or os.environ.get("OVMPK_TEST_GNINA_REAL") != "1",
    reason="Set OVMPK_TEST_GNINA_REAL=1 and ensure gnina is installed to run this.",
)
def test_placeholder_compare_gnina_vs_smina():
    # Placeholder: the real test should dock with smina, rescore with gnina, and compare ranks.
    # We keep this test gated to avoid CI failures on machines without gnina.
    assert True