from __future__ import annotations
import pytest
pytestmark = pytest.mark.order(9)

import os
import pytest

from ovmpk.docking.gnina_wrapper import has_gnina

def test_gnina_presence_does_not_crash():
    # Should always return a boolean and never raise
    ok = has_gnina()
    assert isinstance(ok, bool)