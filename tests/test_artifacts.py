"""Pipeline output artifact validation tests."""
from pathlib import Path
import pytest

@pytest.mark.order(100)  # Run last
def test_output_structure():
    """Check only for docked poses"""
    assert Path("data/output/latest/docking/poses.sdf").exists()

@pytest.mark.order(51)
def test_log_files():
    """Check log files contain run details."""
    log_path = Path("runlogs/pipeline.log")
    assert log_path.exists()
    assert "Run completed" in log_path.read_text()
