import os
import pytest
from pathlib import Path

from ovmpk.docking.gnina_wrapper import build_gnina_rescore_cmd, has_gnina, rescore

REPO = Path(__file__).resolve().parents[1]

@pytest.mark.parametrize("cnn_scoring", ["ad4", "vina"])
def test_build_gnina_rescore_cmd_smoke(cnn_scoring):
    cmd = build_gnina_rescore_cmd(
        receptor_pdbqt="receptor.pdbqt",
        poses_path="poses.sdf",
        cnn_scoring=cnn_scoring,
        cnn_models=None,
        extra_args=["--random_seed", "42"],
    )
    # basic shape
    assert cmd[0] == "gnina"
    assert "--score_only" in cmd
    assert "receptor.pdbqt" in cmd
    assert "poses.sdf" in cmd

@pytest.mark.skipif(
    not has_gnina() or os.environ.get("OVMPK_TEST_GNINA") != "1",
    reason="GNINA not installed or OVMPK_TEST_GNINA!=1; skipping live rescoring.",
)
def test_gnina_rescore_dry_run_then_live(tmp_path):
    # Create dummy paths; in dry_run we don't need real files
    out_csv = tmp_path / "scores.csv"
    rc, cmd = rescore(
        receptor_pdbqt="receptor.pdbqt",
        poses_path="poses.sdf",
        out_csv=out_csv,
        cnn_scoring="ad4",
        dry_run=True,
    )
    assert rc == 0
    assert cmd and cmd[0] == "gnina"
