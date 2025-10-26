import os
os.environ["OVM_DRY_RUN"] = "1"

from ovmpk.cli import pipeline_mvp

def test_pipeline_smoke(tmp_path):
    out = pipeline_mvp("configs/fast_test.yaml", "CYP3A4", "ketoconazole")
    assert "results" in out
