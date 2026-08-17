"""Defect classes: a result assembled from pieces of different protocols; an
artifact outliving its purpose; a step leaving partial state behind."""

import pytest
import yaml

from pcm2 import config
from pcm2.runtime import run_dir, step_output


@pytest.fixture()
def cfg(tmp_path, monkeypatch):
    doc = {"system": {"id": "toy"},
           "data": {"data_root": "/x", "stride": 3, "conditions": [
               {"id": "c", "replicas": [{"id": "r", "topology": "/t.tpr",
                                         "trajectory": "/t.xtc", "lineage": "c/r"}]}]},
           "labels": {"tau_ps": [1000.0], "primary_tau_ps": 1000.0},
           "model": {"seed": 1}}
    p = tmp_path / "cfg.yaml"
    p.write_text(yaml.safe_dump(doc, allow_unicode=True))
    c = config.load(p)
    import pcm2.runtime as rt
    monkeypatch.setattr(rt, "PROJECT_ROOT", tmp_path)
    return c


def test_step_replaces_directory_and_removes_foreign_files(cfg):
    with step_output(cfg, "events") as out:
        (out / "events.parquet").write_text("v1")
    step_dir = run_dir(cfg) / "events"
    (step_dir / "stale_leftover.json").write_text("artifact that outlived its purpose")
    with step_output(cfg, "events") as out:
        (out / "events.parquet").write_text("v2")
    files = sorted(p.name for p in step_dir.iterdir())
    assert "stale_leftover.json" not in files, "a step must delete what it does not produce"
    assert (step_dir / "events.parquet").read_text() == "v2"
    assert (step_dir / "PROVENANCE.json").exists(), "provenance sits next to the output"


def test_failed_step_leaves_previous_output_intact(cfg):
    with step_output(cfg, "events") as out:
        (out / "events.parquet").write_text("good")
    with pytest.raises(RuntimeError):
        with step_output(cfg, "events") as out:
            (out / "events.parquet").write_text("partial")
            raise RuntimeError("failure mid-step")
    step_dir = run_dir(cfg) / "events"
    assert (step_dir / "events.parquet").read_text() == "good", \
        "a failed step must not touch the previous complete file set"
    leftovers = [p for p in step_dir.parent.iterdir() if p.name.startswith(".events.tmp")]
    assert leftovers == []


def test_provenance_records_config_and_code(cfg):
    with step_output(cfg, "labels") as out:
        (out / "labels.parquet").write_text("x")
    import json
    doc = json.loads((run_dir(cfg) / "labels" / "PROVENANCE.json").read_text())
    assert doc["config"]["system"]["id"] == "toy"
    assert "tree_sha256" in doc["code"] and doc["libraries"]["numpy"] != "absent"
