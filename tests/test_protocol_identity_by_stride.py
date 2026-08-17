"""Defect class: a run with one frame stride writing into another stride's
directory — these are different results and must not be compared."""

import yaml

from pcm2 import config
from pcm2.runtime import run_dir


def _cfg(tmp_path, stride, name):
    doc = {"system": {"id": "toy"},
           "data": {"data_root": "/x", "stride": stride, "conditions": [
               {"id": "c", "replicas": [{"id": "r", "topology": "/t.tpr",
                                         "trajectory": "/t.xtc", "lineage": "c/r"}]}]},
           "labels": {"tau_ps": [1000.0], "primary_tau_ps": 1000.0},
           "model": {"seed": 1}}
    p = tmp_path / name
    p.write_text(yaml.safe_dump(doc, allow_unicode=True))
    return config.load(p)


def test_run_dir_carries_stride_in_name(tmp_path):
    a = run_dir(_cfg(tmp_path, 3, "a.yaml"))
    b = run_dir(_cfg(tmp_path, 5, "b.yaml"))
    assert a != b
    assert a.name.endswith("stride3") and b.name.endswith("stride5")
