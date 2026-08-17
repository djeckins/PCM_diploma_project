"""Defect class: if the model layer can open a trajectory, one day it will.

Layer boundaries are drawn along imports and enforced by grepping the sources.
"""

import re
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src" / "pcm2"

MD_RE = re.compile(r"^\s*(import MDAnalysis|from MDAnalysis)", re.M)
ML_RE = re.compile(r"^\s*(import sklearn|from sklearn|import xgboost|from xgboost)", re.M)

# MD reading: only the reading layer and the feature modules that need PBC distances.
MD_ALLOWED = {"io.py", "autodetect.py", "events.py", "labels.py", "coords.py",
              "features/__init__.py"}
# The ML library must not appear in any feature module (nor in the reading layer).
ML_FORBIDDEN_DIRS = {"features"}
ML_FORBIDDEN_FILES = {"io.py", "pore.py", "events.py", "labels.py", "autodetect.py",
                      "inputs.py", "vdw.py", "config.py", "runtime.py"}


def _rel(p: Path) -> str:
    return str(p.relative_to(SRC))


def test_md_library_only_in_reading_layer():
    offenders = [_rel(p) for p in SRC.rglob("*.py")
                 if MD_RE.search(p.read_text()) and _rel(p) not in MD_ALLOWED]
    assert offenders == [], f"MDAnalysis outside the reading layer: {offenders}"


def test_ml_library_not_in_feature_or_reading_layer():
    offenders = []
    for p in SRC.rglob("*.py"):
        rel = _rel(p)
        in_forbidden = (rel in ML_FORBIDDEN_FILES
                        or rel.split("/")[0] in ML_FORBIDDEN_DIRS)
        if in_forbidden and ML_RE.search(p.read_text()):
            offenders.append(rel)
    assert offenders == [], f"ML library in the feature/reading layer: {offenders}"


def test_report_and_viz_do_not_train():
    for name in ("report.py", "viz.py"):
        text = (SRC / name).read_text()
        assert ".fit(" not in text, f"{name} trains a model — reporting must only read artifacts"
