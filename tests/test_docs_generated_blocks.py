"""Defect class: a number typed into a document by hand; the acceptance table
diverging from its source. The test must also check itself: corrupt a block
and verify the divergence is caught."""

import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _run_check(cwd) -> int:
    return subprocess.run([sys.executable, str(ROOT / "tools" / "genblocks.py"),
                           "--check"], cwd=cwd, capture_output=True).returncode


def test_generated_blocks_match_artifacts():
    assert _run_check(ROOT) == 0, "docs/ diverge from artifacts: run python tools/genblocks.py"


def test_corrupted_block_is_detected(tmp_path):
    # Self-check: corruption inside a generated block must be caught.
    sandbox = tmp_path / "proj"
    (sandbox / "tools").mkdir(parents=True)
    shutil.copy(ROOT / "tools" / "genblocks.py", sandbox / "tools" / "genblocks.py")
    shutil.copytree(ROOT / "docs", sandbox / "docs")
    if (ROOT / "runs").exists():
        shutil.copytree(ROOT / "runs", sandbox / "runs",
                        ignore=shutil.ignore_patterns("*.parquet", "*.pkl", "*.png"))
    doc = sandbox / "docs" / "RESULTS.md"
    text = doc.read_text()
    marker = "<!-- BEGIN GENERATED: results -->"
    i = text.index(marker) + len(marker)
    doc.write_text(text[:i] + "\nCORRUPTED 0.999\n" + text[i:])
    assert subprocess.run([sys.executable, str(sandbox / "tools" / "genblocks.py"),
                           "--check"], capture_output=True).returncode == 1
