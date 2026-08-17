"""Shared fixtures. Tests must not write into the published results directory:
a snapshot of runs/ is checked at session end; corruption fails the whole suite."""

from __future__ import annotations

from pathlib import Path

import pytest

RUNS = Path(__file__).resolve().parents[1] / "runs"


def _tree_state(root: Path) -> dict[str, tuple[float, int]]:
    if not root.exists():
        return {}
    out = {}
    for p in sorted(root.rglob("*")):
        if p.is_file():
            out[str(p)] = (p.stat().st_mtime, p.stat().st_size)
    return out


@pytest.fixture(scope="session", autouse=True)
def runs_dir_untouched():
    before = _tree_state(RUNS)
    yield
    after = _tree_state(RUNS)
    changed = {k for k in set(before) | set(after)
               if before.get(k) != after.get(k)}
    assert not changed, (f"tests modified the published results directory: "
                         f"{sorted(changed)[:5]}")
