"""A horizon declared by only part of the systems must not be read as it stands.

The pooled tables are a union over systems, so a label column that only some
configs declare arrives filled for those systems and missing for the rest.
Reading such a column turns the missing half into a value — labels that are
not labels — and the model trains on them without a word. The resolution step
must notice the gap and rebuild the labels instead.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "notebooks"))

import predict_lib as pl  # noqa: E402


def _pool(y_values):
    return pd.DataFrame({
        "system": ["a", "a", "b", "b"],
        "condition": ["c", "c", "c", "c"],
        "replica": ["1", "1", "1", "1"],
        "time_ps": [0.0, 100.0, 0.0, 100.0],
        "y_tau500": y_values,
        "valid_tau500": [True, True, True, True],
    })


def _stub(y, valid):
    def relabel(pool, tau_ps):
        return np.asarray(y), np.asarray(valid)
    return relabel


def test_partially_filled_column_is_rebuilt_not_read():
    pool = _pool([1.0, 0.0, np.nan, np.nan])  # system "b" never declared tau=500
    y, valid, source = pl.resolve_labels(
        pool, 500.0, relabel=_stub([1, 0, 0, 1], [True] * 4))
    assert "relabelled" in source and "part of the systems" in source
    assert list(y) == [1, 0, 0, 1]


def test_complete_column_is_used_as_is():
    pool = _pool([1.0, 0.0, 0.0, 1.0])
    y, valid, source = pl.resolve_labels(
        pool, 500.0, relabel=_stub([9, 9, 9, 9], [True] * 4))
    assert source == "pipeline tables"
    assert list(y) == [1, 0, 0, 1]


def test_absent_column_is_rebuilt():
    pool = _pool([1.0, 0.0, 0.0, 1.0]).drop(columns=["y_tau500", "valid_tau500"])
    _y, _v, source = pl.resolve_labels(
        pool, 500.0, relabel=_stub([0, 1, 0, 1], [True] * 4))
    assert source == "relabelled from events"


def test_non_binary_labels_are_refused():
    pool = _pool([1.0, 0.0, np.nan, np.nan])
    with pytest.raises(ValueError, match="not binary"):
        pl.resolve_labels(pool, 500.0,
                          relabel=_stub([1, 0, -9223372036854775808, 1], [True] * 4))


def test_self_check_reading_the_gap_marks_unknown_frames_negative():
    """The defect this guards against, demonstrated on the same frame.

    Reading the column directly turns the two unknown rows into negatives —
    'no crossing ahead' asserted about frames nothing is known about. The
    result still looks like a valid binary label, so a value check cannot
    catch it; only noticing the gap can.
    """
    pool = _pool([1.0, 0.0, np.nan, np.nan])
    naive = pool["y_tau500"].to_numpy(int)
    assert list(naive) == [1, 0, 0, 0]
    assert np.isin(np.unique(naive), (0, 1)).all()
