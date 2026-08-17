"""Defect classes: preprocessing fitted before the split; window looks forward;
frame-level split; one-sided embargo."""

import numpy as np
import pytest

from pcm2.features._common import backward_window_slice
from pcm2.models.nanpolicy import MissingEncoder
from pcm2.splits import check_no_frame_leak, outer_folds, time_block_folds


def test_imputer_learns_only_on_fold_training_rows():
    X_tr = np.array([[1.0], [3.0]])
    X_all = np.array([[1.0], [3.0], [100.0]])
    enc = MissingEncoder().fit(X_tr)
    filled = enc.transform(np.array([[np.nan]]))[0, 0]
    assert filled == 2.0  # train median, not the whole-matrix median (33.5 would leak)
    enc_bad = MissingEncoder().fit(X_all)
    assert enc_bad.transform(np.array([[np.nan]]))[0, 0] != 2.0  # the defect is caught


def test_backward_window_never_includes_future():
    times = np.arange(0.0, 100.0, 10.0)
    for i in range(len(times)):
        sl = backward_window_slice(times, i, 25.0)
        assert sl.stop == i + 1  # no future frames in the window
        assert times[i] - times[sl.start] <= 25.0


def test_outer_folds_never_split_a_group():
    groups = np.array(["a", "a", "b", "b", "c"])
    for tr, te, g in outer_folds(groups):
        check_no_frame_leak(tr, te, groups)
        assert set(groups[te]) == {g}


def test_frame_level_split_raises():
    groups = np.array(["a", "a", "b"])
    with pytest.raises(RuntimeError, match="leakage"):
        check_no_frame_leak(np.array([0, 2]), np.array([1]), groups)


def test_time_embargo_is_two_sided():
    times = np.arange(0.0, 1000.0, 10.0)
    folds = time_block_folds(times, n_blocks=4, embargo_ps=50.0)
    for tr, te in folds:
        lo, hi = times[te].min(), times[te].max()
        # Training right before the boundary knows about an event in the test
        # block (the label looks forward); right after the boundary the test
        # frame is built from training frames (the window looks backward).
        assert not np.any((times[tr] >= lo - 50.0) & (times[tr] < lo))
        assert not np.any((times[tr] > hi) & (times[tr] <= hi + 50.0))
