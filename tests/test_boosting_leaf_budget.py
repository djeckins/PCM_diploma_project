"""Defect class: a leaf threshold set as a raw number lets a leaf isolate a
fraction of one event; and 1/4 substituted for the measured initial guess."""

import numpy as np

from pcm2.models.trees import count_events, leaf_threshold


def test_events_are_runs_within_trajectory():
    y = np.array([0, 1, 1, 0, 1, 1, 1, 0])
    traj = np.array(["a"] * 8)
    assert count_events(y, traj) == 2
    # A run cut by a trajectory boundary is two events, not one.
    traj2 = np.array(["a", "a", "a", "a", "b", "b", "b", "b"])
    y2 = np.array([0, 0, 1, 1, 1, 1, 0, 0])
    assert count_events(y2, traj2) == 2


def test_threshold_uses_measured_base_score_not_quarter():
    y = np.concatenate([np.ones(10), np.zeros(990)])
    traj = np.array(["a"] * 1000)
    thr, math = leaf_threshold(k_events=2.0, y_tr=y, traj_tr=traj)
    p0 = 0.01
    assert abs(math["base_score_p0"] - p0) < 1e-9
    assert abs(math["h0"] - p0 * (1 - p0)) < 1e-12
    # The "h0 = 1/4" defect would change the threshold severalfold.
    assert abs(thr - 2.0 * 10 * p0 * (1 - p0)) < 1e-9
    assert thr < 2.0 * 10 * 0.25 / 2


def test_threshold_scales_with_frames_per_event():
    traj = np.array(["a"] * 100)
    y_short = np.array([1, 0] * 50)            # 50 one-frame events
    y_long = np.array(([1] * 10 + [0] * 10) * 5)  # 5 ten-frame events
    thr_s, _ = leaf_threshold(2.0, y_short, traj)
    thr_l, _ = leaf_threshold(2.0, y_long, traj)
    assert thr_l > thr_s  # same k_events → larger threshold for longer events
