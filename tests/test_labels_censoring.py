"""Defect classes: an incomplete window silently labelled zero (deflates the
base rate); an event without an observed entry given its exit time in
entrance labels."""

import numpy as np

from pcm2.labels import label_one_trajectory


def test_label_definition_half_open_interval():
    times = np.array([0.0, 10.0, 20.0, 30.0])
    lab = label_one_trajectory(times, np.array([25.0]), np.array([]), 20.0, 1000.0)
    # y=1 ⟺ t < 25 ≤ t+20 → frames 10 and 20; frame 0: 25>20 → 0; frame 30: 25<30 → 0.
    assert lab["y_tau20"].tolist() == [0, 1, 1, 0]


def test_incomplete_negative_window_is_unknown_not_zero():
    times = np.array([0.0, 50.0, 90.0])
    lab = label_one_trajectory(times, np.array([]), np.array([]), 20.0, 100.0)
    assert lab["valid_tau20"].tolist() == [True, True, False]


def test_confirmed_positive_never_censored():
    times = np.array([90.0])
    lab = label_one_trajectory(times, np.array([95.0]), np.array([]), 20.0, 100.0)
    assert lab["y_tau20"].tolist() == [1]
    assert lab["valid_tau20"].tolist() == [True]  # window incomplete, but y=1 is confirmed


def test_entrance_anchor_uses_only_observed_entries():
    times = np.array([0.0, 10.0])
    # Entry unobserved → entry array empty: event excluded, not substituted by its exit.
    lab = label_one_trajectory(times, np.array([15.0]), np.array([]), 20.0, 1000.0)
    assert lab["y_entr_tau20"].sum() == 0
    lab2 = label_one_trajectory(times, np.array([15.0]), np.array([15.0]), 20.0, 1000.0)
    assert lab2["y_entr_tau20"].sum() == 2
