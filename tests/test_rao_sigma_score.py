"""The published Rao structure score must follow the authors' own rules.

Three properties are checked against the author scripts' behaviour: the
nearest-node energy lookup must equal a brute-force minimization over all
grid nodes (their R/Python rule); a clearly hydrophilic-wide point must
contribute nothing while a hydrophobic-narrow point must be flagged with a
positive contour distance; and the score must be the plain sum over flagged
points, insensitive to unflagged company.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pcm2.baselines.published import (RAO_DEWETTING_THRESHOLD_KJMOL,  # noqa: E402
                                      _grid, rao_nearest_energy, rao_sigma_d)


def test_nearest_node_equals_brute_force():
    h_ax, r_ax, E = _grid()
    rng = np.random.default_rng(11)
    h = rng.uniform(h_ax[0] - 0.1, h_ax[-1] + 0.1, size=200)
    r = rng.uniform(r_ax[0] - 0.05, r_ax[-1] + 0.05, size=200)
    ours = rao_nearest_energy(h, r)
    H, R = np.meshgrid(h_ax, r_ax, indexing="ij")
    nodes = np.stack([H.ravel(), R.ravel()], axis=1)
    flat = E.ravel()
    brute = np.array([flat[np.argmin((nodes[:, 0] - hh) ** 2
                                     + (nodes[:, 1] - rr) ** 2)]
                      for hh, rr in zip(h, r)])
    assert np.array_equal(ours, brute)


def test_wet_side_contributes_nothing():
    # A wide, hydrophilic pore point: energy far below the contour.
    sigma, n = rao_sigma_d(np.array([-0.4]), np.array([0.55]))
    assert sigma == 0.0 and n == 0


def test_hydrophobic_narrow_point_is_flagged():
    sigma, n = rao_sigma_d(np.array([0.25]), np.array([0.12]))
    assert n == 1 and sigma > 0
    e = rao_nearest_energy(np.array([0.25]), np.array([0.12]))[0]
    assert e > RAO_DEWETTING_THRESHOLD_KJMOL


def test_score_is_a_sum_over_flagged_points_only():
    h = np.array([0.25, 0.25, -0.4])
    r = np.array([0.12, 0.12, 0.55])
    sigma3, n3 = rao_sigma_d(h, r)
    sigma1, n1 = rao_sigma_d(h[:1], r[:1])
    assert n3 == 2 and n1 == 1
    assert np.isclose(sigma3, 2 * sigma1)


def test_nan_inputs_are_unknown_not_zero():
    sigma, n = rao_sigma_d(np.array([np.nan]), np.array([np.nan]))
    assert np.isnan(sigma) and n == 0
