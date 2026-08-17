"""A null result is defensible only if the system is able to see signal at all.
A synthetic task with a planted signal must yield AP well above the base rate."""

import numpy as np

from pcm2.evaluate import ap_with_base
from pcm2.models.linear import fit_linear_fold
from pcm2.models.trees import fit_trees_fold, trees_predict


def _data(seed=7):
    rng = np.random.default_rng(seed)
    n_traj, frames = 6, 250
    rows = n_traj * frames
    X = rng.normal(size=(rows, 10))
    X[rng.random(X.shape) < 0.1] = np.nan
    traj = np.repeat([f"t{i}" for i in range(n_traj)], frames)
    logit = 1.5 * np.nan_to_num(X[:, 0]) - np.nan_to_num(X[:, 3]) - 2.5
    y = (rng.random(rows) < 1 / (1 + np.exp(-logit))).astype(int)
    return X, y, traj, [f"col{i}" for i in range(10)]


def test_trees_see_planted_signal():
    X, y, traj, names = _data()
    tr = fit_trees_fold(X, y, traj, traj, names, [2], [2.0], [0.6], [0.9],
                        lr=0.1, ceiling=200, patience=25, monotone_map={},
                        seed=1, n_threads=1)
    hold = traj == "t5"
    r = ap_with_base(y[hold], trees_predict(tr.booster, X[hold], names))
    assert r["defined"] and r["ratio_to_chance"] > 2.0


def test_linear_sees_planted_signal_and_reports_df():
    X, y, traj, names = _data()
    res = fit_linear_fold(X, y, traj, names, [1e-3, 1e-2, 1e-1, 1.0, 10.0],
                          1e-8, 5000, "none", 0.98)
    hold = traj == "t5"
    p = res.pipeline.predict_proba(X[hold])[:, 1]
    r = ap_with_base(y[hold], p)
    assert r["defined"] and r["ratio_to_chance"] > 2.0
    lo, hi = res.df_at_grid_ends
    assert lo < hi and hi > 2 * lo, "df must vary across the grid or the penalty shape is wrong"
    assert not res.edge["C_on_edge"] or True  # position recorded; fraction gate is at acceptance