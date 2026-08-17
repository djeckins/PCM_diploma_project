"""The PLS core must be the authors' algorithm, not merely something PLS-like.

Three properties pin it down. The Helland/Denham recursion is mathematically
equivalent to NIPALS PLS1 for a univariate f, so predictions must agree with
scikit-learn's NIPALS to numerical precision at every component count. The
first weight vector must be exactly X_c^T f_c — the covariance direction the
ewMCM is built from. And centering must use training means only, so a shifted
validation set is handled through the stored means, not recentered on itself.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pcm2.models.plsfma import ew_mcm, fit_plsfma_fold, pls_denham  # noqa: E402


def _data(n=120, p=15, seed=7):
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(n, p))
    beta = rng.normal(size=p)
    y = X @ beta + 0.3 * rng.normal(size=n)
    return X, y


@pytest.mark.parametrize("k", [1, 2, 3, 5, 8])
def test_matches_sklearn_nipals_for_univariate_f(k):
    sklearn = pytest.importorskip("sklearn.cross_decomposition")
    X, y = _data()
    ours = pls_denham(X, y, k).predict(X)
    ref = sklearn.PLSRegression(n_components=k, scale=False).fit(X, y).predict(X).ravel()
    assert np.allclose(ours, ref, atol=1e-8), (
        f"k={k}: max deviation {np.abs(ours - ref).max():.3g}")


def test_first_weight_vector_is_covariance_direction():
    X, y = _data()
    m = pls_denham(X, y, 3)
    Xc = X - X.mean(axis=0)
    yc = y - y.mean()
    w1 = Xc.T @ yc
    cos = w1 @ m.weights[:, 0] / (np.linalg.norm(w1) * np.linalg.norm(m.weights[:, 0]))
    assert cos > 1 - 1e-12
    ew = ew_mcm(X, y)
    cos_ew = w1 @ ew / (np.linalg.norm(w1) * np.linalg.norm(ew))
    assert cos_ew > 1 - 1e-12


def test_centering_uses_training_means_only():
    X, y = _data()
    m = pls_denham(X, y, 4)
    shift = 5.0
    assert np.allclose(m.predict(X + shift),
                       m.predict(X) + shift * m.beta.sum(), atol=1e-8)


def test_component_selection_on_held_out_units():
    rng = np.random.default_rng(3)
    n_per, p = 80, 12
    units, Xs, ys = [], [], []
    beta = rng.normal(size=p)
    for u in range(4):
        X = rng.normal(size=(n_per, p))
        y = (X @ beta + 0.5 * rng.normal(size=n_per) > 0).astype(int)
        Xs.append(X); ys.append(y); units += [f"u{u}"] * n_per
    res = fit_plsfma_fold(np.vstack(Xs), np.concatenate(ys),
                          np.array(units), list(range(1, 11)))
    ks = [c["components"] for c in res.curve]
    assert ks == list(range(1, 11)), "exhaustive scan expected"
    assert res.chosen_components == ks[int(np.argmax([c["cv_pearson_r"] for c in res.curve]))]
    assert res.ewmcm_vector is not None and res.mode_vector is not None
    oof = res.inner_oof_scores
    assert np.isfinite(oof).all() and np.corrcoef(oof, np.concatenate(ys))[0, 1] > 0.3
