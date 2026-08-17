"""PLS-FMA benchmark arm: functional mode analysis by partial least squares.

Implements the method of Krivobokova, Briones, Hub, Munk & de Groot (2012),
"Partial least-squares functional mode analysis: application to the membrane
proteins AQP1, Aqy1, and CLC-ec1", Biophys. J. 103:786-796,
doi:10.1016/j.bpj.2012.07.022, following the authors' own implementation of
record: the g_fma trajectory-analysis module distributed by the de Groot lab
(g_fma-beta, fma.cpp / partial_least_squares.cpp by J.H. Peters). The core
here is a line-for-line port of that tool's PLS routine:

  * X and f are mean-centered by the TRAINING means only, with no variance
    scaling of either (the PLS objective maximizes covariance, so column
    standardization would change the extracted mode);
  * Helland (1988) / Denham (1995) PLS: w_j = X^T r, scores t_j = X w_j,
    regression of f on the accumulated scores by least squares, deflation of
    the RESIDUAL of f only (X is never deflated); beta = W q;
  * the number of latent components is scanned exhaustively (the tool's
    -optDim scans 1..dim) and chosen by the Pearson correlation between model
    and data on held-out trajectories, the overfitting control the paper
    makes mandatory;
  * two mode vectors are recorded: the MCM (beta, the regression coefficient
    vector, k-dependent) and the ewMCM (the scaled first weight vector
    w1 = X^T f, SI Eq. 4), which is the basis-independent collective mode the
    paper visualizes.

Declared deviations from the published applications:

  (a) validation is a rotation over whole trajectory units instead of the
      authors' single chronological split, keeping the same independence
      requirement (never shuffled frames; hAQP1 was already validated across
      monomers) and applying it to every unit in turn;
  (b) the functional quantity here is the binary readiness label, whereas
      every published application used a continuous scalar; with a binary f
      the machinery acts as a linear discriminant on coordinates, and the
      results are reported on that basis;
  (c) the input is channel C-alpha coordinates (community practice per the
      authors' FMA pages); the paper's applications used backbone or
      heavy-atom sets.

tools/plsfma_selfcheck.py runs the same core on a continuous structural
quantity of this project's systems, the regime the method was published for.
"""

from __future__ import annotations

import numpy as np

from ..splits import inner_rotation


class PlsDenhamModel:
    """Centered linear model f ≈ y_mean + (X - x_mean) @ beta.

    x_mean and beta have length p = 3·(number of C-alphas): the coordinate
    vector is flat, three Cartesian components per atom in A, in the column order
    of the coords artifact (resid-major, then x, y, z). beta is the MCM, the
    collective mode in that same basis, so an entry of it can be read back as the
    displacement of one C-alpha along one axis.
    weights is [p, k] with the component weight vectors w_j as columns and q is
    the length-k regression of f on the accumulated scores; beta = W q.
    """

    def __init__(self, x_mean: np.ndarray, y_mean: float, beta: np.ndarray,
                 weights: np.ndarray, q: np.ndarray):
        self.x_mean = x_mean
        self.y_mean = y_mean
        self.beta = beta
        self.weights = weights  # columns w_j, in component order
        self.q = q

    def predict(self, X: np.ndarray) -> np.ndarray:
        """Model value of the functional quantity per row of X; an unbounded real.

        X must be superimposed onto the same reference as the training
        coordinates and carry the same columns in the same order: beta is defined
        in that frame, and an unaligned structure would be scored against a rigid
        rotation rather than against the mode.
        """
        return self.y_mean + (np.asarray(X, dtype=np.float64) - self.x_mean) @ self.beta


def pls_denham(X_tr: np.ndarray, y_tr: np.ndarray, k: int) -> PlsDenhamModel:
    """The g_fma PLS core: Helland/Denham with deflation of f's residual only.

    X_tr is [n, p] of superimposed Cartesian coordinates, y_tr the functional
    quantity per frame, k the number of latent components requested. Fewer than k
    are used if the residual runs out first (see the break below), which is not an
    error: the extra components would carry no covariance with f.

    Only f is deflated, X never is. This is the Helland (1988) / Denham (1995)
    form the authors implement; it keeps every w_j expressed in the original
    coordinate basis, which is what makes beta readable as a displacement of the
    channel rather than as a coefficient on a deflated space.
    """
    X = np.asarray(X_tr, dtype=np.float64)
    y = np.asarray(y_tr, dtype=np.float64)
    x_mean = X.mean(axis=0)
    y_mean = float(y.mean())
    Xc = X - x_mean
    yc = y - y_mean
    W, T = [], []
    r = yc
    q = np.zeros(0)
    for _j in range(k):
        # w_j is the covariance of every coordinate with what is left of f: the
        # direction along which the remaining functional variation is carried.
        w = Xc.T @ r
        # Unit-normalizing each weight vector leaves the model invariant
        # (beta = W q absorbs any per-column scale) and keeps the score matrix
        # numerically solvable once the residual — and with it w — has shrunk
        # by many orders of magnitude.
        nw = float(np.linalg.norm(w))
        if not np.isfinite(nw) or nw == 0.0:
            break
        w = w / nw
        W.append(w)
        T.append(Xc @ w)
        Tm = np.column_stack(T)
        # f is regressed on ALL scores accumulated so far, not on the newest one
        # alone: the scores are not orthogonal when X is left undeflated, so a
        # per-component regression would count shared variance twice.
        q, *_ = np.linalg.lstsq(Tm, yc, rcond=None)
        r = yc - Tm @ q
    beta = np.column_stack(W) @ q
    return PlsDenhamModel(x_mean, y_mean, beta, np.column_stack(W), q)


def ew_mcm(X_tr: np.ndarray, y_tr: np.ndarray) -> np.ndarray:
    """Ensemble-weighted MCM (SI Eq. 4): the scaled first weight vector.

    ewMCM ∝ w1 = X_c^T f_c: per-coordinate covariance with the functional
    quantity; unlike the MCM it does not depend on the chosen component count.
    The SI scale factor (f^T X X^T f) / (n · var(f_hat)) is applied so the
    vector carries the paper's units.

    Returns a length-p vector in the coordinate basis of X_tr, one entry per
    Cartesian component. Being independent of the component count, it is the
    vector worth comparing between folds: the MCM of a two-component fit and of a
    nine-component fit are not the same object.
    """
    X = np.asarray(X_tr, dtype=np.float64)
    y = np.asarray(y_tr, dtype=np.float64)
    Xc = X - X.mean(axis=0)
    yc = y - y.mean()
    w1 = Xc.T @ yc
    t1 = Xc @ w1
    denom = len(y) * float(np.var(t1 * (float(t1 @ yc) / float(t1 @ t1))))
    scale = float(yc @ t1) / denom if denom > 0 else 1.0
    return w1 * scale


class PlsFmaFoldResult:
    """The fitted PLS-FMA arm of one fold plus what the fold artifact records.

    curve holds the validation correlation at every admissible component count
    and edge whether the winner sat at a grid end. mode_vector is the MCM (beta of
    the chosen fit) and ewmcm_vector the ensemble-weighted mode; only their norms
    are written to the artifact, the vectors themselves being p-dimensional.
    inner_oof_scores are RAW model values, not probabilities, which is why the
    calibrator for this arm is built with kind="raw".
    """

    def __init__(self):
        self.model: PlsDenhamModel | None = None
        self.chosen_components: int = 0
        self.curve: list[dict] = []
        self.edge: dict = {}
        self.mode_vector: np.ndarray | None = None
        self.ewmcm_vector: np.ndarray | None = None
        self.inner_oof_scores: np.ndarray | None = None
        self.inner_oof_y: np.ndarray | None = None


def fit_plsfma_fold(C_tr: np.ndarray, y_tr: np.ndarray, assemblies_tr: np.ndarray,
                    components_grid: list[int]) -> PlsFmaFoldResult:
    """Fit the PLS-FMA arm on one outer fold: scan the component count, refit, score OOF.

    C_tr is [n_train, p] of superimposed C-alpha coordinates in the row order of
    the model matrix (see models._aligned_coords), y_tr the functional quantity,
    and assemblies_tr the lineage labels defining the rotation. The coordinate
    matrix must be complete: this arm has no missing-value branch, and a NaN would
    propagate through the covariance into every component.

    Raises RuntimeError when no component count in components_grid is admissible.
    """
    res = PlsFmaFoldResult()
    inner = inner_rotation(assemblies_tr)
    n, p = C_tr.shape
    # k < min(n, p) is the rank bound: beyond it the score matrix is singular and
    # the extra components fit residual noise. The floor of 1 is the paper's own
    # lower limit; a zero-component model is the mean of f.
    grid = sorted({k for k in components_grid if 1 <= k < min(n, p)})
    if not grid:
        raise RuntimeError("no admissible PLS component count for this fold")

    def cv_corr(k: int) -> float:
        """Mean Pearson r between model and data over the rotation's held-out units."""
        rs = []
        for tr, va, _g in inner:
            if len(np.unique(y_tr[va])) < 2 or len(np.unique(y_tr[tr])) < 2:
                continue
            pred = pls_denham(C_tr[tr], y_tr[tr], k).predict(C_tr[va])
            # A constant prediction has no correlation defined (zero variance);
            # such a slot is dropped rather than counted as r = 0.
            if np.std(pred) == 0:
                continue
            rs.append(float(np.corrcoef(pred, y_tr[va])[0, 1]))
        # Nothing validated: -inf, so this component count cannot win the argmax.
        return float(np.mean(rs)) if rs else -np.inf

    for k in grid:
        res.curve.append({"components": k, "cv_pearson_r": cv_corr(k)})
    scores = [c["cv_pearson_r"] for c in res.curve]
    # The authors' -optDim rule: the component count with the best validation
    # correlation. The full curve is recorded so the paper's parsimony reading
    # ("lowest dimensionality with adequate predictive power") stays checkable.
    best_i = int(np.argmax(scores))
    res.chosen_components = grid[best_i]
    res.edge = {"on_edge": best_i in (0, len(grid) - 1) and len(grid) > 1,
                "position": best_i, "grid_len": len(grid), "exception2": None}
    if res.edge["on_edge"] and best_i == 0 and grid[0] == 1:
        res.edge["exception2"] = ("lower bound of the admissible component count "
                                  "(one latent component); fewer do not exist")
    # Refit on the whole training part at the chosen component count; the mode
    # vectors are taken from this fit, so they rest on every training trajectory
    # rather than on one rotation slot.
    res.model = pls_denham(C_tr, y_tr, res.chosen_components)
    res.mode_vector = res.model.beta.copy()
    res.ewmcm_vector = ew_mcm(C_tr, y_tr)

    # Inner-OOF at the chosen component count — feeds the calibrator and the
    # threshold, and is fitted again per slot so no held-out frame contributed to
    # the coefficients that score it.
    oof = np.full(n, np.nan)
    for tr, va, _g in inner:
        if len(np.unique(y_tr[tr])) < 2:
            continue
        oof[va] = pls_denham(C_tr[tr], y_tr[tr], res.chosen_components).predict(C_tr[va])
    res.inner_oof_scores = oof
    res.inner_oof_y = y_tr.copy()
    return res


def plsfma_predict(model: PlsDenhamModel, C: np.ndarray) -> np.ndarray:
    """Raw model value per frame; higher means closer to the f=1 end of the fit.

    Not a probability and not bounded to [0, 1]: with the binary readiness label
    the regression can predict outside the label range. Ranking measures use it as
    it is; a probability scale only exists after the raw-mode calibrator.
    """
    return model.predict(C)
