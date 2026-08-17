"""Linear model: ridge-penalized logistic regression.

Ridge rather than selection: with strongly correlated columns a sparsifying
penalty picks a different active set on every fold and turns attributions into
noise. Solver newton-cholesky: exact (full Hessian), and it leaves the
intercept unpenalized, which matters with a rare class. C is selected by log
loss, a strictly proper scoring rule sensitive to shrinkage. Class weighting
is decided by measurement: the likelihood gain is checked AFTER correcting the
prior shift.
"""

from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import log_loss
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from ..relevance import ActivityRail, NearDuplicateReporter
from ..splits import inner_rotation
from .nanpolicy import MissingEncoder


def make_pipeline(C: float, tol: float, max_iter: int, class_weight,
                  n_value_cols: int, corr_threshold: float,
                  feature_names: list[str]) -> Pipeline:
    """The whole linear arm as one unfitted pipeline over n_value_cols columns.

    C is sklearn's inverse penalty strength (λ = 1/C). n_value_cols is the width
    of the incoming descriptor matrix; the encoder doubles it, so every step
    after the encoder works on 2·n_value_cols columns and the rail must be told
    where the value block ends.

    Everything that learns from data is a step here, so a single fit call on the
    training rows of a fold fits the imputation medians, the standardization and
    the rails on those rows alone. Nothing is fitted on the pooled sample.
    """
    # Fixed step order: encoder → standardization → near-duplicates → rail → model.
    # The encoder is first because every later step assumes the [values ++
    # indicators] layout: the reporter names indicator columns by that offset
    # (relevance.NearDuplicateReporter._name), the rail exempts the indicator
    # block from its variance check, and standardization has to reach the
    # indicators too, since an unscaled 0/1 column would carry a penalty of a
    # different magnitude than the scaled values next to it. The rail sits last
    # so that a column it zeroes is zero in the matrix the model actually fits.
    return Pipeline([
        ("encode", MissingEncoder()),
        ("scale", StandardScaler()),
        ("dup", NearDuplicateReporter(threshold=corr_threshold, feature_names=feature_names)),
        # allow_nan=False is safe only here, after the encoder: the encoded
        # matrix has no NaN left, so a column that is constant is genuinely
        # uninformative rather than constant-where-measured.
        ("rail", ActivityRail(n_value_cols=n_value_cols, allow_nan=False)),
        ("clf", LogisticRegression(penalty="l2", C=C, solver="newton-cholesky",
                                   tol=tol, max_iter=max_iter, class_weight=class_weight)),
    ])


def _prior_shift_corrected_proba(pipe: Pipeline, X: np.ndarray, class_weight) -> np.ndarray:
    """Probabilities of a weighted model mapped back to the unweighted prior:
    the logit is shifted by log(w+/w−), which class weighting adds to the
    intercept. The shift is stored in pipe._prior_shift by the fitting code
    (where the class counts are known).

    Returns P(y=1) on the original class prior, one value per row of X. Without
    the correction the two weighting choices could not be compared by log loss
    at all: a balanced fit predicts on an artificially enriched prior, and its
    loss against the true labels would be worse for a reason that has nothing to
    do with how well it separates the classes.
    """
    p = pipe.predict_proba(X)[:, 1]
    if class_weight is None:
        return p
    logit = np.log(np.clip(p, 1e-12, 1 - 1e-12)) - np.log(np.clip(1 - p, 1e-12, 1 - 1e-12))
    logit = logit - getattr(pipe, "_prior_shift", 0.0)
    return 1.0 / (1.0 + np.exp(-logit))


class LinearArmResult:
    """One fitted linear arm plus everything the fold artifact has to record.

    pipeline is fitted on the whole training part of the fold at chosen_C.
    curve holds the log loss at every grid point and edge says whether the
    optimum landed on a grid end (accept.edge_fraction_max is checked against
    it). df_effective and df_at_grid_ends express the penalty as effective
    degrees of freedom, so shrinkage can be read in units of parameters rather
    than in units of C. inner_oof_scores are prior-corrected probabilities on
    the rotation's held-out rows, NaN where a rotation could not be validated,
    and are the only scores the calibrator and the threshold are allowed to see.
    """

    def __init__(self):
        self.pipeline: Pipeline | None = None
        self.chosen_C: float = np.nan
        self.class_weight = None
        self.curve: list[dict] = []
        self.edge: dict = {}
        self.df_effective: float = np.nan
        self.df_at_grid_ends: tuple[float, float] = (np.nan, np.nan)
        self.n_iter: int = 0
        self.inner_oof_scores: np.ndarray | None = None
        self.inner_oof_y: np.ndarray | None = None
        self.weight_measurement: dict = {}


def fit_linear_fold(X_tr: np.ndarray, y_tr: np.ndarray, assemblies_tr: np.ndarray,
                    feature_names: list[str], C_grid: list[float], tol: float,
                    max_iter: int, class_weight_mode: str,
                    corr_threshold: float) -> LinearArmResult:
    """Fit the linear arm on one outer fold: select C, refit, produce inner-OOF scores.

    X_tr is [n_train, p] of descriptor values with NaN for "not measured"; y_tr
    is 0/1 readiness at the horizon; assemblies_tr carries the lineage label of
    each row and is the only thing that defines the inner rotation, so rows of
    one assembly can never sit on both sides of an inner split. feature_names
    must be in the same order as the columns of X_tr — the near-duplicate report
    is written under those names.

    Raises RuntimeError if the solver reaches max_iter, since an unconverged fit
    has no meaningful coefficients or effective degrees of freedom.
    """
    res = LinearArmResult()
    n_vals = X_tr.shape[1]
    inner = inner_rotation(assemblies_tr)

    def cv_logloss(C: float, cw) -> float:
        losses = []
        for tr, va, _g in inner:
            # A rotation slot with one class on either side cannot produce a log
            # loss that says anything about separation; it is skipped rather
            # than scored, so an event-starved assembly does not vote on C.
            if len(np.unique(y_tr[va])) < 2 or len(np.unique(y_tr[tr])) < 2:
                continue
            pipe = make_pipeline(C, tol, max_iter, cw, n_vals, corr_threshold, feature_names)
            pipe.fit(X_tr[tr], y_tr[tr])
            if cw == "balanced":
                # sklearn's "balanced" weight is w_k = n / (2 n_k), so the logit
                # shift log(w+/w-) reduces to log(n_neg/n_pos). It is written in
                # the w_k form to keep the correspondence to class_weight
                # visible; changing one without the other decalibrates the arm.
                n_pos, n_neg = y_tr[tr].sum(), len(tr) - y_tr[tr].sum()
                pipe._prior_shift = float(np.log((len(tr) / (2 * n_pos))
                                                 / (len(tr) / (2 * n_neg))))
            p = _prior_shift_corrected_proba(pipe, X_tr[va], cw)
            losses.append(log_loss(y_tr[va], np.clip(p, 1e-12, 1 - 1e-12), labels=[0, 1]))
        # No usable rotation slot: return infinity so this configuration loses to
        # any configuration that could be validated, instead of winning on an
        # empty average.
        return float(np.mean(losses)) if losses else np.inf

    # Class weighting is decided by measurement: compare likelihood at the
    # mid-grid C.
    # One probe point, not the whole grid: the weighting is settled before the C
    # scan, so weighting and shrinkage are never selected jointly on the same
    # rotation, and the mid-grid C is the least extreme shrinkage on offer.
    weight_choice = None
    if class_weight_mode == "balanced":
        weight_choice = "balanced"
    elif class_weight_mode == "measure":
        c_mid = C_grid[len(C_grid) // 2]
        l_none = cv_logloss(c_mid, None)
        l_bal = cv_logloss(c_mid, "balanced")
        weight_choice = "balanced" if l_bal < l_none else None
        res.weight_measurement = {"C_probe": c_mid, "logloss_unweighted": l_none,
                                  "logloss_balanced_prior_corrected": l_bal,
                                  "chosen": str(weight_choice)}
    res.class_weight = weight_choice

    # Log loss, not average precision, selects C: it is a strictly proper scoring
    # rule, so it is minimized by the honest probability and reacts to over- and
    # under-shrinkage; a ranking measure is blind to the scale the shrinkage sets.
    for C in C_grid:
        loss = cv_logloss(C, weight_choice)
        res.curve.append({"C": C, "logloss": loss})
    losses = [c["logloss"] for c in res.curve]
    best_i = int(np.argmin(losses))
    res.chosen_C = C_grid[best_i]
    res.edge = {"C_on_edge": best_i in (0, len(C_grid) - 1),
                "position": best_i, "grid_len": len(C_grid), "exception2": None}

    pipe = make_pipeline(res.chosen_C, tol, max_iter, weight_choice, n_vals,
                         corr_threshold, feature_names)
    pipe.fit(X_tr, y_tr)
    res.pipeline = pipe
    clf: LogisticRegression = pipe.named_steps["clf"]
    res.n_iter = int(np.max(clf.n_iter_))
    if res.n_iter >= max_iter:
        raise RuntimeError(f"linear model hit the iteration ceiling {max_iter} — "
                           "raise the ceiling instead of silencing the warning")

    # df(λ) in comparable units: singular values of the WEIGHTED matrix
    # W^{1/2}·Z at convergence on the training fold; λ = 1/C in sklearn's form
    # (penalty ½‖w‖² with C·Σloss). Diagnostic: df must vary across the grid.
    # W = p(1−p) are the IRLS working weights of logistic regression at the
    # converged fit, which is what makes df comparable across C: the same
    # ridge trace formula as in the Gaussian case, applied to the weighted
    # design (Hastie, Tibshirani & Friedman, "The Elements of Statistical
    # Learning", 2nd ed., §3.4.1 and §4.4.1).
    Z = pipe[:-1].transform(X_tr)
    p_hat = pipe.predict_proba(X_tr)[:, 1]
    W = p_hat * (1 - p_hat)
    d = np.linalg.svd(np.sqrt(W)[:, None] * Z, compute_uv=False)

    def df_of(C):
        # trace of the ridge hat matrix: sum over singular values d of
        # d²/(d²+λ), which runs from p at λ→0 down to 0 at λ→∞.
        lam = 1.0 / C
        return float(np.sum(d ** 2 / (d ** 2 + lam)))
    res.df_effective = df_of(res.chosen_C)
    res.df_at_grid_ends = (df_of(C_grid[0]), df_of(C_grid[-1]))
    # An optimum on the edge of a PHYSICALLY meaningful range is recorded
    # explicitly; the grid has nothing to extend into. Lower edge: the fit has
    # degenerated to a constant (df≈0, so decreasing C further changes nothing).
    # Upper edge: the unpenalized optimum is reached (df at the edge is
    # indistinguishable from the chosen one).
    if res.edge["C_on_edge"]:
        if res.edge["position"] == 0 and res.df_effective < 0.5:
            res.edge["exception2"] = ("lower edge: fit degenerated to a constant "
                                      f"(df={res.df_effective:.2f}); edge is meaningful")
        elif (res.edge["position"] == len(C_grid) - 1
              and abs(res.df_effective - res.df_at_grid_ends[1]) < 0.1):
            res.edge["exception2"] = ("upper edge: unpenalized optimum reached "
                                      f"(df={res.df_effective:.2f} saturated); "
                                      "edge is meaningful")

    # Inner-OOF at the chosen parameters — feeds the calibrator and the threshold.
    # A separate pass, and not the fits made during selection: those were made at
    # every C in the grid, and a calibrator fitted on them would inherit the
    # selection. Rows left NaN (a rotation slot that could not be fitted) are
    # dropped downstream by the calibrator's own finiteness mask.
    oof_s = np.full(len(y_tr), np.nan)
    for tr, va, _g in inner:
        if len(np.unique(y_tr[tr])) < 2:
            continue
        pp = make_pipeline(res.chosen_C, tol, max_iter, weight_choice, n_vals,
                           corr_threshold, feature_names)
        pp.fit(X_tr[tr], y_tr[tr])
        if weight_choice == "balanced":
            n_pos, n_neg = y_tr[tr].sum(), len(tr) - y_tr[tr].sum()
            pp._prior_shift = float(np.log((len(tr) / (2 * n_pos)) / (len(tr) / (2 * n_neg))))
        oof_s[va] = _prior_shift_corrected_proba(pp, X_tr[va], weight_choice)
    res.inner_oof_scores = oof_s
    res.inner_oof_y = y_tr.copy()
    return res


def linear_contributions(pipe: Pipeline, X: np.ndarray) -> np.ndarray:
    """Per-frame additive contributions IN MODEL SPACE: coef·z after
    preprocessing; (value, indicator) pairs are summed back to the source column.

    Returns [n_rows, p] in logit units on the source columns of X, in the order
    the pipeline was fitted with. The terms sum to the logit minus the intercept,
    which is what makes them comparable with the TreeSHAP contributions of the
    boosted arm; the intercept is deliberately not distributed over columns.
    """
    Z = pipe[:-1].transform(X)
    clf: LogisticRegression = pipe.named_steps["clf"]
    contrib = Z * clf.coef_[0][None, :]
    enc: MissingEncoder = pipe.named_steps["encode"]
    return enc.restore_column_contributions(contrib)
