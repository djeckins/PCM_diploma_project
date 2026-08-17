"""Boosting: xgboost with a wrapper that derives the leaf threshold per fold.

Three task properties → three constraints: the effective sample size is
trajectories (shallow depth, strict leaf); positive frames are grouped into
events (the leaf threshold is parameterized in EVENTS, k_events, and is
computed inside the fold); descriptors are correlated by construction
(column subsampling is tuned).

The learning rate is FIXED; the number of rounds is set by early stopping
with rotation over whole trajectories; every grid point gets its OWN round
count M(g) = median over the rotation; the inner CV compares grid points at
their own M(g); the winner is refit at its own M(g).

base_score is pinned to the base rate of the training rows (and recorded)
instead of the library default, and h0 = p0(1−p0) is derived from it.
"""

from __future__ import annotations

import numpy as np
import xgboost as xgb

from ..runtime import derive_seed
from ..splits import inner_rotation
from .nanpolicy import passthrough_nan_check


def count_events(y: np.ndarray, traj: np.ndarray) -> int:
    """An event is a run of consecutive positive frames within one trajectory.

    Counts the rising edges of y within each trajectory separately and returns at
    least 1. Positive frames are not independent: the label marks every frame
    whose horizon contains the same crossing, so one crossing produces a block of
    consecutive positives. The number of blocks, not the number of positive
    frames, is the count of independent positive observations in the fold.

    Counting per trajectory matters because the concatenated array has no
    boundary marker: a run ending at the last frame of one trajectory and one
    starting at the first frame of the next would merge into a single event.
    The floor of 1 keeps the value usable as the denominator in leaf_threshold.
    """
    n = 0
    for t in np.unique(traj):
        yy = y[traj == t]
        n += int(np.sum((yy == 1) & (np.concatenate([[0], yy[:-1]]) == 0)))
    return max(n, 1)


def leaf_threshold(k_events: float, y_tr: np.ndarray, traj_tr: np.ndarray,
                  class_weight: float = 1.0) -> tuple[float, dict]:
    """threshold = k_events · (positive frames/event on the fold) · h0 · fold weight.

    Returns (min_child_weight for xgboost, the arithmetic behind it for the fold
    artifact). For binary:logistic, min_child_weight is a bound on the sum of
    Hessians of the rows in a leaf, and each row's Hessian is p(1−p) evaluated at
    the current prediction; at the pinned starting value p0 every row contributes
    h0 = p0(1−p0). Multiplying by the mean number of positive frames per event
    converts a requirement stated in EVENTS into that Hessian mass, so k_events
    means the same thing in a fold with 3 crossings and in a fold with 40 —
    a raw threshold would not, and the leaf constraint would drift with the base
    rate rather than with the physics.

    class_weight scales the Hessians when instance weights are in use and so has
    to scale the bound with them; with unweighted rows it is 1.
    """
    n_pos = float(np.sum(y_tr))
    n_events = count_events(y_tr, traj_tr)
    frames_per_event = n_pos / n_events if n_events else 1.0
    # p0 is also the value handed to xgboost as base_score. Clipped away from 0
    # and 1 because a fold with no positive frame would give h0 = 0 and a leaf
    # bound of 0, which admits leaves of a single frame.
    p0 = float(np.clip(np.mean(y_tr), 1e-6, 1 - 1e-6))
    h0 = p0 * (1 - p0)
    thr = k_events * frames_per_event * h0 * class_weight
    return thr, {"n_pos": n_pos, "n_events": n_events,
                 "frames_per_event": frames_per_event, "base_score_p0": p0, "h0": h0}


def _params(depth: int, colsample: float, subsample: float, lr: float,
            min_child_weight: float, base_score: float, seed: int,
            n_threads: int, monotone: tuple[int, ...] | None) -> dict:
    """The xgboost parameter dict for one fit; the only place these keys are set.

    monotone must be one sign per column, in the exact order of the feature_names
    the DMatrix is built with, since xgboost reads the constraint string
    positionally. A tuple misaligned with the columns would silently constrain
    the wrong descriptor.
    """
    p = {"objective": "binary:logistic", "eval_metric": "logloss",
         "max_depth": depth, "eta": lr, "colsample_bytree": colsample,
         "subsample": subsample, "min_child_weight": min_child_weight,
         "base_score": base_score, "seed": seed, "nthread": n_threads,
         "tree_method": "hist"}
    if monotone is not None:
        p["monotone_constraints"] = "(" + ",".join(str(s) for s in monotone) + ")"
    return p


class TreesFoldResult:
    """One fitted booster plus everything the fold artifact has to record.

    chosen is the winning grid point and M the number of rounds it was refit
    with. grid_scores keeps every grid point with its own M and its validation
    loss, so the selection can be re-read without refitting. edges flags a
    parameter that won at a grid end, leaf_math the arithmetic behind
    min_child_weight, hit_ceiling that early stopping never triggered (the round
    count is then a ceiling artifact, not a chosen value), and blocked_inner that
    the time-block fallback replaced the trajectory rotation. inner_oof_scores
    are probabilities on the rotation's held-out rows, NaN where a rotation could
    not be fitted, and are the only scores the calibrator may see.
    """

    def __init__(self):
        self.booster: xgb.Booster | None = None
        self.chosen: dict = {}
        self.rounds_by_rotation: dict = {}
        self.M: int = 0
        self.grid_scores: list[dict] = []
        self.edges: dict = {}
        self.leaf_math: dict = {}
        self.inner_oof_scores: np.ndarray | None = None
        self.inner_oof_y: np.ndarray | None = None
        self.hit_ceiling: bool = False
        self.blocked_inner: bool = False


def _rotation_viable(inner, y: np.ndarray) -> bool:
    """True if at least one rotation slot has both classes on both sides.

    Early stopping needs a validation set that contains positives; one usable
    slot is the minimum for the rotation to say anything, and its absence is what
    triggers the time-block fallback.
    """
    return any(len(np.unique(y[va])) > 1 and len(np.unique(y[tr])) > 1
               for tr, va, _ in inner)


def _blocked_inner(traj_tr: np.ndarray, times_tr: np.ndarray, y_tr: np.ndarray,
                   embargo_ps: float, n_blocks: int = 2):
    """Time-blocked inner rotation for event-starved training pools.

    When fewer than two training trajectories carry both classes, the
    trajectory-level rotation cannot validate. Each event-carrying trajectory
    is then split into contiguous time blocks with the two-sided embargo of
    splits.time_block_folds (feature window + label horizon), and each block
    in turn serves as validation against everything outside its embargo —
    standard blocked time-series validation, outer lineage folds untouched.

    Returns the same (train_idx, val_idx, tag) triples as splits.inner_rotation,
    with indices into the training arrays passed in, so the callers below need no
    branch of their own.
    """
    from ..splits import time_block_folds
    rotations = []
    for t in np.unique(traj_tr):
        m = np.flatnonzero(traj_tr == t)
        if len(np.unique(y_tr[m])) < 2:
            continue
        # Frames of the other trajectories always stay in training: the embargo
        # only has to separate frames that could share a feature window or a
        # label horizon, and those are contiguous in time within one trajectory.
        others = np.flatnonzero(traj_tr != t)
        for b, (btr, bva) in enumerate(
                time_block_folds(times_tr[m], n_blocks, embargo_ps)):
            tr = np.concatenate([m[btr], others])
            rotations.append((tr, m[bva], f"{t}#block{b}"))
    return rotations


def fit_trees_fold(X_tr: np.ndarray, y_tr: np.ndarray, assemblies_tr: np.ndarray,
                   traj_tr: np.ndarray, feature_names: list[str],
                   depth_grid: list[int], k_events_grid: list[float],
                   colsample_grid: list[float], subsample_grid: list[float],
                   lr: float, ceiling: int, patience: int,
                   monotone_map: dict[str, int], seed: int,
                   n_threads: int, times_tr: np.ndarray | None = None,
                   block_embargo_ps: float | None = None) -> TreesFoldResult:
    """Fit the boosted arm on one outer fold: select the grid point, refit, score OOF.

    X_tr is [n_train, p] with NaN kept as NaN (the branch direction for a missing
    value is learned, see nanpolicy). assemblies_tr defines the inner rotation and
    traj_tr the event counting and the time blocks; the two differ whenever one
    assembly contributed several trajectories. monotone_map is keyed by column
    name, so a column absent from feature_names contributes no constraint at all.
    times_tr (ps) and block_embargo_ps are only consulted for the time-block
    fallback and may be None when the trajectory rotation is known to work.

    Raises RuntimeError if no grid point could be validated on any rotation slot,
    which is how an event-starved fold is refused instead of being fitted blind.
    """
    X_tr = passthrough_nan_check(X_tr)
    res = TreesFoldResult()
    monotone = None
    if monotone_map:
        # Columns not named in the config get sign 0, "unconstrained": only signs
        # that follow from physics are declared, and the constraint tuple has to
        # cover every column regardless.
        monotone = tuple(monotone_map.get(f, 0) for f in feature_names)
    inner = inner_rotation(assemblies_tr)
    if (not _rotation_viable(inner, y_tr) and times_tr is not None
            and block_embargo_ps is not None):
        inner = _blocked_inner(traj_tr, times_tr, y_tr, block_embargo_ps)
        res.blocked_inner = True
    grid = [{"depth": d, "k_events": k, "colsample": c, "subsample": s}
            for d in depth_grid for k in k_events_grid
            for c in colsample_grid for s in subsample_grid]

    def one_rotation_fit(g: dict, tr: np.ndarray, va: np.ndarray, tag: str):
        """One grid point on one rotation slot: (booster, per-round val loss, leaf math)."""
        # The leaf bound and the starting probability are recomputed from this
        # slot's training rows only. Taking them from the whole fold would let
        # the base rate of the held-out trajectory into the fit.
        thr, math = leaf_threshold(g["k_events"], y_tr[tr], traj_tr[tr])
        params = _params(g["depth"], g["colsample"], g["subsample"], lr, thr,
                         math["base_score_p0"], derive_seed(seed, tag), n_threads,
                         monotone)
        dtr = xgb.DMatrix(X_tr[tr], label=y_tr[tr], feature_names=feature_names,
                          missing=np.nan)
        dva = xgb.DMatrix(X_tr[va], label=y_tr[va], feature_names=feature_names,
                          missing=np.nan)
        hist: dict = {}
        bst = xgb.train(params, dtr, num_boost_round=ceiling,
                        evals=[(dva, "val")], early_stopping_rounds=patience,
                        evals_result=hist, verbose_eval=False)
        losses = hist["val"]["logloss"]
        return bst, losses, math

    best_g, best_score, best_M = None, np.inf, 0
    for gi, g in enumerate(grid):
        rounds, histories = [], []
        for tr, va, gname in inner:
            if len(np.unique(y_tr[va])) < 2 or len(np.unique(y_tr[tr])) < 2:
                continue
            bst, losses, _ = one_rotation_fit(g, tr, va, f"rot-{gi}-{gname}")
            rounds.append(int(bst.best_iteration) + 1)
            histories.append(losses)
        if not rounds:
            continue
        M = int(np.median(rounds))  # median: one long rotation must not set M
        # Every grid point is compared at ITS OWN M: each rotation's loss is read
        # at round M, or at its last recorded round if it stopped earlier. Reading
        # all grid points at one common M would confound the number of rounds with
        # the parameters, since with a fixed learning rate the two trade off.
        score = float(np.mean([h[min(M, len(h)) - 1] for h in histories]))
        res.grid_scores.append({**g, "M": M, "score": score,
                                "rounds_by_rotation": rounds})
        if score < best_score:
            best_g, best_score, best_M = g, score, M
    if best_g is None:
        raise RuntimeError("no rotation produced both classes — degenerate fold")
    res.chosen = best_g
    res.M = best_M
    res.hit_ceiling = best_M >= ceiling
    # A parameter that won at a grid end is recorded, not corrected here: the
    # accept step compares the fraction of folds doing so against
    # accept.edge_fraction_max. A grid of one value has no interior, so it cannot
    # count as an edge.
    res.edges = {
        "depth": best_g["depth"] in (min(depth_grid), max(depth_grid)) and len(depth_grid) > 1,
        "k_events": best_g["k_events"] in (min(k_events_grid), max(k_events_grid)) and len(k_events_grid) > 1,
        "colsample": best_g["colsample"] in (min(colsample_grid), max(colsample_grid)) and len(colsample_grid) > 1,
        "subsample": best_g["subsample"] in (min(subsample_grid), max(subsample_grid)) and len(subsample_grid) > 1,
    }
    # The winner is refit on all training trajectories at its own M(g). The refit
    # gets no validation set, so the round count has to come from the rotation;
    # holding one trajectory back for early stopping here would spend a whole
    # independent unit on a number that is already measured.
    thr, math = leaf_threshold(best_g["k_events"], y_tr, traj_tr)
    res.leaf_math = {**math, "min_child_weight": thr, "k_events": best_g["k_events"],
                     "base_score_note": "initial value pinned explicitly to the fold "
                                        "base rate and recorded"}
    params = _params(best_g["depth"], best_g["colsample"], best_g["subsample"], lr,
                     thr, math["base_score_p0"], derive_seed(seed, "final"), n_threads,
                     monotone)
    dtr = xgb.DMatrix(X_tr, label=y_tr, feature_names=feature_names, missing=np.nan)
    res.booster = xgb.train(params, dtr, num_boost_round=res.M, verbose_eval=False)

    # Inner-OOF at the chosen parameters — feeds the calibrator and the threshold.
    # Refit from scratch on each rotation slot rather than reusing the boosters
    # from the selection loop: those were fitted with early stopping ON the slot
    # they are predicting, so their held-out scores are not held out at all.
    # A distinct seed per slot keeps the column and row subsampling independent
    # across slots while staying reproducible from the one recorded seed.
    oof = np.full(len(y_tr), np.nan)
    for tr, va, gname in inner:
        if len(np.unique(y_tr[tr])) < 2:
            continue
        thr_i, math_i = leaf_threshold(best_g["k_events"], y_tr[tr], traj_tr[tr])
        p_i = _params(best_g["depth"], best_g["colsample"], best_g["subsample"], lr,
                      thr_i, math_i["base_score_p0"], derive_seed(seed, f"oof-{gname}"),
                      n_threads, monotone)
        d_i = xgb.DMatrix(X_tr[tr], label=y_tr[tr], feature_names=feature_names,
                          missing=np.nan)
        b_i = xgb.train(p_i, d_i, num_boost_round=res.M, verbose_eval=False)
        oof[va] = b_i.predict(xgb.DMatrix(X_tr[va], feature_names=feature_names,
                                          missing=np.nan))
    res.inner_oof_scores = oof
    res.inner_oof_y = y_tr.copy()
    return res


def trees_predict(booster: xgb.Booster, X: np.ndarray,
                  feature_names: list[str]) -> np.ndarray:
    """P(y=1) per row of X, on the booster's own probability scale (uncalibrated).

    feature_names must be the same list, in the same order, that the booster was
    trained with: the DMatrix is checked against the stored names, and the
    monotone constraints are positional. missing=np.nan repeats the training
    convention, so a missing value follows the branch direction learned for it.
    """
    return booster.predict(xgb.DMatrix(passthrough_nan_check(X),
                                       feature_names=feature_names, missing=np.nan))


def trees_contributions(booster: xgb.Booster, X: np.ndarray,
                        feature_names: list[str]) -> np.ndarray:
    """Per-frame additive contributions: exact TreeSHAP (pred_contribs), no bias column.

    Returns [n_rows, p] in the margin (logit) units of the booster, in the order
    of feature_names. xgboost appends the bias term as a last column; it is
    dropped, so the columns sum to the margin minus the bias and are comparable
    with the linear arm's coef·z. TreeSHAP is exact for a tree ensemble (Lundberg
    et al. (2020) Nat. Mach. Intell. 2:56-67, doi:10.1038/s42256-019-0138-9), so
    the mechanism read-off does not rest on a sampling approximation.
    """
    c = booster.predict(xgb.DMatrix(passthrough_nan_check(X),
                                    feature_names=feature_names, missing=np.nan),
                        pred_contribs=True)
    return c[:, :-1]
