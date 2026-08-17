"""Global importances and per-frame mechanistic diagnosis.

The diagnosis is built from the contributions of the ALREADY TRAINED model
(additive, in logits) rather than from a second scoring system alongside it:
the contributions share one unit, and their sum plus the base level yields the
prediction. Importances are aggregated to blocks, since within-block
redistribution among near-duplicate columns is arbitrary.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from .config import Config
from .features._common import ColSpec

# The mechanism dictionary is derived from the domain literature:
# each failure route is its own axis; each axis is backed by column blocks.
AXES: dict[str, list[str]] = {
    "constriction dehydration": ["hydration", "fluctuations"],
    "pathway/lining rearrangement": ["geometry", "symmetry"],
    "ion placement and desolvation": ["occupancy", "named_sites"],
    "water-wire rupture": ["water_wire"],
    "electrostatics": ["electrostatics"],
    "column immobility": ["dynamics"],
    "permeant delivery": ["delivery"],
}


def interpret_run(cfg: Config, oof: pd.DataFrame, contrib: pd.DataFrame,
                  schema: list[ColSpec], model_columns: list[str],
                  merged: pd.DataFrame, out, log) -> None:
    """Write importances.json, mechanism_specificity.json and answers.parquet.

    contrib holds one row per frame and one column per model column, the
    additive contribution of that column to that frame's score in logits, plus
    the three key columns (condition, replica, time_ps). oof supplies the label,
    the validity flag, the fold and the calibrated probability; merged supplies
    the measured descriptor values, which are quoted in the reason string. out
    is the train step's output directory; when folds.json is already written
    there, each frame's verdict is taken at the threshold of the fold that
    predicted that frame.

    Nothing is fitted here. The diagnosis is a reading of the model that already
    exists, so it cannot disagree with the probability it explains.
    """
    primary = cfg["labels.primary_tau_ps"]
    head = cfg["model.head"]
    if contrib.empty:
        (out / "importances.json").write_text(json.dumps({"empty": True}))
        return
    block_of = {c.name: c.block for c in schema}
    feat_cols = [c for c in contrib.columns if c not in ("condition", "replica", "time_ps")]
    C = contrib[feat_cols].to_numpy(float)

    # Global importances: mean |contribution| per frame; aggregated to blocks.
    col_imp = {c: float(np.nanmean(np.abs(C[:, j]))) for j, c in enumerate(feat_cols)}
    block_imp: dict[str, float] = {}
    for c, v in col_imp.items():
        block_imp[block_of[c]] = block_imp.get(block_of[c], 0.0) + v
    (out / "importances.json").write_text(json.dumps(
        {"per_column": col_imp, "per_block": block_imp,
         "note": "attribution kind is recorded in attribution.json; aggregated to blocks"},
        indent=1, ensure_ascii=False))

    # Axis availability: the diagnosis dictionary and the model matrix are
    # different sets.
    axis_cols = {ax: [c for c in feat_cols if block_of[c] in blocks]
                 for ax, blocks in AXES.items()}
    availability = {ax: {"available": bool(cols), "n_columns": len(cols)}
                    for ax, cols in axis_cols.items()}
    unavailable = [ax for ax, a in availability.items() if not a["available"]]
    if unavailable:
        log.say(f"axes unavailable for this protein (no active columns): {unavailable}")

    # Per-frame axis contribution toward "not ready": the sum of NEGATIVE
    # contributions of its columns, taken in absolute value. A missing
    # measurement yields zero contribution, a property of the contributions
    # themselves.
    axis_notready = {}
    for ax, cols in axis_cols.items():
        if not cols:
            continue
        sub = C[:, [feat_cols.index(c) for c in cols]]
        neg = np.where(np.isnan(sub), 0.0, np.minimum(sub, 0.0))
        axis_notready[ax] = -neg.sum(axis=1)
    ax_names = list(axis_notready)
    A = np.vstack([axis_notready[a] for a in ax_names]).T if ax_names else np.empty((len(C), 0))

    head_oof = oof[(oof["arm"] == head) & (oof["tau"] == primary)].copy()
    key = ["condition", "replica", "time_ps"]
    dd = contrib[key].merge(head_oof, on=key, how="left", validate="one_to_one")
    # The scale against which "large" is judged is the base level of the system,
    # in the
    # same logit units as the contributions: log(p0/(1-p0)) at the base rate. A
    # rare-event
    # system starts far below even odds, so a contribution that would be
    # decisive there is
    # noise on a system whose base rate is one half. Comparing against a fixed
    # number of
    # logits instead would make the axis verdict a function of the base rate.
    base_rate = float(dd["y"].mean()) if len(dd) else 0.5
    base_logit = abs(float(np.log(max(base_rate, 1e-6) / max(1 - base_rate, 1e-6))))
    negligible = cfg["interpret.negligible_frac"] * max(base_logit, 1e-6)

    # Axis specificity — computed BEFORE picking the winner, as an odds ratio.
    # An axis
    # that fires on almost every frame will still win the argmax on some of them
    # and read
    # as an explanation, so the discriminating power of each axis is published
    # first and
    # independently: the odds of the axis firing on not-ready frames over the odds on
    # ready frames. Only the valid frames enter, the same set the metrics use.
    fires = A > negligible
    specificity = {}
    y = dd["y"].to_numpy(float)
    valid = dd["valid"].to_numpy(bool)
    for j, ax in enumerate(ax_names):
        f = fires[valid, j]
        yy = y[valid]
        p1 = float(np.mean(f[yy == 1])) if np.any(yy == 1) else np.nan
        p0 = float(np.mean(f[yy == 0])) if np.any(yy == 0) else np.nan
        odds = lambda p: p / (1 - p) if np.isfinite(p) and 0 < p < 1 else np.nan
        specificity[ax] = {"fires_on_ready_frac": p1, "fires_on_notready_frac": p0,
                           "odds_ratio_notready_vs_ready":
                               odds(p0) / odds(p1) if np.isfinite(odds(p0)) and np.isfinite(odds(p1)) else np.nan,
                           "fires_overall_frac": float(np.mean(f))}
    # An axis firing on more than 95% of frames is a constant, whatever its odds
    # ratio:
    # it names itself on ready and not-ready frames alike and carries no
    # diagnosis. The
    # warning is published in the artifact and printed, not silently applied —
    # the axis
    # stays in the vocabulary and the reader decides what to do with it.
    near_universal = [ax for ax, s in specificity.items()
                      if np.isfinite(s["fires_overall_frac"]) and s["fires_overall_frac"] > 0.95]
    (out / "mechanism_specificity.json").write_text(json.dumps(
        {"axes": specificity, "availability": availability,
         "negligible_threshold_logit": negligible,
         "near_universal_axes_warning": near_universal},
        indent=1, ensure_ascii=False))
    if near_universal:
        log.say(f"WARNING: axes named on nearly every frame (they discriminate nothing): "
                f"{near_universal}")

    # Final per-frame answer: probability + READY/NOT_READY + axis + reason.
    # Each fold has its own threshold, chosen inside that fold, and a frame is
    # thresholded by the fold that predicted it. One threshold for the whole
    # table would
    # be selected partly on frames the model was trained on.
    thr_by_fold = {}
    import json as _json
    folds_doc = _json.loads((out / "folds.json").read_text()) if (out / "folds.json").exists() else []
    for fa in folds_doc:
        if fa["tau"] == primary and head in fa.get("arms", {}) and "threshold" in fa["arms"][head]:
            thr_by_fold[fa["fold"]] = fa["arms"][head]["threshold"]
    answers = dd[key].copy()
    answers["tau_ps"] = primary
    answers["prob_ready"] = dd["prob_cal"]
    thr = dd["fold"].map(thr_by_fold).astype(float)
    answers["ready"] = (dd["prob_cal"] >= thr).map({True: "READY", False: "NOT_READY"})
    winner_idx = np.argmax(A, axis=1) if A.shape[1] else np.zeros(len(dd), int)
    winner_val = A[np.arange(len(A)), winner_idx] if A.shape[1] else np.zeros(len(dd))
    axis_lbl = np.array([ax_names[k] for k in winner_idx]) if ax_names else np.array([])
    axis_lbl = np.where(winner_val > negligible, axis_lbl, "unclear")
    answers["mechanism_axis"] = axis_lbl
    reasons = []
    val_lookup = merged.set_index(key)
    for i in range(len(dd)):
        if axis_lbl[i] == "unclear":
            reasons.append("largest contribution is negligible against the base level")
            continue
        ax = axis_lbl[i]
        cols = axis_cols[ax]
        sub = C[i, [feat_cols.index(c) for c in cols]]
        sub = np.where(np.isnan(sub), 0.0, sub)
        # The most negative contribution, i.e. the single column pushing hardest
        # toward
        # not-ready, names the reason; its measured value is quoted next to it so the
        # statement can be checked against the descriptor table.
        top_j = int(np.argmin(sub))
        top_col = cols[top_j]
        v = val_lookup.loc[tuple(dd.loc[i, key]), top_col] if top_col in merged.columns else np.nan
        reasons.append(f"{top_col}={v:.3g} (contribution {sub[top_j]:.3f} logit toward not-ready)"
                       if np.isfinite(sub[top_j]) else "contribution undefined")
    answers["reason"] = reasons
    answers.to_parquet(out / "answers.parquet")
    n_unclear = int(np.sum(axis_lbl == "unclear"))
    log.say(f"per-frame answers: {len(answers)}; 'unclear' on {n_unclear} frames; "
            f"block importances: " + ", ".join(f"{b}={v:.2f}" for b, v in
                                             sorted(block_imp.items(), key=lambda kv: -kv[1])[:5]))
