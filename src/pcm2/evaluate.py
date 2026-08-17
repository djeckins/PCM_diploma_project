"""Metrics, paired comparison, and the monitor over time.

The headline per-frame metric is average precision (area under the PR curve)
paired with its ratio to the base rate; the headline monitor quantity is
event-level: the fraction of warned crossings at the observed empty-alarm
rate, always shown next to a random null comparison. ECE is not computed: its
informative content is carried by the calibration slope and intercept.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from .config import Config
from .events import load_events
from .runtime import derive_seed


def _logit(p):
    """Log-odds of a probability, clipped so that a saturated value stays finite.

    The clip at 1e-7 bounds the axis at about ±16.1. Without it one frame whose
    calibrated probability saturated at exactly 0 or 1 would put an infinity into
    the regressor of the calibration line and the fit would fail.
    """
    p = np.clip(p, 1e-7, 1 - 1e-7)
    return np.log(p / (1 - p))


def ap_with_base(y: np.ndarray, s: np.ndarray) -> dict:
    """Average precision of the ranking s for the labels y, with its base rate.

    y is 0/1, s is any monotone score (higher = readier); they need not be
    probabilities, which is what lets the published-criterion arm be scored on the
    same footing. A random ranker reaches an average precision equal to the base
    rate, so ratio_to_chance = ap / base_rate is the quantity that survives a
    change of prevalence and is the only one compared ACROSS horizons: a longer
    horizon raises the base rate mechanically, and the bare AP with it.

    Frames whose score is not finite are dropped from both arrays — an arm that
    cannot score a frame abstains there; it is not credited with a value.
    """
    m = np.isfinite(s)
    y, s = y[m], s[m]
    # One class only: the ranking cannot be scored at all. Returning zero would
    # enter the record as a bad model instead of an unusable sample, so the caller
    # gets a refusal with the reason and records the fold by name.
    if len(y) == 0 or y.sum() == 0 or y.sum() == len(y):
        return {"defined": False, "reason": "both classes absent — the metric is "
                                            "undefined, not zero"}
    ap = float(average_precision_score(y, s))
    base = float(np.mean(y))
    return {"defined": True, "ap": ap, "base_rate": base, "ratio_to_chance": ap / base,
            "n": int(len(y)), "n_pos": int(y.sum())}


def calibration_line(y: np.ndarray, p: np.ndarray) -> dict:
    """Slope and intercept of the calibration line (logit-logit); reported values.

    The label is regressed on the log-odds of the predicted probability with an
    unpenalized logistic fit: slope 1 and intercept 0 is a perfectly calibrated
    forecast, a slope below 1 means the probabilities are spread too far towards
    the extremes, and the intercept carries the constant bias in the prior. These
    two numbers are reported in place of a binned ECE, which needs an arbitrary
    binning and hides the direction of the error.

    A negative slope means the probability is ordered against the outcome — the
    forecast is worse than uninformative and the flag says so instead of leaving a
    reader to notice a minus sign.
    """
    m = np.isfinite(p)
    if y[m].sum() == 0 or y[m].sum() == m.sum():
        return {"defined": False}
    lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
    lr.fit(_logit(p[m])[:, None], y[m])
    slope = float(lr.coef_[0][0])
    return {"defined": True, "slope": slope, "intercept": float(lr.intercept_[0]),
            "negative_slope_alarm": slope < 0}


def paired_bootstrap(oof: pd.DataFrame, arm_a: str, arm_b: str, tau: float,
                     n_boot: int, alpha: float, seed: int) -> dict:
    """Paired bootstrap OVER TRAJECTORIES: one resample scores both methods.

    Returns the mean AP difference arm_a − arm_b and its two-sided percentile
    interval at level alpha. The resampling unit is the whole trajectory: frames
    inside one trajectory are strongly dependent, and resampling frames would
    treat every near-duplicate as an independent observation and shrink the
    interval to nothing. Both arms are scored on the SAME resample, so the part of
    the spread that comes from which trajectories were drawn cancels in the
    difference; this is the comparison the paired design buys.
    """
    sub = oof[(oof["tau"] == tau) & oof["valid"] & oof["arm"].isin([arm_a, arm_b])]
    # One row per frame with a column per arm: the pairing is built here. Frames one
    # of the two arms could not score are dropped, so both arms are compared on
    # exactly the same set of frames.
    piv = sub.pivot_table(index=["condition", "replica", "time_ps", "y"],
                          columns="arm", values="score_raw").reset_index()
    piv = piv.dropna(subset=[arm_a, arm_b])
    trajs = (piv["condition"] + "/" + piv["replica"]).to_numpy()
    uniq = np.unique(trajs)
    # The stream is derived from the one recorded seed together with the arm names
    # and the horizon, so every comparison is reproducible and no two share a stream.
    rng = np.random.default_rng(derive_seed(seed, "boot", arm_a, arm_b, str(tau)))
    diffs = []
    for _ in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.flatnonzero(trajs == t) for t in pick])
        y = piv["y"].to_numpy()[idx]
        # With few trajectories and a rare positive class a resample can come out
        # single-class; average precision is undefined there, so the resample is
        # skipped and n_boot_effective records how many actually contributed.
        if y.sum() == 0 or y.sum() == len(y):
            continue
        d = (average_precision_score(y, piv[arm_a].to_numpy()[idx])
             - average_precision_score(y, piv[arm_b].to_numpy()[idx]))
        diffs.append(d)
    if not diffs:
        return {"defined": False}
    lo, hi = np.percentile(diffs, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return {"defined": True, "mean_diff": float(np.mean(diffs)),
            "ci": [float(lo), float(hi)], "n_groups": int(len(uniq)),
            "n_boot_effective": len(diffs),
            "method": "percentile paired bootstrap over trajectories",
            "coverage_note": "with few groups the interval coverage is not "
                             "guaranteed; the interval is shown, not hidden"}


def hysteresis_alarm(prob: np.ndarray, on_thr: float, off_thr: float) -> np.ndarray:
    """Armed/not-armed state of the monitor, frame by frame, from a probability series.

    The alarm arms at on_thr and releases only once the probability falls below
    off_thr, with off_thr < on_thr. A single threshold on a fluctuating probability
    chatters: one approach to conduction is cut into many short episodes, which both
    inflates the episode count and makes the empty-alarm rate unreadable. The dead
    band between the two thresholds makes an episode a sustained state rather than
    an artifact of the noise.
    """
    state = np.zeros(len(prob), dtype=bool)
    armed = False
    for i, p in enumerate(prob):
        # A frame the model could not score is no evidence either way: the state is
        # carried over rather than treated as a release.
        if not np.isfinite(p):
            state[i] = armed
            continue
        if not armed and p >= on_thr:
            armed = True
        elif armed and p < off_thr:
            armed = False
        state[i] = armed
    return state


def smooth_backward(x: np.ndarray, times: np.ndarray, window_ps: float) -> np.ndarray:
    """Strictly backward smoothing: future frames are unavailable online.

    Mean of x over the trailing window of window_ps, one value per frame, in the
    units of x. The window is a time and not a frame count, so the smoothing means
    the same thing at any stride. times must be ascending — the left edge j only
    ever moves forward. A frame whose window holds no finite value stays NaN, which
    the alarm reads as "no evidence" rather than as a low probability.
    """
    out = np.full(len(x), np.nan)
    j = 0
    for i in range(len(x)):
        while times[i] - times[j] > window_ps:
            j += 1
        seg = x[j:i + 1]
        seg = seg[np.isfinite(seg)]
        if len(seg):
            out[i] = float(np.mean(seg))
    return out


def alarm_metrics(times: np.ndarray, armed: np.ndarray,
                  exit_times: np.ndarray, tau_ps: float) -> dict:
    """Event-level score of one trajectory's alarm series: warned crossings and lead times.

    times and exit_times in ps on the same clock, armed the state series of
    hysteresis_alarm over the same frames. The headline number is warned_frac, the
    share of the trajectory's crossings for which the alarm was already up inside
    the horizon; it means nothing without duty_cycle beside it, since an alarm that
    stays armed warns every crossing.
    """
    # Rising edges of the state series delimit the episodes. The prepended False
    # makes a trajectory that starts already armed count its first frame as an edge.
    starts = np.flatnonzero(armed & ~np.concatenate([[False], armed[:-1]]))
    n_episodes = int(len(starts))
    warned, leads = 0, []
    for tx in exit_times:
        # The same window that defines the label: a crossing is forewarned if the
        # alarm was up on at least one frame of [tx-τ, tx). Using the label's own
        # window is what keeps the event-level number comparable with the per-frame
        # metrics instead of measuring a differently posed question.
        before = np.flatnonzero((times < tx) & (times >= tx - tau_ps))
        if len(before) and armed[before].any():
            warned += 1
            first = before[armed[before]][0]
            # Lead time is measured from the start of the episode that gave the
            # warning, walking back through the contiguous armed frames — this may
            # reach further back than τ. The quantity of interest is how long the
            # warning had been standing when the ion left, not the position of the
            # first armed frame inside the horizon.
            ep_start = first
            while ep_start > 0 and armed[ep_start - 1]:
                ep_start -= 1
            leads.append(float(tx - times[ep_start]))
    empty = 0
    for s in starts:
        # Walk to the end of this episode: e-1 is its last armed frame.
        e = s
        while e < len(armed) and armed[e]:
            e += 1
        span = (times[s], times[e - 1] if e > 0 else times[s])
        # an episode is empty if no crossing completed during its span + τ
        # The +τ is not slack: an alarm still up on the last frame of the episode is
        # a statement about the following τ, and would be judged unfairly if the
        # crossing it predicted fell just past the release.
        if not np.any((exit_times >= span[0]) & (exit_times <= span[1] + tau_ps)):
            empty += 1
    # A trajectory with no crossings gets warned_frac = NaN, not 0: there was
    # nothing to warn about, which is not the same as failing to warn.
    return {"n_crossings": int(len(exit_times)), "n_warned": warned,
            "warned_frac": warned / len(exit_times) if len(exit_times) else np.nan,
            "n_episodes": n_episodes, "n_empty_episodes": empty,
            "duty_cycle": float(np.mean(armed)),
            "lead_ps": {"median": float(np.median(leads)) if leads else np.nan,
                        "p10": float(np.percentile(leads, 10)) if leads else np.nan,
                        "p90": float(np.percentile(leads, 90)) if leads else np.nan}}


def _random_null(times, armed, exit_times, tau_ps, seed, n_shifts=50) -> dict:
    """Null comparison: the same alarm episodes, randomly shifted cyclically.
    Gives the chance level for the warned-crossing fraction.

    A cyclic shift keeps the number of episodes, their durations and the duty cycle
    and destroys only their alignment with the crossings. The mean warned fraction
    over the shifts is therefore the level an alarm of this shape reaches without
    knowing anything about the trajectory, which is what the reported fraction has
    to be read against: with a high duty cycle a large warned fraction is free."""
    rng = np.random.default_rng(seed)
    fracs, empties = [], []
    for _ in range(n_shifts):
        # Shift of at least one frame: a zero shift would return the observed
        # series and pull the null towards the measurement it is meant to test.
        k = int(rng.integers(1, len(armed)))
        sh = np.roll(armed, k)
        m = alarm_metrics(times, sh, exit_times, tau_ps)
        fracs.append(m["warned_frac"])
        empties.append(m["n_empty_episodes"])
    return {"warned_frac_mean": float(np.nanmean(fracs)),
            "method": "cyclic shift of the alarm series, preserving episode structure"}


def evaluate_run(cfg: Config, oof: pd.DataFrame, merged: pd.DataFrame,
                 fold_arts: list, results: list, out, log) -> None:
    """Score the out-of-fold predictions of a run and write evaluation.json and alarm.json.

    oof carries one row per (frame, horizon, arm) with the raw score, the calibrated
    probability, the label, the validity flag and the transit flag; merged is the
    feature table joined to the labels, used where a metric needs the labels of a
    fold's test rows directly. No predictor is fitted here: the only thing estimated
    on held-out scores is the two-parameter calibration line, which is reported as a
    diagnostic and never applied to a score, and the alarm thresholds, which for
    each trajectory come from the other trajectories.

    The blocks, in order: per-arm metrics at every horizon, the paired comparison
    against the declared head arm at the primary horizon, the resting/transit anchor
    ablation, the per-block ablation, transfer between conditions, and the monitor.
    """
    taus = cfg["labels.tau_ps"]
    primary = cfg["labels.primary_tau_ps"]
    head = cfg["model.head"]
    arms = sorted(oof["arm"].unique())
    evaluation: dict = {"per_arm": {}, "degenerate_folds": [],
                        "headline_note": "headline per-frame metric is pooled AP/base "
                                         "at the primary τ; the headline monitor "
                                         "quantity is event-level"}
    for arm in arms:
        evaluation["per_arm"][arm] = {}
        for tau in taus:
            # Only valid frames enter any metric: the censored tail carries no known
            # label, and scoring it as negative would flatter every arm equally.
            sub = oof[(oof["arm"] == arm) & (oof["tau"] == tau) & oof["valid"]]
            pooled = ap_with_base(sub["y"].to_numpy(), sub["score_raw"].to_numpy())
            per_fold = {}
            for f in sorted(sub["fold"].unique()):
                ff = sub[sub["fold"] == f]
                r = ap_with_base(ff["y"].to_numpy(), ff["score_raw"].to_numpy())
                # A fold whose held-out assembly carries no crossing cannot be
                # scored. It is listed by name so that the pooled number is read
                # knowing which assemblies contributed to it, instead of being
                # dropped quietly or entered as a zero.
                if not r["defined"]:
                    evaluation["degenerate_folds"].append(
                        {"arm": arm, "tau": tau, "fold": f, "reason": r["reason"]})
                per_fold[f] = r
            entry = {"pooled": pooled, "per_fold": per_fold}
            m = np.isfinite(sub["score_raw"].to_numpy())
            yy, ss = sub["y"].to_numpy()[m], sub["score_raw"].to_numpy()[m]
            # AUROC is recorded and deliberately not compared. Where the positive
            # class is rare — base rates of a tenth of a per cent to a few per cent
            # in the KcsA arms — it is governed by the ordering of the large negative
            # class and barely moves when the top of the ranking is wrong. Average
            # precision against the base rate is what the study compares.
            if pooled.get("defined"):
                entry["auroc_report_only"] = float(roc_auc_score(yy, ss))
            # Probabilistic metrics only for arms that expose a probability scale.
            # The published-criterion arms score a frame with a negated energy in
            # kJ/mol; it ranks frames but is not a probability, and a Brier score or
            # a calibration line computed on it would be a number without a meaning.
            if arm not in ("published_rao2019", "published_rao2019_win"):
                p = sub["prob_cal"].to_numpy()
                if np.any(np.isfinite(p)) and pooled.get("defined"):
                    entry["brier"] = float(brier_score_loss(
                        sub["y"].to_numpy()[np.isfinite(p)], p[np.isfinite(p)]))
                    entry["calibration_pooled"] = calibration_line(sub["y"].to_numpy(), p)
                    entry["calibration_per_fold"] = {
                        f: calibration_line(sub[sub["fold"] == f]["y"].to_numpy(),
                                            sub[sub["fold"] == f]["prob_cal"].to_numpy())
                        for f in sorted(sub["fold"].unique())}
            evaluation["per_arm"][arm][str(tau)] = entry
    # Paired comparison at the primary τ: the head arm against every other arm.
    # When the head refused every fold (event-starved system: all folds
    # degenerate) there is nothing to pair against; the refusal is recorded
    # instead of crashing on a missing pivot column.
    evaluation["paired_vs_head"] = {}
    head_scores = oof[(oof["arm"] == head) & (oof["tau"] == primary)]["score_raw"]
    if len(head_scores) and np.isfinite(head_scores.to_numpy(float)).any():
        for arm in arms:
            if arm == head:
                continue
            evaluation["paired_vs_head"][arm] = paired_bootstrap(
                oof, head, arm, primary, cfg["evaluate.n_boot"], cfg["evaluate.ci_alpha"],
                cfg["model.seed"])
    else:
        evaluation["paired_vs_head_note"] = (
            f"head arm '{head}' has no out-of-fold scores at τ={primary}: "
            "every fold refused (degenerate training); paired comparison "
            "is undefined by construction")
    # Anchor ablation: transit frames versus resting frames.
    # On a transit frame the pore already holds the ion whose exit defines the
    # label, so occupancy and hydration descriptors partly describe the crossing
    # itself; on a resting frame they cannot. The two are scored apart so that the
    # claim of the study rests on the half where the question is not circular.
    sub = oof[(oof["arm"] == head) & (oof["tau"] == primary) & oof["valid"]]
    evaluation["anchor_ablation"] = {
        "resting": ap_with_base(sub[~sub["in_transit"]]["y"].to_numpy(),
                                sub[~sub["in_transit"]]["score_raw"].to_numpy()),
        "transit": ap_with_base(sub[sub["in_transit"]]["y"].to_numpy(),
                                sub[sub["in_transit"]]["score_raw"].to_numpy()),
        "note": "the study's conclusion must hold on resting frames"}
    # Block ablation: head-arm AP with each block removed, per fold.
    # The training step already refitted the arm without each descriptor block using
    # that fold's chosen hyperparameters and no new search, so what changes between
    # the full arm and an ablated one is the block and not the tuning. Reported per
    # fold rather than pooled: with a handful of assemblies the spread across folds
    # is the honest measure of how firm a drop is.
    ablation: dict = {}
    for res in results:
        fa = res["fold_art"]
        if fa["tau"] != primary or "block_ablation" not in fa:
            continue
        # The ablated scores were stored in the order of the fold's test rows, so the
        # labels are taken from the same rows of the merged table.
        te = np.asarray(res["te_idx"])
        y = merged[f"y_tau{int(primary)}"].to_numpy(int)[te]
        v = merged[f"valid_tau{int(primary)}"].to_numpy(bool)[te]
        for block, d in fa["block_ablation"].items():
            s = np.asarray(d["scores"])
            r = ap_with_base(y[v], s[v])
            ablation.setdefault(block, {})[fa["fold"]] = r
    evaluation["block_ablation"] = ablation
    # Transfer between conditions.
    # The training step fits the head arm on one condition and predicts another; a
    # score that collapses here says the arm has learned the condition label rather
    # than a property of the frame. Only the ranking is scored — nothing is refitted
    # on the target condition.
    ct_path = out / "condition_transfer.json"
    if ct_path.exists():
        ct = json.loads(ct_path.read_text())
        ct_m = {}
        for k, d in ct.items():
            if "scores" in d:
                y = np.asarray(d["y"])
                v = np.asarray(d["valid"])
                s = np.asarray(d["scores"])
                ct_m[k] = ap_with_base(y[v], s[v])
        evaluation["condition_transfer"] = ct_m

    # Monitor: hysteresis, lead time, random null.
    # This is the operating form of the model: a probability read frame by frame,
    # smoothed backward, turned into an alarm with a dead band, and scored per
    # crossing rather than per frame.
    ev = load_events(cfg)
    # Half the primary horizon if no smoothing time is declared: long enough to
    # damp the frame-to-frame noise, short enough that a warning issued inside the
    # horizon still arrives before the crossing it announces.
    smooth_ps = cfg["evaluate.alarm.smooth_ps"] or primary / 2
    alarm_out: dict = {"smooth_ps": smooth_ps,
                       "thresholds_rule": "quantiles of the calibrated OOF scores of the "
                                          "remaining folds (out of this fold's sample)"}
    head_oof = oof[(oof["arm"] == head) & (oof["tau"] == primary)]
    per_traj, nulls = {}, {}
    for (c, r), g in head_oof.groupby(["condition", "replica"]):
        # Frames must be in time order before anything backward-looking is applied.
        g = g.sort_values("time_ps")
        # Both thresholds come from the OTHER trajectories: a quantile of this
        # trajectory's own probabilities would be a threshold tuned on the sample it
        # is about to be judged on, and it would also force an alarm to fire in a
        # trajectory that never conducts. on_quantile > off_quantile, so the arming
        # level sits above the release level and the dead band is non-empty.
        others = head_oof[~((head_oof["condition"] == c) & (head_oof["replica"] == r))]
        pool = others["prob_cal"].to_numpy()
        pool = pool[np.isfinite(pool)]
        if not len(pool):
            continue
        on_thr = float(np.quantile(pool, cfg["evaluate.alarm.on_quantile"]))
        off_thr = float(np.quantile(pool, cfg["evaluate.alarm.off_quantile"]))
        times = g["time_ps"].to_numpy()
        sm = smooth_backward(g["prob_cal"].to_numpy(), times, smooth_ps)
        armed = hysteresis_alarm(sm, on_thr, off_thr)
        cond = next(cc for cc in cfg["data.conditions"] if cc["id"] == c)
        sel = ev[(ev["condition"] == c) & (ev["replica"] == r)]
        # The same direction filter the labels step applies: the alarm is scored
        # against the crossings it was trained to anticipate.
        if cond.get("direction") is not None:
            sel = sel[(sel["direction"] == cond["direction"]) | (sel["direction"] == 0)]
        exits = sel["t_exit_ps"].to_numpy()
        met = alarm_metrics(times, armed, exits, primary)
        met["on_thr"], met["off_thr"] = on_thr, off_thr
        per_traj[f"{c}/{r}"] = met
        nulls[f"{c}/{r}"] = _random_null(times, armed, exits, primary,
                                         derive_seed(cfg["model.seed"], "null", c, r))
    alarm_out["per_trajectory"] = per_traj
    alarm_out["random_null"] = nulls
    if per_traj:
        # Mean over trajectories, not over crossings: a single long trajectory with
        # many crossings would otherwise decide the headline number on its own.
        # Trajectories without crossings contribute a NaN warned fraction and are
        # left out of the mean rather than counted as failures.
        wf = [m["warned_frac"] for m in per_traj.values() if np.isfinite(m["warned_frac"])]
        nf = [n["warned_frac_mean"] for n in nulls.values() if np.isfinite(n["warned_frac_mean"])]
        alarm_out["headline_event_level"] = {
            "warned_frac_mean": float(np.mean(wf)) if wf else np.nan,
            "random_null_warned_frac_mean": float(np.mean(nf)) if nf else np.nan,
            "empty_episodes_total": int(sum(m["n_empty_episodes"] for m in per_traj.values())),
            "episodes_total": int(sum(m["n_episodes"] for m in per_traj.values()))}
    (out / "evaluation.json").write_text(json.dumps(evaluation, indent=1,
                                                    ensure_ascii=False, default=str))
    (out / "alarm.json").write_text(json.dumps(alarm_out, indent=1,
                                               ensure_ascii=False, default=str))
    hp = evaluation["per_arm"].get(head, {}).get(str(primary), {}).get("pooled", {})
    if hp.get("defined"):
        log.say(f"head arm {head} at τ={primary}: pooled AP={hp['ap']:.3f} "
                f"at base rate {hp['base_rate']:.3f} "
                f"(×{hp['ratio_to_chance']:.2f} over chance)")
    if alarm_out.get("headline_event_level"):
        h = alarm_out["headline_event_level"]
        log.say(f"monitor: warned fraction {h['warned_frac_mean']:.2f} of crossings "
                f"(random null {h['random_null_warned_frac_mean']:.2f}); "
                f"empty episodes {h['empty_episodes_total']}/{h['episodes_total']}")
