"""Figure catalog: one figure answers one question; this list is the source of truth.

Drawing rules hard-coded here: PR is a step function with markers on the
attainable points (no marker on the library's padded endpoint); the horizon
axis shows only the ratio to chance, with a unity line; calibration axes span
the observed range, the diagonal is drawn last, the zoom is stated in the
title; tick labels are horizontal and sorted by median; the base-rate line is
mandatory. A figure outside the catalog is never saved: the attempt must fail.
"""

from __future__ import annotations

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve

from .config import Config
from .runtime import StepLog, run_dir, step_output


class CatalogError(Exception):
    """Raised when a figure is written under a name the catalog does not declare."""


def _save(fig, out_dir, name: str) -> None:
    """Write one catalog figure as out_dir/<name>.png and release the canvas.

    The membership check is where the catalog rule is enforced. The figures step
    replaces its whole directory on each run, so a plot saved under an
    undeclared name would exist until the next run and then vanish without a
    trace in catalog.json; failing here turns that into a visible error.
    """
    if name not in FIGURES:
        raise CatalogError(f"figure {name} is not in the catalog — the next run would delete it")
    fig.savefig(out_dir / f"{name}.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


def _load(cfg: Config):
    """Read the four train-step artifacts every figure draws from.

    Returns (oof, evaluation, importances, attribution): the per-frame
    out-of-fold scores and labels, the metric document, the per-column and
    per-block importances, and the record of which attribution method produced
    them. Nothing here is recomputed — the figures are a view of the train
    step, so a figure and the report always quote the same numbers.
    """
    troot = run_dir(cfg) / "train"
    oof = pd.read_parquet(troot / "oof.parquet")
    evaluation = json.loads((troot / "evaluation.json").read_text())
    imp = json.loads((troot / "importances.json").read_text())
    attribution = json.loads((troot / "attribution.json").read_text())
    return oof, evaluation, imp, attribution


# The thesis comparison set: exactly these five methods plus the base-rate
# (chance) line appear in every cross-method figure, in this order. The other
# benchmark arms (fitted controls, the Rao energy-heuristic variants) stay in
# the benchmark tables but out of the figures.
COMPARE_ARMS = [
    # (arm, label, fixed color) — the color belongs to the METHOD, not to the
    # plot order, so a system where some arm is absent keeps every other
    # method's color identical to the rest of the thesis.
    ("published_rao2019_sigma", "Rao 2019 Σd (published structure rule)", "#1f77b4"),
    ("plsfma_coords", "PLS-FMA (g_fma port)", "#ff7f0e"),
    ("linear", "Linear (ridge-logistic)", "#2ca02c"),
    ("trees", "GBT (this work)", "#d62728"),
    ("clock", "Clock (temporal null)", "#9467bd"),
]


def _sigma_frame_scores(cfg):
    """Per-frame ranking score of the published Rao structure rule (−Σd),
    or None for feature tables built before the classifier was added."""
    import pyarrow.parquet as pq
    p = run_dir(cfg) / "features" / "features.parquet"
    if "bl_rao_sigma_d" not in pq.read_schema(p).names:
        return None
    f = pd.read_parquet(p, columns=["condition", "replica", "time_ps",
                                    "bl_rao_sigma_d"])
    f["score_raw"] = -f["bl_rao_sigma_d"]
    return f.drop(columns=["bl_rao_sigma_d"])


def _arm_frames(oof, sigma, arm, tau):
    """Valid frames with a score for one comparison arm at one horizon.

    Fitted arms come straight from the OOF table; the Rao Σd rule is not
    fitted, so its score is merged onto the head arm's frame set (same frames,
    same labels — only the score column differs)."""
    if arm == "published_rao2019_sigma":
        if sigma is None:
            return None
        # Frames and labels are shared across arms; the carrier must be the
        # arm with the WIDEST coverage (the head may cover only some folds on
        # event-starved systems, and a narrow carrier would silently shrink
        # the published rule's evaluation to those frames).
        cov = {a: int(((oof["arm"] == a) & (oof["tau"] == tau)).sum())
               for a in ("trees", "linear", "clock")}
        carrier = max(cov, key=cov.get)
        if cov[carrier] == 0:
            return None
        sub = oof[(oof["arm"] == carrier) & (oof["tau"] == tau) & oof["valid"]]
        sub = sub.merge(sigma, on=["condition", "replica", "time_ps"],
                        suffixes=("_fit", ""), validate="one_to_one")
    else:
        sub = oof[(oof["arm"] == arm) & (oof["tau"] == tau) & oof["valid"]]
    return sub


def fig_pr_primary(cfg, out_dir, oof, evaluation, imp, attribution):
    """Ranking quality: step-function PR for the six-way comparison, five
    methods plus the mandatory base-rate line."""
    primary = cfg["labels.primary_tau_ps"]
    sigma = _sigma_frame_scores(cfg)
    fig, ax = plt.subplots(figsize=(7, 5))
    base = None
    for arm, label, color in COMPARE_ARMS:
        sub = _arm_frames(oof, sigma, arm, primary)
        if sub is None:
            continue
        m = np.isfinite(sub["score_raw"].to_numpy())
        y, s = sub["y"].to_numpy()[m], sub["score_raw"].to_numpy()[m]
        # Precision and recall are undefined with one class only: with no
        # positive there is nothing to recall, with no negative precision is 1
        # everywhere. Such an arm is left out of the panel rather than drawn flat.
        if len(y) == 0 or y.sum() in (0, len(y)):
            continue
        prec, rec, thr = precision_recall_curve(y, s)
        base = float(np.mean(y))
        # The library appends a padded endpoint — draw the steps without it.
        ax.step(rec[:-1], prec[:-1], where="post", label=label, color=color)
        # Few distinct thresholds means few attainable operating points, and the
        # steps between them are interpolation rather than measurement; the
        # markers say where the measured points are. Above ~60 they merge into a
        # band and are dropped.
        if len(thr) < 60:
            ax.plot(rec[:-1], prec[:-1], "o", ms=3, color=color)
    if base is None:
        # No arm could rank a single positive on its covered frames; the
        # catalog records the figure as not produced.
        plt.close(fig)
        return
    ax.axhline(base, ls="--", c="gray", label="base rate (chance level)")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.set_title(f"Per-frame PR, τ={int(primary)} ps (pooled OOF)")
    ax.legend(fontsize=7, loc="upper right")
    _save(fig, out_dir, "fig_pr_primary")


def fig_tau_ratio(cfg, out_dir, oof, evaluation, imp, attribution):
    """Horizon dependence: ratio to chance only, for the same six-way
    comparison as the PR figure.

    A longer horizon admits more crossings into every window, so the base rate
    rises with tau and the raw average precision rises with it for free. Only
    AP divided by the base rate of the same frame set can be read across
    horizons, which is why no absolute metric appears on this axis.
    """
    from .evaluate import ap_with_base
    fig, ax = plt.subplots(figsize=(7, 4.5))
    taus = cfg["labels.tau_ps"]
    sigma = _sigma_frame_scores(cfg)
    for k, (arm, label, color) in enumerate(COMPARE_ARMS):
        xs, ys = [], []
        for tau in taus:
            if arm == "published_rao2019_sigma":
                sub = _arm_frames(oof, sigma, arm, tau)
                if sub is None:
                    continue
                s = sub["score_raw"].to_numpy(float)
                m = np.isfinite(s)
                e = ap_with_base(sub["y"].to_numpy(int)[m], s[m])
                r = e.get("ratio_to_chance") if e.get("defined") else None
            else:
                e = evaluation["per_arm"].get(arm, {}).get(str(tau), {}).get("pooled", {})
                r = e.get("ratio_to_chance") if e.get("defined") else None
            if r is not None:
                xs.append(tau * (1 + 0.02 * (k - len(COMPARE_ARMS) / 2)))  # dodge points horizontally
                ys.append(r)
        if xs:
            ax.plot(xs, ys, "o-", label=label, color=color, alpha=0.8)
    if not ax.lines:
        plt.close(fig)
        return
    ax.axhline(1.0, ls="--", c="gray", label="chance level")
    ax.set_xscale("log")
    ax.set_xlabel("Horizon τ, ps (log scale)")
    ax.set_ylabel("AP / base rate")
    ax.set_title("Horizon dependence: the raw metric is incomparable, the ratio is comparable")
    ax.legend(fontsize=7)
    _save(fig, out_dir, "fig_tau_ratio")


def _first_fitted_arm(oof, cfg, tau, mask=None):
    """The head arm when it produced out-of-fold probabilities, else the first
    fitted branch. On event-starved systems the head refuses degenerate folds,
    and the figures then show whichever branch was fitted, labelled with its
    name."""
    head = cfg["model.head"]
    for arm in (head, "linear", "clock"):
        sub = oof[(oof["arm"] == arm) & (oof["tau"] == tau)]
        if mask is not None:
            sub = sub[mask(sub)]
        if len(sub) and np.isfinite(sub["prob_cal"].to_numpy(float)).any():
            note = ("" if arm == head else
                    f" [{head} refused: degenerate folds; showing {arm}]")
            return arm, sub, note
    return None, None, ""


def fig_timeline(cfg, out_dir, oof, evaluation, imp, attribution):
    """Behaviour over time: readiness series plus true crossings.

    One held-out trajectory, time in nanoseconds. The red verticals sit at
    `t_exit_ps`, the moment a crossing completes, because that is the instant
    the label refers to: y(t)=1 means a completion falls in (t, t+tau], so a
    useful monitor rises to the LEFT of each line.
    """
    from .events import load_events
    primary = cfg["labels.primary_tau_ps"]
    ev = load_events(cfg)
    counts = ev.groupby(["condition", "replica"]).size()
    if len(counts) == 0:
        return
    # Prefer the most event-rich unit that a fitted arm covers with calibrated
    # probabilities; on event-starved systems the richest unit is the one whose
    # fold refused, so a covered unit is drawn instead and noted as such.
    best = counts.idxmax()
    units = [u for u, _n in counts.sort_values(ascending=False).items()]
    all_reps = {(cc["id"], rr["id"]) for cc in cfg["data.conditions"]
                for rr in cc["replicas"]}
    units += sorted(all_reps - set(units))
    arm = g = None
    for c, r in units:
        arm, g, note = _first_fitted_arm(
            oof, cfg, primary,
            mask=lambda s: (s["condition"] == c) & (s["replica"] == r))
        if arm is not None:
            break
    if arm is None:
        return
    if (c, r) != best:
        note += (f" [no fitted fold covers {best[0]}/{best[1]} "
                 f"(most crossings); showing {c}/{r}]")
    g = g.sort_values("time_ps")
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(g["time_ps"] / 1000, g["prob_cal"], lw=0.7, label=f"P(ready), {arm}")
    sel = ev[(ev["condition"] == c) & (ev["replica"] == r)]
    for tx in sel["t_exit_ps"]:
        ax.axvline(tx / 1000, c="crimson", alpha=0.35, lw=0.8)
    ax.set_xlabel("Time, ns")
    ax.set_ylabel("Readiness probability")
    ax.set_title(f"Held-out trajectory {c}/{r}: readiness series and true crossings "
                 f"(red){note}", fontsize=10)
    ax.legend(fontsize=7)
    _save(fig, out_dir, "fig_timeline")


def fig_calibration(cfg, out_dir, oof, evaluation, imp, attribution):
    """Calibration of the reported probability; axes span the observed range.

    Bins are deciles of the predicted probability, not a fixed grid on [0, 1]:
    on a system with a base rate of a few per cent every prediction lands in the
    lowest fixed bin and the plot says nothing. Equal-count bins put the same
    number of frames behind each point, so the scatter of the points is
    comparable along the curve.
    """
    primary = cfg["labels.primary_tau_ps"]
    arm, sub, note = _first_fitted_arm(oof, cfg, primary,
                                       mask=lambda s: s["valid"])
    if arm is None:
        return
    p = sub["prob_cal"].to_numpy()
    y = sub["y"].to_numpy()
    m = np.isfinite(p)
    p, y = p[m], y[m]
    if len(p) == 0:
        return
    qs = np.quantile(p, np.linspace(0, 1, 11))
    # Bins are half-open (a, b]; nudging the first edge below the minimum keeps
    # the smallest probability inside bin 1 instead of falling out of every bin.
    qs[0] -= 1e-9
    fig, ax = plt.subplots(figsize=(5.5, 5))
    xs, ys = [], []
    for a, b in zip(qs[:-1], qs[1:]):
        mm = (p > a) & (p <= b)
        if mm.sum():
            xs.append(float(np.mean(p[mm])))
            ys.append(float(np.mean(y[mm])))
    ax.plot(xs, ys, "o-", label="observed positive fraction")
    lo = min(min(xs), min(ys))
    hi = max(max(xs), max(ys))
    pad = 0.05 * (hi - lo + 1e-9)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "--", c="gray",
            label="ideal calibration")  # the diagonal is drawn last
    ax.set_xlabel("Predicted probability")
    ax.set_ylabel("Positive fraction")
    ax.set_title(f"Calibration {arm}, τ={int(primary)} ps — ZOOMED to the observed "
                 f"range [{lo:.3f}, {hi:.3f}], not the full square{note}",
                 fontsize=9)
    ax.legend(fontsize=7)
    _save(fig, out_dir, "fig_calibration")


def fig_fold_spread(cfg, out_dir, oof, evaluation, imp, attribution):
    """Per-fold spread behind the pooled mean; horizontal, sorted by median."""
    primary = str(cfg["labels.primary_tau_ps"])
    rows = []
    for arm, taus in evaluation["per_arm"].items():
        e = taus.get(primary, {})
        for fold, r in e.get("per_fold", {}).items():
            if r.get("defined"):
                rows.append((arm, r["ratio_to_chance"]))
    if not rows:
        return
    df = pd.DataFrame(rows, columns=["arm", "ratio"])
    order = df.groupby("arm")["ratio"].median().sort_values().index.tolist()
    fig, ax = plt.subplots(figsize=(7, 0.5 * len(order) + 1.5))
    for i, arm in enumerate(order):
        vals = df[df["arm"] == arm]["ratio"]
        ax.plot(vals, [i] * len(vals), "o", alpha=0.65)
        ax.plot(vals.median(), i, "k|", ms=14)
    ax.axvline(1.0, ls="--", c="gray")  # chance line: ratio 1 is no signal
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order, fontsize=8)  # horizontal labels, sorted by median
    ax.set_xlabel("AP / base rate (a point below 1 is inverted ranking, not a weak signal)")
    ax.set_title(f"Per-fold spread, τ={primary} ps")
    _save(fig, out_dir, "fig_fold_spread")


def fig_block_importance(cfg, out_dir, oof, evaluation, imp, attribution):
    """Descriptor attribution; the label names the attribution kind used.

    Bars are sums over the columns of a block of the mean absolute per-frame
    contribution, in logits. Blocks rather than columns because near-duplicate
    columns inside one block share their credit arbitrarily, while the block
    total is stable; the axis label carries the attribution method so a
    TreeSHAP figure is never mistaken for a coefficient one.
    """
    per_block = imp.get("per_block", {})
    if not per_block:
        return
    items = sorted(per_block.items(), key=lambda kv: kv[1])
    fig, ax = plt.subplots(figsize=(7, 0.4 * len(items) + 1.5))
    ax.barh([k for k, _ in items], [v for _, v in items])
    ax.set_xlabel(f"Σ of mean |contributions| over block columns, logits; "
                  f"attribution: {attribution.get('kind', '?')}")
    ax.set_title("Block importance (aggregated to blocks)")
    _save(fig, out_dir, "fig_block_importance")


FIGURES = {
    "fig_pr_primary": fig_pr_primary,
    "fig_tau_ratio": fig_tau_ratio,
    "fig_timeline": fig_timeline,
    "fig_calibration": fig_calibration,
    "fig_fold_spread": fig_fold_spread,
    "fig_block_importance": fig_block_importance,
}


def run_step(cfg: Config) -> None:
    """Draw the whole catalog into the run's figures/ directory.

    catalog.json records the declared list next to the list actually produced.
    The two differ legitimately: a figure whose ingredients are missing on an
    event-starved system returns without saving. The report step compares the
    two lists, so a silently absent figure becomes a failed acceptance gate
    rather than an omission nobody notices.
    """
    data = _load(cfg)
    with step_output(cfg, "figures") as out:
        log = StepLog(out)
        for name, fn in FIGURES.items():
            fn(cfg, out, *data)
        produced = sorted(p.name for p in out.glob("*.png"))
        (out / "catalog.json").write_text(json.dumps(
            {"catalog": sorted(FIGURES), "produced": produced}, indent=1, ensure_ascii=False))
        log.say(f"figures: {produced}")
        log.close()
