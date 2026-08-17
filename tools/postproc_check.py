"""Measure the fast post-processing candidates on STORED scores — no refits.

Candidates: causal (backward) score smoothing; GBT+Linear rank ensemble;
multi-horizon rank aggregation; threshold-rule headroom (max-MCC vs the
F1 rule). Each is scored exactly like the thesis headline: average precision
over the base rate, on valid frames. Zero-shot smoothing is measured on the
canonical generalisation scores as well, per held protein.

    python tools/postproc_check.py         # writes runs/postproc-check/results.json

The point of measuring these before adopting any of them is that all four look
like free improvements and are cheap to apply, so the only honest ground for
leaving them out is a number. A frame-level gain here is necessary but not
sufficient: smoothing in particular was rejected on the event-level test at
matched time-in-warning, which lives in event_level_check.py.

Requirements: the committed out-of-fold tables of the three trainable systems and
the committed cross-protein score tables. Nothing is refitted and no trajectory
is opened, so this runs on the published tree with no trajectory data present.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pcm2 import config  # noqa: E402
from pcm2.evaluate import ap_with_base  # noqa: E402
from pcm2.runtime import run_dir  # noqa: E402
from sklearn.metrics import matthews_corrcoef  # noqa: E402

SHOW = {"gramicidin": 2000.0, "mthk": 2000.0, "kcsa": 2000.0}
SMOOTH_PS = [250.0, 500.0, 1000.0]
OUT = ROOT / "runs" / "postproc-check"


def ratio(y, s):
    """AP over the base rate of these frames, or None when it is undefined.

    The single figure of merit used throughout this script, so every candidate is
    compared on the metric the thesis reports rather than on one that flatters it.
    """
    m = np.isfinite(s)
    e = ap_with_base(y[m].astype(int), s[m])
    return round(e["ratio_to_chance"], 3) if e.get("defined") else None


def causal_smooth(df, col, window_ps):
    """Backward rolling mean of one column, per trajectory, in an array aligned to df.

    Smoothing is per trajectory because a window spanning two trajectories would
    average across a discontinuity in time, and strictly backward because an online
    monitor cannot see later frames. window_ps is converted to a frame count through
    each trajectory's own median frame spacing, so the physical window is the same
    on systems written at different strides. The caller must pass a df with a
    positional 0..n-1 index — the results are scattered back by index.
    """
    out = np.full(len(df), np.nan)
    for _u, g in df.groupby(["condition", "replica"], sort=False):
        g = g.sort_values("time_ps")
        dt = float(np.median(np.diff(g["time_ps"]))) if len(g) > 1 else 50.0
        k = max(1, int(round(window_ps / dt)))
        sm = g[col].rolling(k, min_periods=1).mean()  # backward-only window
        out[g.index.to_numpy()] = sm.to_numpy()
    return out


def pct_rank(s):
    """Percentile ranks, NaN preserved.

    Two arms are combined by averaging ranks rather than scores: a boosted logit
    and a ridge-logistic probability live on different scales, and averaging them
    directly would let whichever has the wider spread decide the ensemble.
    """
    r = pd.Series(s).rank(pct=True).to_numpy()
    r[~np.isfinite(s)] = np.nan
    return r


def main():
    """Measure every candidate on both settings and write results.json."""
    OUT.mkdir(exist_ok=True)
    report: dict = {"in_system": {}, "zero_shot_smoothing": {}}

    for sysid, tau in SHOW.items():
        cfg = config.load(f"configs/{sysid}.yaml")
        oof = pd.read_parquet(run_dir(cfg) / "train" / "oof.parquet")
        key = ["condition", "replica", "time_ps"]
        arms = {}
        for arm in ("trees", "linear"):
            a = oof[(oof["arm"] == arm) & (oof["tau"] == tau)].drop_duplicates(key)
            arms[arm] = a.sort_values(key).reset_index(drop=True)
        t = arms["trees"]
        v = t["valid"].to_numpy(bool)
        y = t["y"].to_numpy(int)
        base_score = t["score_raw"].to_numpy(float)
        entry = {"baseline_trees": ratio(y[v], base_score[v])}

        for w in SMOOTH_PS:
            sm = causal_smooth(t, "score_raw", w)
            entry[f"smooth_{int(w)}ps"] = ratio(y[v], sm[v])

        lin = arms["linear"]
        if len(lin) == len(t):
            ens = np.nanmean(np.stack([pct_rank(base_score),
                                       pct_rank(lin["score_raw"].to_numpy(float))]),
                             axis=0)
            entry["ensemble_trees+linear"] = ratio(y[v], ens[v])
            best_w = max(SMOOTH_PS,
                         key=lambda w: entry[f"smooth_{int(w)}ps"] or 0)
            td = t.copy()
            td["ens"] = ens
            ens_sm = causal_smooth(td, "ens", best_w)
            entry[f"ensemble_smoothed_{int(best_w)}ps"] = ratio(y[v], ens_sm[v])

        taus = sorted(oof[oof["arm"] == "trees"]["tau"].unique())
        ranks = []
        for tt in taus:
            a = (oof[(oof["arm"] == "trees") & (oof["tau"] == tt)]
                 .drop_duplicates(key).sort_values(key).reset_index(drop=True))
            if len(a) == len(t):
                ranks.append(pct_rank(a["score_raw"].to_numpy(float)))
        if len(ranks) > 1:
            entry[f"multi_tau_ranks({len(ranks)})"] = ratio(
                y[v], np.nanmean(np.stack(ranks), axis=0)[v])

        # Threshold headroom: max-MCC over the OOF sweep vs the stored F1 rule.
        # This is an
        # upper bound, not a candidate rule: the maximum is taken on the same
        # frames it is
        # read from, so it says how much a better threshold could possibly buy,
        # and the
        # pipeline's threshold stays the one chosen inside each fold.
        p = t["prob_cal"].to_numpy(float)
        mfin = np.isfinite(p) & v
        if mfin.sum() and 0 < y[mfin].sum() < mfin.sum():
            cand = np.quantile(p[mfin], np.linspace(0.5, 0.999, 60))
            mccs = [matthews_corrcoef(y[mfin], (p[mfin] >= c).astype(int))
                    for c in cand]
            entry["mcc_headroom_max"] = round(float(np.max(mccs)), 3)
        report["in_system"][sysid] = entry
        print(sysid, json.dumps(entry))

    # Zero-shot: causal smoothing of the canonical generalisation scores.
    HELD = {"gramicidin": ("lopo_gramicidin", 2000),
            "mthk": ("lopo_mthk", 2000),
            "kcsa_family": ("lopo_kcsa_family", 2000),
            "cx43": ("negative_control", 4000)}
    for held, (exp, tau) in HELD.items():
        sc = pd.read_parquet(ROOT / "runs" / f"monitor-generalisation-tau{tau}"
                             / "scores.parquet")
        g = sc[(sc["experiment"] == exp) & (sc["variant"] == "full")]
        if held == "cx43":
            g = g[g["unit"].str.startswith("cx43:")]
        # The score table keys frames by unit, not by condition/replica. Naming
        # the unit as
        # the condition makes causal_smooth group by trajectory as it does
        # in-system; the
        # reset_index is what gives it the positional index it scatters results
        # back into.
        g = g.assign(condition=g["unit"],
                     replica="1").reset_index(drop=True)
        entry = {}
        for w in [None] + SMOOTH_PS:
            s = (g["prob"].to_numpy(float) if w is None
                 else causal_smooth(g, "prob", w))
            per = []
            for u, gu in g.groupby("unit"):
                idx = gu.index.to_numpy()
                m = gu["valid"].to_numpy(bool)
                r = ratio(gu["y"].to_numpy(int)[m], s[idx][m])
                if r is not None:
                    per.append(r)
            tag = "baseline" if w is None else f"smooth_{int(w)}ps"
            entry[tag] = round(float(np.median(per)), 3) if per else None
        report["zero_shot_smoothing"][held] = entry
        print("zero-shot", held, json.dumps(entry))

    (OUT / "results.json").write_text(json.dumps(report, indent=1))
    print("written:", OUT / "results.json")


if __name__ == "__main__":
    main()
