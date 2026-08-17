"""Event-level test of score smoothing at MATCHED time-in-warning.

The reviewer's requirement (firing-power literature, Snyder-style chance
control): a smoothing gain counts as skill only if, at the SAME fraction of
time spent in warning, the smoothed monitor anticipates more crossings than
the raw one — both compared against a structure-preserving random alarm of
identical duty. Everything runs on stored zero-shot scores; no refits.

    python tools/event_level_check.py      # writes runs/postproc-check/event_level_lopo4.json

Why matched duty is the whole point: smoothing raises the warned fraction simply
by keeping the alarm on longer, and an alarm that is armed half the time warns of
almost everything. Fixing the duty cycle first, by taking the threshold at the
(1 - duty) quantile of each series, removes that route and leaves only whether the
warnings sit in better places. This test is what rejected smoothing, which had
looked like a gain at frame level in postproc_check.py.

Requirements: the committed score tables of the four-protein-pool runs and the
committed events tables. Nothing is refitted and no trajectory is opened, so this
runs on the published tree with no trajectory data present.
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
from pcm2.evaluate import _random_null, alarm_metrics  # noqa: E402
from pcm2.runtime import run_dir, derive_seed  # noqa: E402

# (system, experiment in the score table, unit filter, horizon ps, smoothing
# window ps).
# A unit filter of None takes every unit of that experiment and averages over
# them; the
# other two name one unit each — the conducting E71A arm of the KcsA family, and
# the only
# Cx43 unit that carries positive frames.
CASES = [
    ("cx43", "lopo_cx43", "cx43:200mV/R2", 4000, 500.0),
    ("gramicidin", "lopo_gramicidin", None, 2000, 250.0),
    ("kcsa", "lopo_kcsa_family", "kcsa:E71A/1", 2000, 250.0),
]
# Fractions of total time the alarm is allowed to be armed. Spanning 2% to 20% shows
# whether an advantage survives across operating points or exists only at one of
# them.
DUTIES = [0.02, 0.05, 0.10, 0.20]


def unit_exits(sysid: str, unit: str) -> np.ndarray:
    """Sorted completion times, in ps, of the crossings this unit is judged against.

    Filtered to the condition's conduction direction, as the labels are: back
    crossings under an applied field are a different and rarer process, and counting
    them would credit the alarm for warning of something the label never marked.
    Direction 0 marks an event whose sign was never recorded and is kept.
    """
    cfg = config.load(f"configs/{sysid}.yaml")
    ev = pd.read_parquet(run_dir(cfg) / "events" / "events.parquet")
    cond, rep = unit.split(":")[1].split("/")
    sel = ev[(ev["condition"] == cond) & (ev["replica"] == rep)]
    direction = next(c.get("direction") for c in cfg["data.conditions"]
                     if c["id"] == cond)
    if direction is not None:
        sel = sel[(sel["direction"] == direction) | (sel["direction"] == 0)]
    return np.sort(sel["t_exit_ps"].to_numpy(float))


def smooth_series(times, p, window_ps):
    """Backward rolling mean over one trajectory; the window is given in ps.

    Converted to a frame count through this trajectory's own median spacing, so the
    physical window does not depend on the stride the system was written at.
    """
    dt = float(np.median(np.diff(times))) if len(times) > 1 else 50.0
    k = max(1, int(round(window_ps / dt)))
    return pd.Series(p).rolling(k, min_periods=1).mean().to_numpy()


def one_unit(times, p_raw, p_sm, exits, tau):
    """One row per duty target: warned fraction of raw, of smoothed, and of the null.

    A crossing counts as warned when the alarm was armed at some frame inside the
    tau window before its completion; the median lead time is reported beside it, in
    ps, over the crossings that were warned.
    """
    rows = []
    for duty in DUTIES:
        row = {"duty_target": duty}
        for tag, p in (("raw", p_raw), ("smoothed", p_sm)):
            # Threshold as a quantile of the series itself, so both series are
            # armed for
            # the same share of time whatever their scales are. This is the
            # matching that
            # the comparison rests on: a fixed probability threshold would give the
            # smoothed series, whose values are pulled toward the mean, a
            # different duty.
            thr = float(np.quantile(p[np.isfinite(p)], 1 - duty))
            armed = p >= thr
            m = alarm_metrics(times, armed, exits, tau)
            # The null keeps this alarm's episode structure and shifts it
            # cyclically in
            # time: same number of episodes, same lengths, same duty, positions
            # random.
            # That is the chance level an alarm has to beat, and it is not the
            # base rate —
            # a few long episodes warn of many crossings by covering the record.
            null = (_random_null(times, armed, exits, tau,
                                 derive_seed(20260816, "null", tag))
                    if armed.any() and len(exits) else
                    {"warned_frac_mean": float("nan")})
            row[f"{tag}_warned"] = round(m["warned_frac"], 3)
            row[f"{tag}_null"] = round(null["warned_frac_mean"], 3)
            row[f"{tag}_lead_med_ps"] = (round(float(np.median(m["lead_times_ps"])))
                                         if m.get("lead_times_ps") else None)
        rows.append(row)
    return rows


def main():
    """Run every case, print the per-duty table and write event_level_lopo4.json."""
    out = {}
    for sysid, exp, unit_filter, tau, window in CASES:
        sc = pd.read_parquet(ROOT / "runs" / f"monitor-generalisation-tau{tau}-cx43pool"
                             / "scores.parquet")
        g = sc[(sc["experiment"] == exp) & (sc["variant"] == "full")]
        units = ([unit_filter] if unit_filter
                 else sorted(u for u in g["unit"].unique()))
        agg: dict = {}
        for unit in units:
            gu = g[g["unit"] == unit].sort_values("time_ps")
            exits = unit_exits(sysid, unit)
            if not len(exits):
                continue
            times = gu["time_ps"].to_numpy(float)
            p = gu["prob"].to_numpy(float)
            rows = one_unit(times, p, smooth_series(times, p, window), exits, tau)
            for r in rows:
                d = agg.setdefault(r["duty_target"], {"raw": [], "smoothed": [],
                                                      "null": [], "n_units": 0})
                d["raw"].append(r["raw_warned"])
                d["smoothed"].append(r["smoothed_warned"])
                d["null"].append(np.nanmean([r["raw_null"], r["smoothed_null"]]))
                d["n_units"] += 1
        out[sysid] = {
            "tau_ps": tau, "window_ps": window, "n_units": len(units),
            "by_duty": {str(k): {"raw_warned_mean": round(float(np.mean(v["raw"])), 3),
                                 "smoothed_warned_mean": round(float(np.mean(v["smoothed"])), 3),
                                 "random_null_mean": round(float(np.nanmean(v["null"])), 3)}
                        for k, v in sorted(agg.items())}}
        print(sysid, json.dumps(out[sysid]["by_duty"]))
    d = ROOT / "runs" / "postproc-check"
    d.mkdir(exist_ok=True)
    (d / "event_level_lopo4.json").write_text(json.dumps(out, indent=1))
    print("written:", d / "event_level_lopo4.json")


if __name__ == "__main__":
    main()
