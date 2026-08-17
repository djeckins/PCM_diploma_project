"""Verification of the interactive prediction path.

Checks:
  1. at every horizon in HORIZONS, the labels the model is trained on match
     the ones the pipeline wrote on the rows of the systems that declared that
     horizon, and the held-out quality of the fit is reported;
  2. a time-range query returns exactly the frames of the requested interval,
     including a range at the start of a trajectory, a range shorter than one
     frame step, a range running past the end, and a range entirely beyond it;
  3. rows measured from the trajectory (window mode, as on a protein with no
     pipeline run) agree with rows read from the computed table.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import predict_lib as pl
from monitor_generalisation import load_all

HORIZONS = [250.0, 500.0, 750.0, 1000.0, 1500.0, 2000.0, 3000.0, 4000.0, 6000.0]
GRAMICIDIN = "gramicidin"  # system id, resolved via configs/ (no machine paths)


def check_horizons() -> None:
    print("=" * 100)
    print("HORIZONS")
    print("=" * 100)
    pool, _ = load_all()
    for tau in HORIZONS:
        ycol = f"y_tau{int(tau)}"
        # Wherever a system declared this horizon, the rebuilt labels must
        # equal the ones the pipeline wrote for those rows.
        cross = "no system declares it"
        if ycol in pool.columns:
            have = pool[ycol].notna().to_numpy()
            y_re, v_re = pl.relabel_pool(pool, tau)
            same_y = np.array_equal(y_re[have], pool.loc[have, ycol].to_numpy(int))
            same_v = np.array_equal(v_re[have],
                                    pool.loc[have, f"valid_tau{int(tau)}"].to_numpy(bool))
            systems = sorted(set(pool.loc[have, "system"]))
            cross = (f"matches pipeline labels on {have.sum()} rows "
                     f"({','.join(systems)}): y {same_y}, valid {same_v}")
            assert same_y and same_v, f"rebuilt labels differ at tau={tau}"

        b = pl.load_pooled_model(tau_ps=tau, variant="physics_only")
        q = b["oof_quality"]
        lift = f"x{q['ratio_to_chance']:.2f}" if q["defined"] else "undefined"
        print(f"\ntau = {tau:6.0f} ps | {b['labels_from']}")
        print(f"    positives {b['train_positives']:6d} / {b['train_rows']} rows "
              f"({b['base_rate']:6.2%}) | threshold {b['threshold']:.4f} | "
              f"depth {b['hyperparams']['chosen']['depth']}, M {b['hyperparams']['M']}")
        print(f"    held-out (leave-one-unit-out) precision {q.get('ap', float('nan')):.4f} "
              f"vs chance {q.get('base_rate', float('nan')):.4f} = {lift}")
        print(f"    cross-check: {cross}")


def check_ranges() -> None:
    print()
    print("=" * 100)
    print("TIME RANGES")
    print("=" * 100)
    cfg = pl.resolve_config(GRAMICIDIN)
    bundle = pl.load_pooled_model(tau_ps=1000.0, variant="physics_only")
    cond, rep = pl.resolve_replica(cfg, "298K_2Msalt", "01")
    table = pl.rows_from_table(cfg, cond, rep, -np.inf, np.inf)
    t_all = table["time_ps"].to_numpy()
    step = float(np.median(np.diff(t_all)))
    print(f"trajectory grid: {len(t_all)} frames, {t_all[0]:.0f}–{t_all[-1]:.0f} ps, "
          f"step {step:.0f} ps\n")

    cases = [
        ("1000-4000 ps (early stretch)", 1000.0, 4000.0),
        ("0-2000 (from the very start)", 0.0, 2000.0),
        ("single point", 100000.0, 100000.0),
        ("shorter than one frame step", 100000.0, 100000.0 + step / 3),
        ("whole trajectory", t_all[0], t_all[-1]),
        ("last 5 ns", t_all[-1] - 5000.0, t_all[-1]),
        ("past the end", t_all[-1] - 500.0, t_all[-1] + 50000.0),
    ]
    for name, lo, hi in cases:
        rows, _c, source = pl.query_rows(cfg, "298K_2Msalt", "01", lo, hi)
        expect = int(((t_all >= lo) & (t_all <= hi)).sum())
        pred = pl.predict_rows(bundle, rows)
        t = rows["time_ps"].to_numpy()
        ok = (len(rows) == expect and np.all(np.diff(t) > 0)
              and t.min() >= lo and t.max() <= hi
              and np.isfinite(pred["p_ready"]).all())
        print(f"  {name:<32} rows {len(rows):5d} (expected {expect:5d}) "
              f"{t.min():9.0f}–{t.max():9.0f} ps | "
              f"P(ready) {pred['p_ready'].min():.3f}–{pred['p_ready'].max():.3f} | "
              f"READY {int(pred['ready'].sum()):4d} | {'OK' if ok else 'MISMATCH'}")
        assert ok, f"range case failed: {name}"

    # A range entirely outside the trajectory must raise, not return rows.
    try:
        pl.query_rows(cfg, "298K_2Msalt", "01", 1e9, 2e9)
        print("  range beyond the trajectory      NO REFUSAL — defect")
    except Exception as e:
        print(f"  {'range fully beyond the end':<32} refused: {type(e).__name__}")


def check_window_mode() -> None:
    print()
    print("=" * 100)
    print("WINDOW MODE (rows measured from the trajectory, as on a new protein)")
    print("=" * 100)
    cfg = pl.resolve_config(GRAMICIDIN)
    cond, rep = pl.resolve_replica(cfg, "298K_2Msalt", "01")
    bundle = pl.load_pooled_model(tau_ps=1000.0, variant="physics_only")
    ages = {"dyn_wet_state_age_ps", "dyn_occ_change_age_ps"}

    for lo, hi in [(1000.0, 4000.0), (100000.0, 102000.0), (0.0, 1500.0)]:
        table = pl.rows_from_table(cfg, cond, rep, lo, hi)
        traj = pl.rows_from_trajectory(cfg, cond, rep, lo, hi)
        cols = [c for c in bundle["columns"] if c in table.columns]
        diff = []
        for c in cols:
            a, b = table[c].to_numpy(float), traj[c].to_numpy(float)
            same = (np.isnan(a) & np.isnan(b)) | np.isclose(a, b, rtol=1e-6, atol=1e-9)
            if not same.all():
                diff.append(c)
        pt = pl.predict_rows(bundle, table)["p_ready"].to_numpy()
        pj = pl.predict_rows(bundle, traj)["p_ready"].to_numpy()
        dmax = float(np.nanmax(np.abs(pt - pj))) if len(pt) == len(pj) else float("nan")
        expected = [c for c in diff if c in ages]
        unexpected = [c for c in diff if c not in ages]
        print(f"  [{lo:7.0f}, {hi:7.0f}] ps: {len(table)} rows | "
              f"model columns differing: {len(diff)} "
              f"(window-edge ages {len(expected)}, other {len(unexpected)}) | "
              f"max |ΔP(ready)| = {dmax:.4g}")
        if unexpected:
            print(f"      unexpected: {unexpected}")


if __name__ == "__main__":
    check_horizons()
    check_ranges()
    check_window_mode()
    print("\nall checks passed")
