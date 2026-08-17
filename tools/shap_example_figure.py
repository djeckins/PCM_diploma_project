"""Figure: how a single READY verdict is explained (TreeSHAP -> mechanism axes).

Left panel: signed axis contributions (log-odds) for one READY frame of the
held-out gramicidin showcase unit, scored by the pooled three-protein monitor;
the top column behind each axis is annotated with its measured value.
Right panel: the same frame's snapshot quantities — what the monitor says is
happening in the pore right now. Output: runs/FINAL/figG_shap_example.png

    python tools/shap_example_figure.py

The frame is not chosen by hand: it is the highest-probability frame inside a
declared time window of the showcase unit, so the figure shows the monitor at its
most confident rather than at its most convenient. Contributions are TreeSHAP
values in log-odds, grouped into mechanism axes by feature block and therefore
additive: they sum, with the base level, to the score behind the printed P(ready).

Requirements: the notebooks/predict_lib.py helper, the pooled monitor at this
horizon (taken from the cache in runs/pooled-predict/, refitted from the
committed feature tables if the cache is absent), and the showcase system's own
feature table. No trajectory is opened.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "notebooks"))
sys.path.insert(0, str(ROOT / "tools"))

plt.rcParams.update({
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "figure.titlesize": 10, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.18,
    "grid.linewidth": 0.4, "axes.titlelocation": "left",
})

import predict_lib as pl  # noqa: E402

TAU = 2000.0
UNIT_CFG = "gramicidin"
COND_ID, REP_ID = "330K_1Msalt", "04"
WINDOW_NS = (140.0, 166.0)  # the showcase unit's late crossing cluster


def main(unit_cfg=None, cond_id=None, rep_id=None, window=None, out=None,
         title_note=None):
    """Draw the worked example; every argument overrides one module default.

    The overrides exist so the same figure can be produced for another protein,
    another unit or another time window without a second copy of this script.
    window is (t_lo, t_hi) in PICOSECONDS, unlike WINDOW_NS above, which records
    the same default window in nanoseconds.
    """
    global UNIT_CFG, COND_ID, REP_ID
    UNIT_CFG = unit_cfg or UNIT_CFG
    COND_ID = cond_id or COND_ID
    REP_ID = rep_id or REP_ID
    bundle = pl.load_pooled_model(tau_ps=TAU, variant="full")
    cfg = pl.resolve_config(UNIT_CFG)
    cond = next(c for c in cfg["data.conditions"] if c["id"] == COND_ID)
    rep = next(r for r in cond["replicas"] if r["id"] == REP_ID)
    # The showcase frame: the highest-P(ready) frame inside the unit's late
    # crossing cluster — the moment the monitor is most confident about.
    lo, hi = window or (140000.0, 166000.0)
    win = pl.rows_from_table(cfg, cond, rep, lo, hi)
    pw = pl.predict_rows(bundle, win)
    best = int(pw["p_ready"].idxmax())
    rows = win.iloc[[best]].reset_index(drop=True)
    pred = pw.iloc[[best]].reset_index(drop=True)
    reasons = pl.explain_rows(bundle, rows, top_k=7)[0]
    p = float(pred["p_ready"].iloc[0])
    verdict = "READY" if bool(pred["ready"].iloc[0]) else "NOT READY"

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(9.5, 3.6), gridspec_kw={"width_ratios": [1.25, 1]})
    names = [r["axis"] for r in reasons][::-1]
    vals = [r["logodds"] for r in reasons][::-1]
    cols = ["#1f5fa8" if v > 0 else "#b3593b" for v in vals]
    ax.barh(range(len(vals)), vals, color=cols, height=0.6)
    span = max(abs(v) for v in vals) or 1.0
    for k, (v, r) in enumerate(zip(vals, reasons[::-1])):
        cv = r["column_value"]
        note = f"{r['top_column']}" + (f" = {cv:.2f}" if cv is not None else "")
        ax.annotate(note, (v + 0.02 * span * (1 if v >= 0 else -1), k),
                    fontsize=5.5, color="#5a6270", va="center",
                    ha="left" if v >= 0 else "right")
    ax.set_xlim(min(0, min(vals)) - 0.35 * span, max(0, max(vals)) + 0.35 * span)
    ax.set_yticks(range(len(names)), names)
    ax.axvline(0, color="#22262b", lw=0.8)
    ax.set_xlabel("TreeSHAP contribution to the verdict, log-odds")
    ax.set_title(f"TreeSHAP by mechanism axis · P(ready) = {p:.2f} ({verdict})\n"
                 f"{UNIT_CFG} {COND_ID}/{REP_ID} · "
                 f"t = {float(rows['time_ps'].iloc[0])/1000:.1f} ns · τ = {int(TAU)} ps")

    snap = pl.snapshot_text(rows.iloc[0], {TAU: (p, bool(pred["ready"].iloc[0]))}).splitlines()
    ax2.axis("off")
    ax2.set_title("Snapshot readout, same frame", loc="left")
    ax2.text(0.0, 0.96, "\n".join(snap), family="monospace", fontsize=5.6,
             va="top", transform=ax2.transAxes)
    fig.tight_layout()
    out = out or ROOT / "runs" / "FINAL" / "figG_shap_example.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    print("written:", out, "| verdict:", verdict, "P:", round(p, 3))
    for r in reasons[:3]:
        print("  axis:", r["axis"], "weight:", round(r["weight"], 2),
              "top:", r["top_column"])


if __name__ == "__main__":
    main()
