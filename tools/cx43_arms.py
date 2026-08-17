"""Comparison arms for connexin-43, which has no in-system model.

Both classes are present in only one of the two Cx43 assemblies, so the fold
that fits has nothing to rank. The published rules still apply frame by frame,
and the fitted arms can be applied zero-shot from the conducting pool, which is
how the monitor itself is scored on this protein.

    python tools/cx43_arms.py        # writes runs/FINAL/cx43_arms.json

The file it writes is read by tools/final_results.py, which fills the Cx43 group of
the method matrices from it and leaves that group empty when it is missing. Run this
first, then final_results.py.

Requirements: the committed Cx43 benchmark table and the committed feature and label
tables of the conducting systems, plus scikit-learn. The linear arms are fitted here;
no boosted model and no trajectory is involved. PLS-FMA has no entry at all — its
input is one protein's superimposed Cartesian coordinates, which have no counterpart
in another protein.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tools"))

from pcm2 import config                                   # noqa: E402
from pcm2.baselines.clock import clock_columns            # noqa: E402
from pcm2.baselines.linear_control import control_columns  # noqa: E402
from pcm2.evaluate import ap_with_base                    # noqa: E402
from pcm2.models.linear import fit_linear_fold            # noqa: E402
from pcm2.runtime import run_dir                          # noqa: E402
from monitor_generalisation import load_all               # noqa: E402

OUT = ROOT / "runs" / "FINAL" / "cx43_arms.json"
TAUS = (2000.0, 4000.0)
SEED = 20260816
POOL = ("gramicidin", "mthk")


def published(tau: float) -> dict:
    """The unfitted rules, scored in their published orientation.

    Read straight out of the Cx43 benchmark table, which was written at that run's
    primary horizon — 4000 ps — so this is only called for that tau; at any other
    horizon the labels in the table would not match the one asked for.
    """
    cfg = config.load("configs/cx43.yaml")
    b = pd.read_parquet(run_dir(cfg) / "benchmark" / "benchmark.parquet")
    y = b["truth_crossing_within_tau"].to_numpy(int)
    m = b["valid"].to_numpy(bool)
    out = {}
    for name, col in (("Rao 2019 energy", "published_rao2019_energy_kJmol"),
                      ("Rao 2019 sigma", "published_rao2019_sigma_score")):
        if col not in b.columns:
            continue
        # Both published quantities grow toward the non-conducting state: the
        # dewetting
        # energy above its ~1 RT contour, and sigma_d with the contour distance
        # past that
        # line. Negating puts them in the same orientation as every fitted arm,
        # where a
        # high score means ready, so one metric can be applied to all of them.
        s = -b[col].to_numpy(float)          # high score = wetted = ready
        k = m & np.isfinite(s)
        r = ap_with_base(y[k], s[k])
        out[name] = None if r.get("ratio_to_chance") is None else round(r["ratio_to_chance"], 3)
    return out


def zero_shot(tau: float) -> dict:
    """Fitted arms trained on the conducting pool and applied to Cx43.

    Scored on cx43:200mV/R2 alone, the only Cx43 unit that carries positive frames;
    with no positive there is nothing to rank and the metric would be undefined.
    Empty when this horizon is not declared in the label tables, or when the training
    rows carry no positive at all.
    """
    pool, _schema = load_all()
    ycol, vcol = f"y_tau{int(tau)}", f"valid_tau{int(tau)}"
    if ycol not in pool.columns:
        return {}
    train = pool["system"].isin([p for p in POOL]) | (
        (pool["system"] == "kcsa") & (pool["condition"] == "E71A"))
    target = (pool["system"] == "cx43") & (pool["replica"] == "R2")
    keep = pool[vcol].fillna(False).to_numpy(bool)
    # Keys, labels and the bl_ baseline columns are excluded: the bl_ columns
    # hold the
    # published rules' own outputs, and feeding them to a fitted arm would make
    # it a model
    # of those rules rather than an independent comparison.
    model_cols = [c for c in pool.columns
                  if c not in ("system", "condition", "replica", "time_ps", "unit")
                  and not c.startswith(("y_tau", "valid_tau", "bl_"))]
    Xtr = pool.loc[train & keep, model_cols].to_numpy(np.float64)
    ytr = pool.loc[train & keep, ycol].to_numpy(int)
    Xte = pool.loc[target & keep, model_cols].to_numpy(np.float64)
    yte = pool.loc[target & keep, ycol].to_numpy(int)
    if not len(yte) or ytr.sum() == 0:
        return {}

    arms = {"Clock (temporal null)": clock_columns(model_cols),
            "Linear control (structural)": control_columns(model_cols),
            "Linear (ridge-logistic)": model_cols}
    out = {}
    # The inner rotation for selecting C runs over training UNITS, so no penalty
    # is chosen
    # on frames of the assembly it was fitted on. Cx43 takes no part in that
    # selection.
    assemblies = pool.loc[train & keep, "unit"].to_numpy()
    for name, cols in arms.items():
        idx = [model_cols.index(c) for c in cols if c in model_cols]
        if not idx:
            continue
        names = [model_cols[i] for i in idx]
        # Linear-arm settings are borrowed from a per-system config rather than
        # restated,
        # so this comparison uses the same C grid, tolerance and class weighting
        # as the
        # linear arm of every pipeline run.
        lin = config.load("configs/gramicidin.yaml")["model.linear"]
        res = fit_linear_fold(Xtr[:, idx], ytr, assemblies, names, lin["C_grid"],
                              lin["tol"], lin["max_iter"], lin["class_weight_mode"], 0.999)
        p = res.pipeline.predict_proba(Xte[:, idx])[:, 1]
        r = ap_with_base(yte, p)
        out[name] = None if r.get("ratio_to_chance") is None else round(r["ratio_to_chance"], 3)
    return out


def main() -> None:
    """Write cx43_arms.json: one entry per horizon, one number per arm."""
    result = {"note": "Rao arms are unfitted rules on all valid Cx43 frames; the fitted arms are "
                      "trained on the three conducting proteins and applied to cx43:200mV/R2, the "
                      "only Cx43 unit carrying positives. PLS-FMA has no entry: its input is the "
                      "superimposed Cartesian coordinates of one protein, so it cannot transfer.",
              "seed": SEED, "by_tau": {}}
    for tau in TAUS:
        row = published(tau) if tau == 4000.0 else {}
        row.update(zero_shot(tau))
        result["by_tau"][str(int(tau))] = row
        print(int(tau), row)
    OUT.write_text(json.dumps(result, indent=1))
    print("written:", OUT)


if __name__ == "__main__":
    main()
