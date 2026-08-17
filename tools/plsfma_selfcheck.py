"""PLS-FMA exercised in its published regime: a continuous structural quantity.

Every application in Krivobokova et al. 2012 models a continuous scalar: a
channel-opening distance, a solvent-accessible surface, a permeability. The
benchmark arm of this project hands the method a binary readiness label, which
is outside that scope. Here the same PLS core models a continuous quantity of
this project's systems — the pore constriction radius, the analog of Aqy1's
channel-opening distance — validated on held-out trajectory units.

Measured quantity: the Pearson correlation between model and data on the
held-out units, per component count, with the correlation of the benchmark
arm's binary task available for comparison.

    python tools/plsfma_selfcheck.py       # writes runs/plsfma-selfcheck/

The purpose is to separate two possible readings of a weak PLS-FMA score in the
benchmark: a faulty port, or a method used outside its scope. A high held-out
correlation here settles it in favour of the second, and the benchmark number
becomes a statement about the forecasting task rather than about the code.

Requirements: coords/coords.parquet of each system, written by the coords step,
which streams the trajectories. Those tables are the heaviest artifact of the
project and are NOT part of the published repository, so on a fresh clone every
system is reported as skipped and no output is produced. Regenerate them with
`python -m pcm2.run coords --config configs/<system>.yaml` once the trajectory
data is in place.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pcm2 import config  # noqa: E402
from pcm2.models.plsfma import pls_denham  # noqa: E402
from pcm2.runtime import PROJECT_ROOT, run_dir  # noqa: E402
from pcm2.splits import inner_rotation  # noqa: E402

TARGET = "geo_r_constriction_A"
K_MAX = 40  # the paper's applications used up to 30-40 components
OUT = PROJECT_ROOT / "runs" / "plsfma-selfcheck"


def selfcheck(cfg_path: str) -> dict | None:
    """Component-count curve for one system, or None when its coords table is absent.

    Returns the frame and coordinate counts, the number of validation units and the
    curve: for each component count k, R_m on the rows the model was built from and
    R_c on the held-out units. R_c is the criterion — R_m rises with k by
    construction, so only the held-out correlation can say where the mode stops
    being structure and starts being noise. Validation rotates over whole
    trajectory units, never over frames, since adjacent frames are near-duplicates
    and a frame-level split would report a correlation that means nothing.

    Also writes selfcheck_<system>.png, the same curve with the chosen k marked.
    """
    cfg = config.load(cfg_path)
    root = run_dir(cfg)
    cpath, fpath = root / "coords" / "coords.parquet", root / "features" / "features.parquet"
    if not (cpath.exists() and fpath.exists()):
        return None
    coords = pd.read_parquet(cpath)
    feats = pd.read_parquet(fpath)[["condition", "replica", "time_ps", TARGET]]
    m = coords.merge(feats, on=["condition", "replica", "time_ps"], validate="one_to_one")
    # A frame whose constriction radius could not be measured has no target value to
    # model. Dropping such frames is right here and would be wrong in the
    # pipeline, where
    # a missing measurement is itself informative; this is a regression on the
    # radius.
    m = m[np.isfinite(m[TARGET])]
    ccols = [c for c in m.columns if c.startswith("ca")]
    X = m[ccols].to_numpy(np.float64)
    # Union-aligned coords: columns absent from one condition's topology come
    # in as missing; the analysis runs on the shared C-alpha set.
    keep = np.isfinite(X).all(axis=0)
    n_dropped = int((~keep).sum())
    X = X[:, keep]
    f = m[TARGET].to_numpy(np.float64)
    units = (m["condition"] + "/" + m["replica"]).to_numpy()
    inner = inner_rotation(units)

    curve = []
    # At most min(n_frames, n_coordinates) - 1 components exist; beyond that the
    # deflation has nothing left to extract, so the scan stops there even if
    # K_MAX allows
    # more.
    for k in range(1, min(K_MAX, min(X.shape) - 1) + 1):
        r_m, r_c = [], []
        for tr, va, _g in inner:
            model = pls_denham(X[tr], f[tr], k)
            pm, pv = model.predict(X[tr]), model.predict(X[va])
            # A constant prediction has zero variance and its correlation with
            # anything is
            # 0/0. Such a rotation is left out of the mean rather than entering
            # it as a
            # zero, which would read as a measured absence of correlation.
            if np.std(pm) > 0:
                r_m.append(float(np.corrcoef(pm, f[tr])[0, 1]))
            if np.std(pv) > 0:
                r_c.append(float(np.corrcoef(pv, f[va])[0, 1]))
        curve.append({"components": k,
                      "R_m": float(np.mean(r_m)) if r_m else None,
                      "R_c": float(np.mean(r_c)) if r_c else None})
    rc = [c["R_c"] for c in curve]
    best = int(np.argmax(rc))
    result = {"system": cfg["system.id"], "target": TARGET,
              "n_frames": int(len(f)), "n_coordinates": int(X.shape[1]),
              "n_union_columns_dropped": n_dropped,
              "n_units": int(len(set(units))), "curve": curve,
              "best_components": curve[best]["components"],
              "best_R_c": rc[best]}

    fig, ax = plt.subplots(figsize=(5.2, 3.4), dpi=120)
    ks = [c["components"] for c in curve]
    ax.plot(ks, [c["R_m"] for c in curve], "o-", ms=3.5, lw=1.4,
            color="#9ca3af", label="$R_m$ (training)")
    ax.plot(ks, rc, "o-", ms=3.5, lw=1.6, color="#2563eb",
            label="$R_c$ (held-out units)")
    ax.axvline(curve[best]["components"], color="#2563eb", lw=0.8, ls=":")
    ax.set_xlabel("PLS components")
    ax.set_ylabel("Pearson correlation")
    ax.set_title(f"{cfg['system.id']}: PLS-FMA on {TARGET}", fontsize=10)
    ax.set_ylim(0, 1.02)
    ax.grid(color="#e5e7eb", lw=0.5)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(fontsize=8, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(OUT / f"selfcheck_{cfg['system.id']}.png")
    plt.close(fig)
    return result


def main() -> None:
    """Run the self-check on every system that has a coords table and write the summary.

    A system without one is reported and skipped, not treated as a failure: the
    coords step is the one soft dependency of the pipeline and its tables are absent
    from the published tree.
    """
    OUT.mkdir(parents=True, exist_ok=True)
    results = []
    for p in ["configs/gramicidin.yaml", "configs/mthk.yaml",
              "configs/kcsa.yaml", "configs/kcsa_na.yaml", "configs/cx43.yaml"]:
        if not Path(p).exists():
            continue
        r = selfcheck(p)
        if r is None:
            print(f"{p}: tables not on disk yet — skipped")
            continue
        results.append(r)
        print(f"{r['system']:>10}: best R_c = {r['best_R_c']:.3f} at "
              f"{r['best_components']} components "
              f"({r['n_frames']} frames, {r['n_coordinates']} coords, "
              f"{r['n_units']} units)")
    (OUT / "selfcheck.json").write_text(json.dumps(results, indent=1))
    lines = ["# PLS-FMA self-check on a continuous functional quantity", "",
             "The published method models continuous structural scalars; here the "
             f"same PLS core models `{TARGET}` with validation on held-out "
             "trajectory units. High $R_c$ means the implementation delivers what "
             "the method promises in its own regime; the benchmark's readiness "
             "scores are then a statement about the forecasting task.", ""]
    lines += ["| system | frames | coords | units | best $R_c$ | components |",
              "|---|---|---|---|---|---|"]
    for r in results:
        lines.append(f"| {r['system']} | {r['n_frames']} | {r['n_coordinates']} | "
                     f"{r['n_units']} | {r['best_R_c']:.3f} | {r['best_components']} |")
    (OUT / "SELFCHECK.md").write_text("\n".join(lines) + "\n")
    print(f"written: {OUT}")


if __name__ == "__main__":
    main()
