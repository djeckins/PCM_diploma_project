"""Pooled leave-one-channel-out: train on several systems, test on an unseen one.

Trains a single model on the pooled descriptor tables of all available systems
minus one held-out unit (a KcsA filter mutant), and evaluates on that unit.
This is the cross-architecture generalisation experiment; pooling is possible
because the schema is a union over column names.

Constraints:
  * the common horizon is the only tau computed for every system;
  * the headline metric is ranking WITHIN the held-out unit against that unit's
    own base rate. Pooled AP across systems whose base rates differ by 5x is
    not comparable; within-unit ranking is insensitive to cross-system score
    offsets;
  * columns whose per-system value ranges do not even overlap are reported
    before training: they are potential system labels, and a pooled model
    could exploit them to memorise the system instead of the physics;
  * a coordinate-based arm (PLS-FMA) is impossible here by construction, since
    the Cartesian spaces of different proteins are incommensurable.

Reads only published artifacts of the per-system runs (features, labels,
schema); writes its own results directory atomically.

    python tools/pooled_loco.py            # writes runs/pooled-loco-tau1000/

This is the earlier, unit-level version of the cross-protein experiment and has
been superseded by monitor_generalisation.py, which holds out whole PROTEINS
rather than single units and trains only where conduction exists. It is kept
because its output directory is committed evidence of what was measured before
that change; nothing in the write-up depends on rerunning it.

Requirements: the committed features, labels and schema artifacts of the four
systems it names, plus xgboost and scikit-learn. Models ARE fitted here — one
boosted and two linear arms per held-out unit — so the run costs real compute. No
trajectory is opened.
"""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pcm2 import config  # noqa: E402
from pcm2.evaluate import ap_with_base  # noqa: E402
from pcm2.features._common import ColSpec  # noqa: E402
from pcm2.models import Calibrator, f1_threshold  # noqa: E402
from pcm2.models.linear import fit_linear_fold  # noqa: E402
from pcm2.models.trees import fit_trees_fold, trees_predict  # noqa: E402
from pcm2.baselines.clock import clock_columns  # noqa: E402
from pcm2.runtime import PROJECT_ROOT, run_dir, write_provenance  # noqa: E402

TAU_PS = 1000.0  # the only horizon computed for every system
SYSTEM_CONFIGS = ["configs/gramicidin.yaml", "configs/mthk.yaml",
                  "configs/kcsa.yaml", "configs/kcsa_na.yaml"]
# Units held out one at a time: every KcsA filter variant, in both ion arms.
HELD_OUT_PREFIXES = ("kcsa:", "kcsa_na:")

# Search only the parameters the per-system flatness measurements left open;
# row/column subsampling was criterion-flat on both measured systems.
DEPTH_GRID, K_GRID = [2, 3, 4], [2.0, 4.0, 8.0]
COLSAMPLE, SUBSAMPLE = [0.7], [0.8]


def load_system(cfg_path: str) -> tuple[str, pd.DataFrame, list[ColSpec]]:
    """(system id, features joined to labels, schema) for one system.

    The unit key is built from the replica's LINEAGE, not from its replica id: two
    replicas that descend from one equilibrated assembly are not independent, and
    holding out one of them while training on the other would leak. Where each
    replica is its own assembly the two coincide.
    """
    cfg = config.load(cfg_path)
    root = run_dir(cfg)
    feats = pd.read_parquet(root / "features" / "features.parquet")
    labels = pd.read_parquet(root / "labels" / "labels.parquet")
    doc = json.loads((root / "features" / "schema.json").read_text())
    schema = [ColSpec(**d) for d in doc["columns"]]
    merged = feats.merge(labels, on=["condition", "replica", "time_ps"],
                         validate="one_to_one")
    merged.insert(0, "system", cfg["system.id"])
    lineage = {}
    for cond in cfg["data.conditions"]:
        for rep in cond["replicas"]:
            lineage[(cond["id"], rep["id"])] = f"{cfg['system.id']}:{rep['lineage']}"
    merged["unit"] = [lineage[(c, r)] for c, r in
                      zip(merged["condition"], merged["replica"])]
    return cfg["system.id"], merged, schema


def main() -> None:
    """Hold out each KcsA filter-mutant unit in turn and write the results directory."""
    tables, schemas = [], {}
    for p in SYSTEM_CONFIGS:
        sysid, t, schema = load_system(p)
        tables.append(t)
        schemas[sysid] = [c.name for c in schema]
        print(f"{sysid}: {len(t)} rows")
    names0 = list(schemas.values())[0]
    # Pooling only means anything if the column lists are identical AND in the
    # same order:
    # the model matrix is built by position, so a permuted schema would silently
    # feed each
    # system's values into another's column. Compared as lists for exactly that
    # reason.
    for sysid, names in schemas.items():
        assert names == names0, f"{sysid}: schema names differ — union broken"
    pool = pd.concat(tables, ignore_index=True)

    schema = [ColSpec(**d) for d in json.loads(
        (run_dir(config.load(SYSTEM_CONFIGS[0])) / "features" / "schema.json")
        .read_text())["columns"]]
    model_cols = [c.name for c in schema if c.to_model]

    # System-label screen: columns whose per-system ranges do not overlap. Such
    # a column
    # can be thresholded to recover which system a frame came from, and a pooled
    # model is
    # free to use it to reproduce each system's base rate instead of measuring
    # physics.
    # Here the list is only reported, before training, so the reader knows what
    # the pooled
    # model had available; the later experiment refits without them and measures the
    # difference.
    sys_labels = []
    for c in model_cols:
        ranges = []
        for s, g in pool.groupby("system"):
            v = g[c].dropna()
            if len(v):
                ranges.append((float(v.min()), float(v.max())))
        if len(ranges) >= 2:
            lo = max(r[0] for r in ranges)
            hi = min(r[1] for r in ranges)
            if lo > hi:
                sys_labels.append(c)
    print(f"potential system-label columns (disjoint per-system ranges): {sys_labels}")

    ycol, vcol = f"y_tau{int(TAU_PS)}", f"valid_tau{int(TAU_PS)}"
    y = pool[ycol].to_numpy(int)
    valid = pool[vcol].to_numpy(bool)
    X = pool[model_cols].to_numpy(float)
    units = pool["unit"].to_numpy()
    held_units = sorted(u for u in np.unique(units)
                        if u.startswith(HELD_OUT_PREFIXES))
    seed = 20260814

    results: dict = {"tau_ps": TAU_PS, "held_out": {},
                     "system_label_columns": sys_labels,
                     "n_pool_rows": int(len(pool)),
                     "note_plsfma": "coordinate-based FMA cannot pool across proteins: "
                                    "Cartesian spaces are incommensurable"}
    oof_rows = []
    for held in held_units:
        # Training rows must be valid: a right-censored frame has an unknown
        # label, and
        # treating it as negative would teach the model that the end of every
        # trajectory is
        # quiet. Test rows are kept whole and the validity mask is applied when
        # scoring.
        te = np.flatnonzero(units == held)
        tr = np.flatnonzero((units != held) & valid)
        base = float(np.mean(y[te][valid[te]]))
        entry: dict = {"n_train": int(len(tr)), "n_test": int(len(te)),
                       "test_base_rate": base, "arms": {}}
        print(f"\n=== held out: {held} (train {len(tr)} rows, "
              f"test base rate {base:.3f}) ===")

        arms: dict[str, np.ndarray] = {}
        r_tr = fit_trees_fold(X[tr], y[tr], units[tr], units[tr], model_cols,
                              DEPTH_GRID, K_GRID, COLSAMPLE, SUBSAMPLE,
                              lr=0.05, ceiling=2000, patience=50, monotone_map={},
                              seed=seed, n_threads=1)
        cal = Calibrator("sigmoid", 200).fit(r_tr.inner_oof_scores, r_tr.inner_oof_y)
        thr = f1_threshold(cal.predict(r_tr.inner_oof_scores), r_tr.inner_oof_y)
        raw = trees_predict(r_tr.booster, X[te], model_cols)
        arms["trees"] = raw
        entry["arms"]["trees"] = {
            "chosen": r_tr.chosen, "M_rounds": r_tr.M,
            "eval": ap_with_base(y[te][valid[te]], raw[valid[te]]),
            "verdict_accuracy": float(np.mean(
                (cal.predict(raw)[valid[te]] >= thr) == y[te][valid[te]]))}

        lin = fit_linear_fold(X[tr], y[tr], units[tr], model_cols,
                              [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0],
                              1e-8, 5000, "measure", 0.98)
        arms["linear"] = lin.pipeline.predict_proba(X[te])[:, 1]
        entry["arms"]["linear"] = {
            "chosen_C": lin.chosen_C,
            "eval": ap_with_base(y[te][valid[te]], arms["linear"][valid[te]])}

        ccols = clock_columns(model_cols)
        csub = [model_cols.index(c) for c in ccols]
        clk = fit_linear_fold(X[np.ix_(tr, csub)], y[tr], units[tr], ccols,
                              [1e-4, 1e-3, 1e-2, 1e-1, 1.0, 10.0, 100.0],
                              1e-8, 5000, "measure", 0.98)
        arms["clock"] = clk.pipeline.predict_proba(X[np.ix_(te, csub)])[:, 1]
        entry["arms"]["clock"] = {
            "eval": ap_with_base(y[te][valid[te]], arms["clock"][valid[te]])}

        for arm, scores in arms.items():
            e = entry["arms"][arm].get("eval", {})
            if e.get("defined"):
                print(f"  {arm:8s} AP={e['ap']:.3f} (x{e['ratio_to_chance']:.2f} "
                      f"over the unit's own chance)")
        for k, i in enumerate(te):
            oof_rows.append((held, pool["condition"].iat[i], pool["time_ps"].iat[i],
                             float(arms["trees"][k]), float(arms["linear"][k]),
                             float(arms["clock"][k]), int(y[i]), bool(valid[i])))
        results["held_out"][held] = entry

    out_name = "pooled-loco-tau1000"
    final = PROJECT_ROOT / "runs" / out_name
    tmp = PROJECT_ROOT / "runs" / f".{out_name}.tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    (tmp / "results.json").write_text(json.dumps(results, indent=1,
                                                 ensure_ascii=False, default=str))
    pd.DataFrame(oof_rows, columns=["held_out", "condition", "time_ps",
                                    "trees_score", "linear_score", "clock_score",
                                    "y", "valid"]).to_parquet(tmp / "oof.parquet")
    lines = [f"# Pooled leave-one-channel-out (tau = {int(TAU_PS)} ps)", "",
             "Trained on the pooled tables of all systems minus one held-out KcsA "
             "filter mutant; evaluated by ranking within the held-out unit against "
             "its own base rate.", "",
             "| held-out unit | base rate | trees | linear | clock |", "|---|---|---|---|---|"]
    for held, e in results["held_out"].items():
        # One arm's cell: the ratio to the held unit's own chance, or an em dash
        # where the
        # metric is undefined because the unit carries a single class.
        def r(a):
            ev = e["arms"][a].get("eval", {})
            return f"x{ev['ratio_to_chance']:.2f}" if ev.get("defined") else "—"
        lines.append(f"| {held} | {e['test_base_rate']:.3f} | {r('trees')} | "
                     f"{r('linear')} | {r('clock')} |")
    lines += ["", f"Potential system-label columns (disjoint per-system ranges): "
              f"{', '.join(sys_labels) if sys_labels else 'none'}.",
              "", results["note_plsfma"] + "."]
    (tmp / "POOLED.md").write_text("\n".join(lines))
    write_provenance(tmp, config.load(SYSTEM_CONFIGS[-1]), "transfer",
                     extra={"experiment": "pooled leave-one-channel-out",
                            "systems": [str(p) for p in SYSTEM_CONFIGS]})
    if final.exists():
        shutil.rmtree(final)
    tmp.rename(final)
    print(f"\nwritten: {final}/POOLED.md")


if __name__ == "__main__":
    main()
