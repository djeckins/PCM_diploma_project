"""Conducting-channels generalisation: train only where conduction exists.

Two experiments, both with the same unit roles:

  1. Leave-one-protein-out over the three CONDUCTING channels (gramicidin,
     MthK, KcsA-E71A): train on the other two proteins' units, predict the
     held-out protein's frames. This is the cross-architecture claim measured
     without event-starved units polluting the training pool.

  2. Negative-control application: train on all three conducting proteins,
     apply to the engineered non-conducting / rare arms (E71A-Na, G77A-K/Na,
     T75A-K/Na). A monitor that has never seen these channels should read them
     as almost never READY; the G77A-Na arm, with its residual Na leak, is the
     one rare-positive test where ranking is definable.

Each fit is run twice: on the full model-column set and on a physics-only set
that drops columns whose per-system value ranges do not overlap. Those columns
identify the system, so the difference between the two fits bounds how much of
the transfer they account for.

    PCM2_GEN_TAU=2000 python tools/monitor_generalisation.py

Writes runs/monitor-generalisation-tau<TAU>/: results.json (per-unit and median
AP ratios, the dropped column list, the chosen hyperparameters),
scores.parquet (per-frame probability, label and validity for every
out-of-training unit), one readiness timeline per unit, GENERALISATION.md and
PROVENANCE.json. Switches: PCM2_GEN_TAU picks the horizon; PCM2_POOL_ADD_CX43=1
runs the pool-inclusion ablation into a "-cx43pool" directory; PCM2_INNER set to
"protein" runs the rejected selection variant into a "-pinner" one.

Requirements: the committed features, labels, schema and events artifacts of all
five systems, plus xgboost. No trajectory is opened — every measurement was made
by the per-system pipeline — but models ARE fitted here, so this is the one
cross-protein step that costs real compute. It replaces its output directory
wholesale, so rerunning it overwrites the committed results with the results of
the machine it runs on.
"""

from __future__ import annotations

import os

import json
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pcm2 import config  # noqa: E402
from pcm2.evaluate import ap_with_base  # noqa: E402
from pcm2.labels import label_one_trajectory  # noqa: E402
from pcm2.models import Calibrator, f1_threshold  # noqa: E402
from pcm2.models.trees import fit_trees_fold, trees_predict  # noqa: E402
from pcm2.models.linear import fit_linear_fold  # noqa: E402
from pcm2.baselines.clock import clock_columns  # noqa: E402
from pcm2.runtime import PROJECT_ROOT, run_dir, write_provenance  # noqa: E402
from pcm2.features._common import ColSpec  # noqa: E402

TAU_PS = float(os.environ.get("PCM2_GEN_TAU", "1000"))
SYSTEM_CONFIGS = ["configs/gramicidin.yaml", "configs/mthk.yaml",
                  "configs/kcsa.yaml", "configs/kcsa_na.yaml",
                  "configs/cx43.yaml"]
# The three conducting channels. The G77A Na leak is held out of training:
# LOPO medians are unchanged by its 43 positive frames (gA x1.87->1.83,
# MthK x1.30->1.33 between pool variants), and holding it out keeps it
# available as a zero-shot test.
CONDUCTING = {("gramicidin", None), ("mthk", None), ("kcsa", "E71A")}
# Pool-inclusion experiment: PCM2_POOL_ADD_CX43=1 admits the cx43 units into
# the training pool, LOPO gains a fourth held protein, and the output moves to
# a "-cx43pool" directory so the canonical three-pool runs stay untouched.
POOL_ADD_CX43 = os.environ.get("PCM2_POOL_ADD_CX43") == "1"
if POOL_ADD_CX43:
    CONDUCTING = CONDUCTING | {("cx43", None)}
PROTEIN_OF = {"gramicidin": "gramicidin", "mthk": "mthk",
              "kcsa": "kcsa_family", "kcsa_na": "kcsa_family",
              "cx43": "cx43"}


def load_all():
    """Every system's descriptor table joined to its labels, stacked into one pool.

    Returns (pool, schema). Pooling across proteins is only meaningful because
    the feature schema is a union over column names, identical for all systems;
    the schema is taken from the first system read and the columns of the others
    line up with it by construction. Each row gains a "system" column and a
    "unit" key of the form "<system>:<condition>/<replica>", which is the level
    every split in this file is made at.
    """
    frames = []
    schema = None
    for p in SYSTEM_CONFIGS:
        cfg = config.load(p)
        root = run_dir(cfg)
        t = pd.read_parquet(root / "features" / "features.parquet")
        lab = pd.read_parquet(root / "labels" / "labels.parquet")
        m = t.merge(lab, on=["condition", "replica", "time_ps"], validate="one_to_one")
        m.insert(0, "system", cfg["system.id"])
        m["unit"] = m["system"] + ":" + m["condition"] + "/" + m["replica"]
        frames.append(m)
        if schema is None:
            schema = [ColSpec(**d) for d in json.loads(
                (root / "features" / "schema.json").read_text())["columns"]]
    return pd.concat(frames, ignore_index=True), schema


def is_conducting(row_system: str, row_condition: str) -> bool:
    """Whether this row may enter a training pool.

    A system-wide entry (system, None) admits all its conditions; a
    (system, condition) entry admits that condition only, which is how the
    conducting E71A arm of KcsA is separated from its non-conducting siblings in
    the same config.
    """
    return ((row_system, None) in CONDUCTING
            or (row_system, row_condition) in CONDUCTING)


def pool_labels(pool: pd.DataFrame, tau_ps: float) -> tuple[np.ndarray, np.ndarray]:
    """Labels at this horizon for every pool row.

    Taken straight from the label tables when every system declares the
    horizon; otherwise rebuilt from the recorded crossings with the pipeline's
    own labelling rule (same core as the prediction notebook, where the rebuild
    is bit-identical to the tables). The fallback exists because reading a
    partially declared column directly would turn unknown frames into
    negatives."""
    ycol, vcol = f"y_tau{int(tau_ps)}", f"valid_tau{int(tau_ps)}"
    if ycol in pool.columns and not pool[ycol].isna().any():
        return pool[ycol].to_numpy(int), pool[vcol].to_numpy(bool)
    y = np.zeros(len(pool), dtype=int)
    valid = np.zeros(len(pool), dtype=bool)
    for path in SYSTEM_CONFIGS:
        cfg = config.load(path)
        sysid = cfg["system.id"]
        ev = pd.read_parquet(run_dir(cfg) / "events" / "events.parquet")
        for cond in cfg["data.conditions"]:
            direction = cond.get("direction")
            for rep in cond["replicas"]:
                m = ((pool["system"] == sysid)
                     & (pool["condition"] == cond["id"])
                     & (pool["replica"] == rep["id"])).to_numpy()
                if not m.any():
                    continue
                times = pool.loc[m, "time_ps"].to_numpy()
                sel = ev[(ev["condition"] == cond["id"])
                         & (ev["replica"] == rep["id"])]
                if direction is not None:
                    sel = sel[(sel["direction"] == direction)
                              | (sel["direction"] == 0)]
                lab = label_one_trajectory(
                    times, sel["t_exit_ps"].to_numpy(),
                    sel.loc[sel["entrance_observed"], "t_entry_ps"].to_numpy(),
                    tau_ps, float(times[-1]))
                idx = np.flatnonzero(m)
                y[idx] = lab[ycol].to_numpy()
                valid[idx] = lab[vcol].to_numpy()
    return y, valid


# PCM2_INNER=protein makes hyperparameter selection leave-one-PROTEIN-out
# inside the training pool (train on one protein, validate on another),
# aligning the selection criterion with the cross-protein task; the default
# is a within-pool unit rotation. Either way the held-out protein is unseen.
INNER = os.environ.get("PCM2_INNER", "unit")


def fit_and_score(X, y, units, tr_mask, te_mask, cols_idx, names, seed,
                  inner_groups=None):
    """Fit one boosted monitor on tr_mask and score te_mask with it.

    Returns (raw scores, calibrated probabilities, decision threshold, chosen
    hyperparameters) for the rows selected by te_mask, in the order those rows
    appear in X. Only cols_idx columns are shown to the model, which is how the
    full and physics-only variants differ.

    Depth and the leaf constraint k_events are searched; column and row
    subsampling are fixed at the single values the per-system flatness
    measurements found the selection criterion insensitive to, and the learning
    rate is fixed by design with the round count left to early stopping (ceiling
    2000, patience 50). Calibration and the threshold see ONLY the inner
    rotation's held-out scores, never the held-out protein, so nothing about the
    test units enters the probability scale.
    """
    tr = np.flatnonzero(tr_mask)
    te = np.flatnonzero(te_mask)
    rot = (inner_groups[tr] if (INNER == "protein" and inner_groups is not None)
           else units[tr])
    r = fit_trees_fold(X[np.ix_(tr, cols_idx)], y[tr], rot, units[tr],
                       names, [2, 3, 4], [2.0, 4.0, 8.0], [0.7], [0.8],
                       lr=0.05, ceiling=2000, patience=50, monotone_map={},
                       seed=seed, n_threads=1)
    cal = Calibrator("sigmoid", 200).fit(r.inner_oof_scores, r.inner_oof_y)
    thr = f1_threshold(cal.predict(r.inner_oof_scores), r.inner_oof_y)
    raw = trees_predict(r.booster, X[np.ix_(te, cols_idx)], names)
    prob = cal.predict(raw)
    return raw, prob, thr, {"chosen": r.chosen, "M": r.M}


def _draw_timelines(score_df: pd.DataFrame, out_dir: Path) -> None:
    """One readiness timeline per out-of-training unit: the cross-trained
    monitor's P(ready) over time, true crossing completions marked in red."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sub = score_df[score_df["variant"] == "full"]
    for (exp, unit), g in sub.groupby(["experiment", "unit"]):
        g = g.sort_values("time_ps")
        t_ns = g["time_ps"].to_numpy() / 1000.0
        yv = g["y"].to_numpy().astype(bool)
        fig, ax = plt.subplots(figsize=(9, 3))
        ax.plot(t_ns, g["prob"].to_numpy(), lw=0.6,
                label="P(ready), monitor trained without this unit")
        ends = np.flatnonzero(yv & ~np.concatenate([yv[1:], [False]]))
        for e in ends:
            ax.axvline(t_ns[e], color="crimson", alpha=0.5, lw=0.9)
        n_ev = len(ends)
        ax.set_xlabel("Time, ns")
        ax.set_ylabel("P(ready)")
        ax.set_title(f"{unit} — {exp}; red lines: {n_ev} true crossing completions")
        ax.legend(fontsize=7, loc="upper right")
        safe = unit.replace(":", "_").replace("/", "-")
        fig.savefig(out_dir / f"timeline_{exp}_{safe}.png", dpi=150,
                    bbox_inches="tight")
        plt.close(fig)


def main():
    """Run both experiments at PCM2_GEN_TAU and write the output directory."""
    pool, schema = load_all()
    model_cols = [c.name for c in schema if c.to_model]
    y, valid = pool_labels(pool, TAU_PS)
    X = pool[model_cols].to_numpy(np.float32)
    units = pool["unit"].to_numpy()
    conducting = np.array([is_conducting(s, c) for s, c in
                           zip(pool["system"], pool["condition"])])
    protein = pool["system"].map(PROTEIN_OF).to_numpy()

    # System-label screen on the conducting pool only. A column whose value
    # ranges do not
    # even overlap between two systems can be thresholded to recover which
    # system a frame
    # came from, and a pooled model is free to use it that way: it would then
    # reproduce
    # each system's base rate instead of measuring physics. The screen runs on the
    # training pool alone, since a range taken from the held-out protein would be
    # information the model is not allowed to have.
    sys_label_cols = []
    cpool = pool[conducting]
    for c in model_cols:
        # Disjoint iff the largest per-system minimum lies above the smallest
        # per-system
        # maximum. Columns measured on no frame of a system are skipped there
        # rather than
        # counted as an empty range.
        ranges = [(g[c].min(), g[c].max()) for _s, g in cpool.groupby("system")
                  if g[c].notna().any()]
        if len(ranges) >= 2 and max(r[0] for r in ranges) > min(r[1] for r in ranges):
            sys_label_cols.append(c)
    physics_cols = [c for c in model_cols if c not in sys_label_cols]
    variants = {"full": model_cols, "physics_only": physics_cols}
    print(f"system-label columns dropped in physics_only ({len(sys_label_cols)}): "
          f"{sys_label_cols}")

    results = {"tau_ps": TAU_PS, "lopo": {}, "negative_control": {},
               "system_label_columns": sys_label_cols}
    score_rows = []
    cond_arr = pool["condition"].to_numpy()
    time_arr = pool["time_ps"].to_numpy()

    # ---- Experiment 1: leave-one-protein-out over conducting channels.
    helds = ["gramicidin", "mthk", "kcsa_family"] + (["cx43"] if POOL_ADD_CX43 else [])
    for held in helds:
        tr_mask = conducting & valid & (protein != held)
        results["lopo"][held] = {}
        for vname, cols in variants.items():
            idx = [model_cols.index(c) for c in cols]
            te_units = sorted(set(units[conducting & (protein == held)]))
            te_mask = np.isin(units, te_units)
            raw, prob, thr, hp = fit_and_score(X, y, units, tr_mask, te_mask,
                                               idx, cols, seed=20260815,
                                               inner_groups=protein)
            per_unit = {}
            te = np.flatnonzero(te_mask)
            for k, i in enumerate(te):
                score_rows.append((f"lopo_{held}", vname, units[i], cond_arr[i],
                                   float(time_arr[i]), float(prob[k]), int(y[i]),
                                   bool(valid[i])))
            # Scored WITHIN each held-out trajectory against that trajectory's
            # own base
            # rate. Pooling the held protein's units into one average precision
            # would mix
            # base rates that differ several-fold, and the pooled number would
            # then move
            # with the mixture rather than with the ranking; a within-unit ratio
            # is also
            # insensitive to a constant score offset between units.
            for u in te_units:
                m = units[te] == u
                e = ap_with_base(y[te][m & valid[te]], raw[m & valid[te]])
                per_unit[u] = (round(e["ratio_to_chance"], 2) if e.get("defined")
                               else None)
            ratios = [v for v in per_unit.values() if v is not None]
            results["lopo"][held][vname] = {
                "per_unit_ratio": per_unit,
                "median_ratio": float(np.median(ratios)) if ratios else None,
                "hyper": hp}
            print(f"LOPO held={held:12s} [{vname}]: median x"
                  f"{results['lopo'][held][vname]['median_ratio']}")

    # ---- Experiment 2: train on all conducting; apply to the mutant arms.
    # The threshold comes from the training pool's inner rotation, so the READY
    # fractions
    # below are read at the operating point the monitor would actually be
    # deployed at.
    # Choosing a threshold on the control arms themselves would define silence into
    # existence. Every non-conducting unit is a test unit; several carry no
    # positive frame
    # at all, and there the AP ratio is left undefined and only the silence is
    # quoted.
    tr_mask = conducting & valid
    test_units = sorted(set(units[~conducting]))
    for vname, cols in variants.items():
        idx = [model_cols.index(c) for c in cols]
        te_mask = ~conducting
        raw, prob, thr, hp = fit_and_score(X, y, units, tr_mask, te_mask,
                                           idx, cols, seed=20260816,
                                           inner_groups=protein)
        te = np.flatnonzero(te_mask)
        for k, i in enumerate(te):
            score_rows.append(("negative_control", vname, units[i], cond_arr[i],
                               float(time_arr[i]), float(prob[k]), int(y[i]),
                               bool(valid[i])))
        # Reference READY-rate on the conducting E71A frames (in-training;
        # for scale only).
        entry = {"threshold": thr, "per_unit": {}}
        for u in test_units:
            m = units[te] == u
            pv = prob[m]
            e = ap_with_base(y[te][m & valid[te]], raw[m & valid[te]])
            entry["per_unit"][u] = {
                "n_frames": int(m.sum()),
                "base_rate": float(np.mean(y[te][m & valid[te]])) if (m & valid[te]).any() else 0.0,
                "ready_fraction": float(np.nanmean(pv >= thr)),
                "mean_prob": float(np.nanmean(pv)),
                "ap_ratio": round(e["ratio_to_chance"], 2) if e.get("defined") else None,
            }
        results["negative_control"][vname] = entry
        print(f"negative-control [{vname}]: " + "; ".join(
            f"{u.split(':')[1]}({'Na' if 'na' in u.split(':')[0] else 'K'}):"
            f" READY {v['ready_fraction']:.2f}"
            for u, v in entry["per_unit"].items()))

    out_name = (f"monitor-generalisation-tau{int(TAU_PS)}"
                + ("-cx43pool" if POOL_ADD_CX43 else "")
                + ("-pinner" if INNER == "protein" else ""))
    final = PROJECT_ROOT / "runs" / out_name
    # Same discipline as the pipeline steps: assemble in a hidden temporary
    # directory and
    # rename over the old one at the end, so an interrupted run leaves the previous
    # results intact instead of a half-written mixture of two runs.
    tmp = PROJECT_ROOT / "runs" / f".{out_name}.tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    (tmp / "results.json").write_text(json.dumps(results, indent=1,
                                                 ensure_ascii=False, default=str))
    # Per-frame scores of every out-of-training unit -> parquet + readiness
    # timelines with the unit's true crossings marked ("prediction vs reality").
    score_df = pd.DataFrame(score_rows,
                            columns=["experiment", "variant", "unit", "condition",
                                     "time_ps", "prob", "y", "valid"])
    score_df.to_parquet(tmp / "scores.parquet")
    _draw_timelines(score_df, tmp)
    pool_desc = ("gramicidin + MthK + KcsA-E71A + Cx43" if POOL_ADD_CX43
                 else "gramicidin + MthK + KcsA-E71A")
    lines = [f"# Conducting-channels generalisation (tau = {int(TAU_PS)} ps)", "",
             f"Training pool: {pool_desc} (conducting channels only).",
             "", "## Leave-one-protein-out (median AP over the held protein's own chance)",
             "", "| held-out protein | full columns | physics-only |", "|---|---|---|"]
    for held, d in results["lopo"].items():
        lines.append(f"| {held} | x{d['full']['median_ratio']} | "
                     f"x{d['physics_only']['median_ratio']} |")
    lines += ["", "## Negative control: monitor applied to unseen mutant arms", "",
              "| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |",
              "|---|---|---|---|---|---|"]
    for u, v in results["negative_control"]["full"]["per_unit"].items():
        lines.append(f"| {u} | {v['n_frames']} | {v['base_rate']:.4f} | "
                     f"{v['ready_fraction']:.3f} | {v['mean_prob']:.3f} | "
                     f"{v['ap_ratio'] if v['ap_ratio'] is not None else 'undefined'} |")
    lines += ["", f"Columns dropped in physics-only: {len(sys_label_cols)} "
              f"(listed in results.json)."]
    (tmp / "GENERALISATION.md").write_text("\n".join(lines))
    write_provenance(tmp, config.load(SYSTEM_CONFIGS[0]), "transfer",
                     extra={"experiment": "conducting-channels generalisation"})
    if final.exists():
        shutil.rmtree(final)
    tmp.rename(final)
    print(f"\nwritten: {final}/GENERALISATION.md")


if __name__ == "__main__":
    main()
