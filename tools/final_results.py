"""Final thesis deliverable: runs/FINAL — the three-block experiment summary.

Block A: the four proteins, each on its own terms (passport table + a 2x2
         panel of the six-way method comparison).
Block B: train-on-three / zero-shot-on-each (LOPO over four proteins):
         the central table, a four-panel prediction-vs-reality figure with
         each protein at its own primary horizon, and a per-unit median plot.
Block C: the final three-protein monitor applied to the non-conducting
         mutant arms (silence table + figure), plus the measured
         three-vs-four pool comparison behind the canonical-pool choice.

Blocks D, E and F were added later and carry the method-by-protein matrices, the
per-model training passport and the standard ML views; they are produced by the
same command.

    python tools/final_results.py          # writes runs/FINAL/

Everything is read from canonical run artifacts; nothing is refitted here. What
must already exist: the five per-system pipeline runs, the cross-protein runs at
all three horizons in both pool variants (monitor-generalisation-tau{1000,2000,
4000} and their -cx43pool siblings), and runs/FINAL/cx43_arms.json from
tools/cx43_arms.py — the Cx43 comparison arms are skipped, not estimated, when
that file is absent. No trajectory is opened, so this runs on the published tree
with no trajectory data present.
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pcm2 import config  # noqa: E402
from pcm2.runtime import run_dir  # noqa: E402
from pcm2.viz import COMPARE_ARMS, _arm_frames, _sigma_frame_scores  # noqa: E402
from sklearn.metrics import (average_precision_score,  # noqa: E402
                             balanced_accuracy_score, brier_score_loss,
                             f1_score, matthews_corrcoef,
                             precision_recall_curve, roc_auc_score, roc_curve)


# One compact style for every figure: small consistent type set via rcParams
# (no per-call font sizes), light grid, no top/right spines.
plt.rcParams.update({
    "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6,
    "figure.titlesize": 10, "axes.spines.top": False,
    "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.18,
    "grid.linewidth": 0.4, "legend.framealpha": 0.85,
    "axes.titlelocation": "left",
})

OUT = ROOT / "runs" / "FINAL"

# The index this script leaves in runs/FINAL. It is written from here rather than
# maintained by hand, so that re-running the command cannot silently replace a
# corrected index with a shorter one.
FINAL_README = """# FINAL: the three-block thesis experiment

Block A — the four proteins on their own terms: `tableA_proteins`,
`figA_methods_per_protein.png` (six-way method comparison per protein).

Block B — train on three, zero-shot on the fourth (all four rotations):
`tableB_lopo4`, `figB1_zeroshot_timelines.png` (prediction vs reality at each
protein's own horizon), `figB2_lopo_medians.png` (per-unit spread).

Block C — the final three-protein monitor on the non-conducting mutant arms:
`tableC_negative_controls`, `figC_controls.png` (silence + the one real leak),
and `tableC2_pool_choice` (measured three-vs-four pool comparison behind the
canonical-pool decision).

Blocks D, E and F were added later and come from the same command. D is the
method-by-protein matrix: `figD_method_matrix.png` with each protein at its own
horizon, `figD_method_matrix_tau{1000,2000,4000}.png` at the three common
horizons, and `tableD_method_matrix`. E is `tableE_ml_quality`, the training
passport of each in-system model. F is `figF_training_quality.png`, the same
models in the standard ML views (ROC with AUC, precision-recall with PR-AUC,
reliability with the Brier score).

All numbers are read from canonical run artifacts (104-column canon); sources:
runs/monitor-generalisation-tau*{,-cx43pool}, per-system runs/<id>/train and
events. Every `table*` above is written twice, as `.csv` and as `.md`, from the
same rows.

## What writes what

| artifact | written by |
|---|---|
| `tableA_proteins`, `figA_methods_per_protein.png` | `tools/final_results.py`, block A |
| `tableB_lopo4`, `figB1_zeroshot_timelines.png`, `figB2_lopo_medians.png` | `tools/final_results.py`, block B |
| `tableC_negative_controls`, `figC_controls.png`, `tableC2_pool_choice` | `tools/final_results.py`, block C |
| `tableD_method_matrix`, `figD_method_matrix.png`, `figD_method_matrix_tau{1000,2000,4000}.png` | `tools/final_results.py`, block D |
| `tableE_ml_quality` | `tools/final_results.py`, block E |
| `figF_training_quality.png` | `tools/final_results.py`, block F |
| `cx43_arms.json` | `tools/cx43_arms.py`, which has to run first: `final_results.py` reads it and leaves the Cx43 group of the matrices empty when it is missing |
| `figG_shap_example.png` | `tools/shap_example_figure.py` |
| `figG2_shap_mthk.png` | the same script with the system, condition, replica, time window and output path passed to its `main()`; no command line in the repository reproduces it |
| `figure_2_1_workflow.png` | `tools/workflow_figure.py`, whose panel text is literal and does not follow the pipeline by itself |
| `figH_methods_showcase.png` | exported by hand from `notebooks/methods_showcase.ipynb` |
| `tableT5_one_moment.txt` | exported by hand from the same notebook |
"""

SYSTEMS = ["gramicidin", "mthk", "kcsa", "cx43"]
# Held-protein key in the LOPO-4 results -> (system id, showcase unit, own tau)
# The three K-channel proteins are shown at 2000 ps, the horizon common to all
# of them —
# for KcsA that is not the 1000 ps its own config declares primary — and Cx43 at
# 4000 ps,
# matched to its measured median transit through a far wider pore. Every block
# reads this
# same mapping, so a number carried from one table to another refers to the same
# horizon.
HELD = {
    "gramicidin": ("gramicidin", "gramicidin:330K_1Msalt/04", 2000),
    "mthk": ("mthk", "mthk:200MV/2", 2000),
    "kcsa_family": ("kcsa", "kcsa:E71A/1", 2000),
    "cx43": ("cx43", "cx43:200mV/R2", 4000),
}
TAUS = (1000, 2000, 4000)
PRETTY = {"gramicidin": "Gramicidin A", "mthk": "MthK",
          "kcsa_family": "KcsA-E71A (family)", "kcsa": "KcsA-E71A", "cx43": "Cx43"}


def pool4(tau: int) -> dict:
    """results.json of the ablation run whose training pool includes Cx43."""
    return json.loads((ROOT / "runs" / f"monitor-generalisation-tau{tau}-cx43pool"
                       / "results.json").read_text())


def pool3(tau: int) -> dict:
    """results.json of the canonical three-conducting-protein run at this horizon."""
    return json.loads((ROOT / "runs" / f"monitor-generalisation-tau{tau}"
                       / "results.json").read_text())


def scores4(tau: int) -> pd.DataFrame:
    """Per-frame scores of the four-protein-pool ablation run."""
    return pd.read_parquet(ROOT / "runs" /
                           f"monitor-generalisation-tau{tau}-cx43pool" / "scores.parquet")


def scores3(tau: int) -> pd.DataFrame:
    """Per-frame scores of the canonical three-protein-pool run."""
    return pd.read_parquet(ROOT / "runs" /
                           f"monitor-generalisation-tau{tau}" / "scores.parquet")


# Cx43 is excluded from every pool; basis in tableC2. Pool members are scored
# by monitors trained on the other TWO conducting proteins; cx43 is scored
# zero-shot by the full three-protein monitor (the negative-control fit of the
# canonical runs).
def held_numbers(held: str, tau: int, variant: str):
    """(median AP ratio over the held protein's units, per-unit ratios).

    Both come from the canonical three-protein run. For the three pool members
    that is its leave-one-protein-out block; for Cx43, which is in no pool, it is
    the negative-control block of the same run, where the full three-protein
    monitor is applied to units it never saw. Median is None when no unit of the
    held protein carries both classes.
    """
    r = pool3(tau)
    if held == "cx43":
        pu = {u: d["ap_ratio"] for u, d in
              r["negative_control"][variant]["per_unit"].items()
              if u.startswith("cx43:")}
        ratios = [v for v in pu.values() if v]
        return (float(np.median(ratios)) if ratios else None), pu
    d = r["lopo"][held][variant]
    return d.get("median_ratio"), d.get("per_unit_ratio", {})


def held_scores(held: str, tau: int, variant: str = "full") -> pd.DataFrame:
    """Per-frame scores of the held protein's units, from the same run as
    held_numbers, so a figure and the table beside it cannot disagree."""
    sc = scores3(tau)
    if held == "cx43":
        return sc[(sc["experiment"] == "negative_control")
                  & (sc["variant"] == variant)
                  & sc["unit"].str.startswith("cx43:")]
    return sc[(sc["experiment"] == f"lopo_{held}") & (sc["variant"] == variant)]


# ---------------------------------------------------------------- block A
def block_a():
    """Block A: the protein passport table plus the per-protein method panels.

    One row per protein (Cx43 twice, see below) with the architecture, the frame
    count, the horizon it is read at, the base rate there, the crossing count and
    median transit, the in-system model's AP over chance and the event-level
    alarm reading. Writes tableA_proteins.{csv,md} and
    figA_methods_per_protein.png.
    """
    rows = []
    for s in SYSTEMS:
        cfg = config.load(f"configs/{s}.yaml")
        root = run_dir(cfg)
        evt = pd.read_parquet(root / "events" / "events.parquet")
        if s == "kcsa":  # the passport row is the conducting E71A unit
            evt = evt[evt["condition"] == "E71A"]
        # Transit time only over events whose entrance was actually seen: an ion
        # already
        # inside the pore on the first frame has no entrance, and exit minus the
        # start of
        # the trajectory would be a property of when recording began. The median
        # is the
        # quantity the declared horizons are cross-checked against.
        obs = evt[evt["entrance_observed"]]
        transit = (float((obs["t_exit_ps"] - obs["t_entry_ps"]).median())
                   if len(obs) else float("nan"))
        alarm = json.loads((root / "train" / "alarm.json").read_text())
        he = alarm.get("headline_event_level", {})
        evaluation = json.loads((root / "train" / "evaluation.json").read_text())
        # Showcase horizon: same per-protein tau as blocks B/D (KcsA at 2000,
        # not its config primary 1000), so every block reads consistently.
        held_key = next(k for k, v in HELD.items() if v[0] == s)
        head = cfg["model.head"]
        lab = pd.read_parquet(root / "labels" / "labels.parquet")
        warned = (f"{he['warned_frac_mean']:.2f} vs "
                  f"{he['random_null_warned_frac_mean']:.2f} random"
                  if he and np.isfinite(he.get("warned_frac_mean", np.nan))
                  else "n/a")
        # Cx43 appears twice: at the common horizon for comparability, and at
        # the horizon matched to its measured transit.
        variants = ([(2000, "Cx43 (common τ)"), (4000, "Cx43 (transit-matched τ)")]
                    if s == "cx43" else [(HELD[held_key][2], PRETTY.get(s, s))])
        for primary, prot_label in variants:
            pooled = (evaluation["per_arm"].get(head, {})
                      .get(str(float(primary)), {}).get("pooled", {}))
            ycol = f"y_tau{primary}"
            base = float(lab[ycol].mean())
            in_sys = (f"x{pooled['ratio_to_chance']:.2f}" if pooled.get("defined")
                      else "refused (degenerate folds)")
            rows.append({
                "Protein": prot_label,
                "Architecture": {"gramicidin": "beta-helix dimer, single-file",
                                 "mthk": "tetramer, RCK-gated K+",
                                 "kcsa": "tetramer, filter mutants (K+/Na+ arms)",
                                 "cx43": "12-chain gap junction, 2 membranes"}[s],
                "Frames analysed": len(lab),
                "Horizon τ (ps)": primary,
                "Base rate at τ": round(base, 4),
                "Ion crossings": int(len(evt)),
                "Median transit (ps)": round(transit) if np.isfinite(transit) else "n/a",
                "In-system model, AP/chance": in_sys,
                "Alarm: crossings warned vs random": warned,
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tableA_proteins.csv", index=False)
    (OUT / "tableA_proteins.md").write_text(df.to_markdown(index=False))

    # The method-comparison panels cover the three proteins with an in-system
    # model matrix; cx43 has none (every fold refused), so its method
    # comparison appears in blocks B/D as the zero-shot construction.
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))
    for ax, s in zip(axes.ravel(), [x for x in SYSTEMS if x != "cx43"]):
        cfg = config.load(f"configs/{s}.yaml")
        root = run_dir(cfg)
        oof = pd.read_parquet(root / "train" / "oof.parquet")
        held_key = next(k for k, v in HELD.items() if v[0] == s)
        primary = float(HELD[held_key][2])
        sigma = _sigma_frame_scores(cfg)
        base = None
        for arm, label, color in COMPARE_ARMS:
            sub = _arm_frames(oof, sigma, arm, primary)
            if sub is None:
                continue
            m = np.isfinite(sub["score_raw"].to_numpy())
            yv, sv = sub["y"].to_numpy()[m], sub["score_raw"].to_numpy()[m]
            if len(yv) == 0 or yv.sum() in (0, len(yv)):
                continue
            prec, rec, thr = precision_recall_curve(yv, sv)
            base = float(np.mean(yv))
            ap = float(average_precision_score(yv, sv))
            ax.step(rec[:-1], prec[:-1], where="post", color=color, lw=1.1,
                    label=f"{label} — PR-AUC {ap:.3f}")
        if base is not None:
            ax.axhline(base, ls="--", c="gray", lw=0.9, label="base rate")
            ax.set_title(f"{PRETTY.get(s, s)} — τ={int(primary)} ps")
        else:
            # In-system training refused (event-starved): show the zero-shot
            # monitor instead, trained on the other three proteins and ranked
            # on this protein's own frames, plus the published rule, which
            # needs no fitting.
            held = next(k for k, v in HELD.items() if v[0] == s)
            own = HELD[held][2]
            sc = scores4(own)
            g = sc[(sc["experiment"] == f"lopo_{held}")
                   & (sc["variant"] == "full") & sc["valid"]]
            yv, pv = g["y"].to_numpy(int), g["prob"].to_numpy(float)
            if yv.sum():
                prec, rec, _ = precision_recall_curve(yv, pv)
                ax.step(rec[:-1], prec[:-1], where="post", color="#d62728",
                        lw=1.1, label="GBT zero-shot (trained on other three)")
            if sigma is not None:
                g = g.assign(replica=g["unit"].str.split("/").str[-1])
                m2 = g[["condition", "replica", "time_ps"]].merge(
                    sigma, on=["condition", "replica", "time_ps"])
                sv = m2["score_raw"].to_numpy(float)
                if len(sv) == len(yv) and 0 < yv.sum() < len(yv):
                    prec, rec, _ = precision_recall_curve(yv, sv)
                    ax.step(rec[:-1], prec[:-1], where="post", color="#1f77b4",
                            lw=1.1, label="Rao 2019 Σd (published rule)")
            ax.axhline(float(np.mean(yv)), ls="--", c="gray", lw=0.9,
                       label="base rate")
            ax.set_title(f"{PRETTY.get(s, s)} — zero-shot, τ = {own} ps")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.legend(loc="upper right")
    fig.suptitle("Per-frame precision–recall by method, pooled OOF")
    fig.tight_layout()
    fig.savefig(OUT / "figA_methods_per_protein.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- block B
def block_b():
    """Block B: train on three proteins, score the fourth it never saw.

    The table gives the median held-out AP over chance for each rotation at each
    of the three horizons, in the full and physics-only column variants side by
    side; the difference between the two columns is what bounds the contribution
    of the system-identifying descriptors. Writes tableB_lopo4.{csv,md},
    figB1_zeroshot_timelines.png and figB2_lopo_medians.png.
    """
    rows = []
    for held, (sysid, _u, own) in HELD.items():
        row = {"Held-out protein": PRETTY[held], "Own horizon τ (ps)": own}
        for tau in TAUS:
            for v in ("full", "physics_only"):
                mr, _pu = held_numbers(held, tau, v)
                row[("τ=%d ps, full set" % tau) if v == "full" else ("τ=%d ps, physics-only" % tau)] = f"x{mr:.2f}" if mr else "undefined"
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tableB_lopo4.csv", index=False)
    (OUT / "tableB_lopo4.md").write_text(df.to_markdown(index=False))

    # Four-panel prediction-vs-reality, each held protein at its own horizon.
    # The heavy line is a backward rolling mean of the per-frame probability
    # over ~1 ns, the alarm's smoothing scale.
    fig, axes = plt.subplots(4, 1, figsize=(11, 12), sharex=False)
    for ax, (held, (sysid, unit, own)) in zip(axes, HELD.items()):
        g = held_scores(held, own)
        g = g[g["unit"] == unit].sort_values("time_ps")
        t_ns = g["time_ps"].to_numpy() / 1000.0
        p = g["prob"].to_numpy()
        yv = g["y"].to_numpy().astype(bool)
        # The window is in picoseconds, so it must be converted to frames
        # through this
        # unit's own frame spacing: the systems run at strides of 100 to 250 ps and a
        # fixed frame count would mean a different physical window on each panel.
        dt = float(np.median(np.diff(g["time_ps"]))) if len(g) > 1 else 50.0
        w = max(3, int(round(1000.0 / dt)))
        smooth = pd.Series(p).rolling(w, min_periods=1).mean().to_numpy()
        ax.plot(t_ns, p, lw=0.4, color="#9db8d6", label="P(ready), per frame")
        ax.plot(t_ns, smooth, lw=1.3, color="#1f5fa8",
                label="P(ready), 1-ns rolling mean")
        # Last frame of each run of positive labels. Since y(t)=1 means a
        # completion falls
        # in (t, t+tau], the end of a run is the last frame from which the
        # crossing is
        # still ahead, i.e. one frame before it completes. Reading the crossings
        # off the
        # label column this way keeps the markers on exactly the events the score was
        # judged against, which re-reading the events table would not guarantee.
        ends = np.flatnonzero(yv & ~np.concatenate([yv[1:], [False]]))
        for e in ends:
            ax.axvline(t_ns[e], color="crimson", alpha=0.5, lw=1.0)
        med, _ = held_numbers(held, own, "full")
        src = ("trained on the three conducting proteins" if held == "cx43"
               else "trained on the other two conducting proteins")
        ax.set_title(f"{PRETTY[held]} · {unit.split(':')[1]} · τ = {own} ps")
        ax.set_ylabel("P(ready)")
        ax.legend(loc="upper right")
    axes[-1].set_xlabel("Time, ns")
    fig.suptitle("Held-out P(ready) time series and true crossing completions")
    fig.tight_layout()
    fig.savefig(OUT / "figB1_zeroshot_timelines.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    # Median + per-unit dots, full vs physics-only, at each protein's own tau.
    # The ratio ceiling (1/base rate) differs by orders of magnitude between
    # proteins, so it is drawn per group and each protein's ratios are read
    # against its own ceiling.
    fig, ax = plt.subplots(figsize=(9, 5.5))
    xt, xl = [], []
    for i, (held, (sysid, _u, own)) in enumerate(HELD.items()):
        sub = held_scores(held, own)
        sub = sub[sub["valid"]]
        bases = sub.groupby("unit")["y"].mean()
        bases = bases[bases > 0]
        if len(bases):
            # A perfect ranker reaches AP = 1, so the largest attainable ratio
            # to chance is
            # 1/base rate. On a unit with a base rate of 0.1% that ceiling is
            # 1000 and on
            # one at 50% it is 2, which is why the ratios of different proteins
            # must not be
            # compared to each other without their own ceiling drawn beside them.
            ceil = float(1.0 / bases.median())
            ax.hlines(ceil, i - 0.32, i + 0.32, color="#b3593b", lw=1.4, ls=":",
                      label="ceiling 1/base (median unit)" if i == 0 else None)
            ax.annotate(f"{ceil:.0f}", (i + 0.34, ceil),
                        color="#b3593b", va="center")
        for j, (v, color) in enumerate((("full", "#1f5fa8"), ("physics_only", "#2ca02c"))):
            x = i + (j - 0.5) * 0.3
            _med, pu = held_numbers(held, own, v)
            ratios = [r for r in pu.values() if r]
            if not ratios:
                continue
            ax.scatter([x] * len(ratios), ratios, s=22, alpha=0.55, color=color,
                       label=v if i == 0 else None)
            ax.hlines(float(np.median(ratios)), x - 0.12, x + 0.12, color=color, lw=2.5)
        xt.append(i)
        xl.append(f"{PRETTY[held]}\nτ={own} ps")
    ax.axhline(1.0, ls="--", c="gray", lw=1, label="chance")
    ax.set_yscale("log")
    ax.set_yticks([1, 2, 3, 5, 10, 15], ["1", "2", "3", "5", "10", "15"])
    ax.set_xticks(xt, xl)
    ax.set_ylabel("AP / held unit's own chance (log scale)")
    ax.set_title("Held-out AP / base rate per unit, with medians and ceilings")
    ax.legend()
    fig.savefig(OUT / "figB2_lopo_medians.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- block C
def block_c():
    """Block C: the three-protein monitor on the non-conducting mutant arms.

    The claim being tested is silence: a monitor that has learned physics rather
    than the appearance of a K channel must read arms that do not conduct as
    almost never ready. The one exception is the Na+ G77A-E71A arm, whose
    residual leak makes it the single control where ranking is definable at all.
    Writes tableC_negative_controls.{csv,md}, figC_controls.png and
    tableC2_pool_choice.{csv,md}, the measured basis for keeping Cx43 out of the
    training pool.
    """
    rows = []
    arms = None
    for tau in TAUS:
        nc = pool3(tau)["negative_control"]["full"]
        arms = [u for u in nc["per_unit"] if not u.startswith("cx43:")]
        for u in arms:
            d = nc["per_unit"][u]
            rows.append({"Mutant arm": u, "Horizon τ (ps)": tau, "Frames": d["n_frames"],
                         "True base rate": round(d["base_rate"], 4),
                         "READY fraction": round(d["ready_fraction"], 4),
                         "Mean P(ready)": round(d["mean_prob"], 3),
                         "AP / chance": d["ap_ratio"] if d["ap_ratio"] else "undefined"})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tableC_negative_controls.csv", index=False)
    (OUT / "tableC_negative_controls.md").write_text(df.to_markdown(index=False))

    fig = plt.figure(figsize=(11, 7))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.1, 1])
    ax = fig.add_subplot(gs[0])
    width = 0.25
    for k, tau in enumerate(TAUS):
        nc = pool3(tau)["negative_control"]["full"]["per_unit"]
        vals = [nc[a]["mean_prob"] for a in arms]
        xs = np.arange(len(arms)) + (k - 1) * width
        ax.bar(xs, vals, width, label=f"mean P(ready), τ={tau} ps")
        # Chance reference: a perfectly calibrated monitor's mean probability
        # equals the arm's true base rate at this horizon. So a bar standing
        # above its own
        # black line is over-confident on an arm that does not conduct, and that
        # gap, not
        # the bar height, is the reading of this panel.
        FLOOR = 2e-3  # log-axis floor; a zero base rate is drawn at the floor
        bases = [max(nc[a]["base_rate"], FLOOR) for a in arms]
        ax.hlines(bases, xs - width * 0.55, xs + width * 0.55,
                  color="#22262b", lw=1.8, zorder=3,
                  label="chance = true base rate" if k == 0 else None)
    ax.set_xticks(np.arange(len(arms)),
                  [a.replace("kcsa_na:", "Na: ").replace("kcsa:", "K: ")
                   for a in arms])
    ax.axhline(0.5, ls="--", c="gray", lw=0.8)
    ax.text(0.01, 0.55, "READY territory", color="gray",
            transform=ax.get_yaxis_transform())
    ax.set_yscale("log")
    ax.set_ylim(2e-3, 1.0)
    ax.set_ylabel("mean P(ready), log scale")
    ax.set_title("Non-conducting mutant arms: mean P(ready) by horizon")
    ax.legend()

    ax2 = fig.add_subplot(gs[1])
    sc = scores3(1000)
    g = sc[(sc["experiment"] == "negative_control") & (sc["variant"] == "full")
           & (sc["unit"] == "kcsa_na:G77AE71A/1")].sort_values("time_ps")
    t_ns = g["time_ps"].to_numpy() / 1000.0
    yv = g["y"].to_numpy().astype(bool)
    p = g["prob"].to_numpy()
    dt = float(np.median(np.diff(g["time_ps"]))) if len(g) > 1 else 50.0
    # 10-ns smoothing window, set by the 2-us length of this trajectory: at the
    # panel's
    # width a 1-ns mean is a solid band, and the question here is whether the
    # probability
    # stays low over microseconds, not what it does frame by frame.
    w = max(3, int(round(10000.0 / dt)))
    smooth = pd.Series(p).rolling(w, min_periods=1).mean().to_numpy()
    ax2.plot(t_ns, smooth, lw=1.0, color="#1f5fa8",
             label="P(ready), 10-ns rolling mean")
    ax2.legend(loc="upper right")
    ends = np.flatnonzero(yv & ~np.concatenate([yv[1:], [False]]))
    for e in ends:
        ax2.axvline(t_ns[e], color="crimson", alpha=0.5, lw=1.0)
    nc = pool3(1000)["negative_control"]["full"]["per_unit"]["kcsa_na:G77AE71A/1"]
    ax2.set_title(f"G77A\u2013E71A / Na+ · P(ready), 10-ns rolling mean · τ = 1000 ps")
    ax2.set_xlabel("Time, ns")
    ax2.set_ylabel("P(ready)")
    fig.tight_layout()
    fig.savefig(OUT / "figC_controls.png", dpi=160, bbox_inches="tight")
    plt.close(fig)

    # Pool choice: three-protein vs four-protein medians side by side.
    rows = []
    for held in ("gramicidin", "mthk", "kcsa_family"):
        row = {"Held-out protein": PRETTY[held]}
        for tau in TAUS:
            p3 = pool3(tau)["lopo"][held]["full"].get("median_ratio")
            p4 = pool4(tau)["lopo"][held]["full"].get("median_ratio")
            row[f"τ={tau} ps, 3-protein pool"] = f"x{p3:.2f}"
            row[f"τ={tau} ps, +Cx43 pool"] = f"x{p4:.2f}"
        rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tableC2_pool_choice.csv", index=False)
    (OUT / "tableC2_pool_choice.md").write_text(df.to_markdown(index=False))


# ---------------------------------------------------------------- block D
METHOD_LABEL = {
    "published_rao2019_sigma": "Rao 2019 Σd",
    "plsfma_coords": "PLS-FMA",
    "linear": "Linear",
    "trees": "GBT (in-system)",
    "clock": "Clock",
}
METHOD_COLOR = {
    "published_rao2019_sigma": "#1f77b4", "plsfma_coords": "#ff7f0e",
    "linear": "#2ca02c", "trees": "#d62728", "clock": "#9467bd",
    "zero_shot": "#7a4b2b",
}


def _method_ratios(sysid: str, tau: float) -> dict:
    """AP/chance of every comparison method on one system at one horizon,
    computed uniformly from the stored OOF scores (so any declared tau works,
    not only the primary)."""
    from pcm2.evaluate import ap_with_base
    cfg = config.load(f"configs/{sysid}.yaml")
    oof = pd.read_parquet(run_dir(cfg) / "train" / "oof.parquet")
    sigma = _sigma_frame_scores(cfg)
    out = {}
    for arm, _label, _c in COMPARE_ARMS:
        sub = _arm_frames(oof, sigma, arm, tau)
        if sub is None:
            continue
        s = sub["score_raw"].to_numpy(float)
        m = np.isfinite(s)
        y = sub["y"].to_numpy(int)
        if m.sum() and 0 < y[m].sum() < m.sum():
            e = ap_with_base(y[m], s[m])
            if e.get("defined"):
                out[arm] = e["ratio_to_chance"]
    return out


CX43_ARMS = OUT / "cx43_arms.json"


def _cx43_arms(tau: float) -> dict:
    """Comparison arms for Cx43, written by tools/cx43_arms.py.

    Both classes sit in one of its two assemblies, so the fold that fits has
    nothing to rank and no in-system model exists. The published rules still
    apply frame by frame, and the fitted arms are applied zero-shot from the
    conducting pool, which is how the monitor itself is scored here."""
    if not CX43_ARMS.exists():
        return {}
    by_tau = json.loads(CX43_ARMS.read_text())["by_tau"]
    row = by_tau.get(str(int(tau)), {})
    label_to_arm = {"Rao 2019 energy": "published_rao2019",
                    "Rao 2019 sigma": "published_rao2019_sigma",
                    "Clock (temporal null)": "clock",
                    "Linear control (structural)": "linear_control_structural"}
    return {label_to_arm[k]: v for k, v in row.items()
            if k in label_to_arm and v is not None}


def _draw_method_matrix(tag: str, tau_of: dict, fname: str, title: str):
    """One method-by-protein bar chart. The Cx43 group is zero-shot or unfitted
    throughout: it has no in-system model."""
    fig, ax = plt.subplots(figsize=(10, 5))
    methods = [a for a, _l, _c in COMPARE_ARMS] + ["zero_shot"]
    width = 0.13
    for k, mth in enumerate(methods):
        vals, xs = [], []
        for i, s in enumerate(SYSTEMS):
            tau = tau_of[s]
            if mth == "zero_shot":
                held = next(kk for kk, v in HELD.items() if v[0] == s)
                v, _pu = held_numbers(held, tau, "full")
            elif s == "cx43":
                v = _cx43_arms(float(tau)).get(mth)
            else:
                v = _method_ratios(s, float(tau)).get(mth)
            if v is not None and np.isfinite(v):
                vals.append(v)
                xs.append(i + (k - (len(methods) - 1) / 2) * width)
        label = ("GBT zero-shot (trained on the other conducting proteins)"
                 if mth == "zero_shot"
                 else dict((a, l) for a, l, _ in COMPARE_ARMS)[mth])
        color = (METHOD_COLOR["zero_shot"] if mth == "zero_shot"
                 else dict((a, c) for a, _l, c in COMPARE_ARMS)[mth])
        ax.bar(xs, vals, width * 0.92, color=color, label=label)
    ax.axhline(1.0, ls="--", c="gray", lw=1)
    ax.text(0.005, 1.03, "chance", color="gray",
            transform=ax.get_yaxis_transform())
    ax.set_xticks(range(len(SYSTEMS)),
                  [f"{PRETTY.get(s, s)}\nτ={tau_of[s]} ps" for s in SYSTEMS])
    ax.set_ylabel("AP / base rate (per-frame ranking)")
    ax.set_title(title)
    ax.legend(ncols=2)
    fig.savefig(OUT / fname, dpi=160, bbox_inches="tight")
    plt.close(fig)


def block_d():
    """Method-by-protein matrices: the showcase chart (each protein at its
    own horizon) plus one chart per common horizon, and the summary table."""
    show = {s: HELD[next(k for k, v in HELD.items() if v[0] == s)][2]
            for s in SYSTEMS}
    _draw_method_matrix("showcase", show, "figD_method_matrix.png",
                        "AP / base rate by method and protein, per-protein horizons")
    for tau in TAUS:
        _draw_method_matrix(str(tau), {s: tau for s in SYSTEMS},
                            f"figD_method_matrix_tau{tau}.png",
                            f"AP / base rate by method and protein, τ = {tau} ps")

    label_of = dict((a, l) for a, l, _ in COMPARE_ARMS)
    rows = []
    for a, l, _c in COMPARE_ARMS:
        row = {"Method": l}
        for s in SYSTEMS:
            if s == "cx43":
                row[PRETTY.get(s, s)] = "—"
            else:
                v = _method_ratios(s, float(show[s])).get(a)
                row[PRETTY.get(s, s)] = (f"x{v:.2f}" if v is not None else
                                         "refused/undefined")
        rows.append(row)
    zrow = {"Method": "GBT zero-shot (other conducting proteins)"}
    for s in SYSTEMS:
        held = next(k for k, v in HELD.items() if v[0] == s)
        v, _pu = held_numbers(held, show[s], "full")
        zrow[PRETTY.get(s, s)] = f"x{v:.2f}" if v else "undefined"
    rows.append(zrow)
    rows.append({"Method": "Base rate (chance)",
                 **{PRETTY.get(s, s): "x1.00" for s in SYSTEMS}})
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tableD_method_matrix.csv", index=False)
    (OUT / "tableD_method_matrix.md").write_text(df.to_markdown(index=False))


# ---------------------------------------------------------------- block F
def block_f():
    """Training quality of the in-system models in the standard ML views: ROC
    with AUC, PR with PR-AUC (the primary view under class imbalance, Saito &
    Rehmsmeier 2015), reliability diagram with the Brier score
    (Niculescu-Mizil & Caruana 2005). Out-of-fold, head arm, primary horizon,
    only the systems whose in-system training succeeded."""
    trainable = []
    for s in SYSTEMS:
        cfg = config.load(f"configs/{s}.yaml")
        oof = pd.read_parquet(run_dir(cfg) / "train" / "oof.parquet")
        primary = cfg["labels.primary_tau_ps"]
        sub = oof[(oof["arm"] == cfg["model.head"]) & (oof["tau"] == primary)
                  & oof["valid"]]
        sc = sub["score_raw"].to_numpy(float)
        m = np.isfinite(sc)
        yv = sub["y"].to_numpy(int)[m]
        if len(yv) and 0 < yv.sum() < len(yv):
            trainable.append((s, int(primary), yv, sc[m],
                              sub["prob_cal"].to_numpy(float)[m]))
    fig, axes = plt.subplots(len(trainable), 3,
                             figsize=(13, 4.2 * len(trainable)))
    axes = np.atleast_2d(axes)
    for row, (s, primary, yv, sc, pc) in zip(axes, trainable):
        fpr, tpr, _ = roc_curve(yv, sc)
        auc = roc_auc_score(yv, sc)
        row[0].plot(fpr, tpr, color="#d62728", lw=1.2,
                    label=f"ROC-AUC {auc:.3f}")
        row[0].plot([0, 1], [0, 1], "--", c="gray", lw=0.9, label="chance")
        row[0].set_xlabel("False positive rate")
        row[0].set_ylabel("True positive rate")
        row[0].set_title(f"{PRETTY.get(s, s)} — ROC (OOF, τ={primary} ps)")
        prec, rec, _ = precision_recall_curve(yv, sc)
        ap = average_precision_score(yv, sc)
        base = float(np.mean(yv))
        row[1].step(rec[:-1], prec[:-1], where="post", color="#d62728", lw=1.2,
                    label=f"PR-AUC {ap:.3f}")
        row[1].axhline(base, ls="--", c="gray", lw=0.9,
                       label=f"chance = base rate {base:.3f}")
        row[1].set_xlabel("Recall")
        row[1].set_ylabel("Precision")
        row[1].set_title("Precision–recall")
        mf = np.isfinite(pc)
        qs = np.quantile(pc[mf], np.linspace(0, 1, 11))
        qs[0] -= 1e-9
        xs, ys = [], []
        for a, b in zip(qs[:-1], qs[1:]):
            mm = (pc[mf] > a) & (pc[mf] <= b)
            if mm.sum():
                xs.append(float(np.mean(pc[mf][mm])))
                ys.append(float(np.mean(yv[mf][mm])))
        brier = brier_score_loss(yv[mf], pc[mf])
        row[2].plot(xs, ys, "o-", color="#d62728", ms=4, lw=1.1,
                    label=f"Brier {brier:.4f}")
        lim = [min(xs + ys) - 0.02, max(xs + ys) + 0.02]
        row[2].plot(lim, lim, "--", c="gray", lw=0.9, label="ideal")
        row[2].set_xlabel("Predicted P(ready)")
        row[2].set_ylabel("Observed positive fraction")
        row[2].set_title("Reliability, decile bins")
        for ax in row:
            ax.legend()
    fig.suptitle("Out-of-fold ROC, precision\u2013recall and calibration")
    fig.tight_layout()
    fig.savefig(OUT / "figF_training_quality.png", dpi=160, bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- block E
def block_e():
    """ML training-quality passport per in-system model: data volume, class
    balance, fold coverage vs refusals, chosen hyperparameters, boosting
    rounds, OOF ranking quality, probability calibration (Brier), and the
    acceptance-gate verdict."""
    rows = []
    for s in SYSTEMS:
        cfg = config.load(f"configs/{s}.yaml")
        root = run_dir(cfg)
        held_key = next(k for k, v in HELD.items() if v[0] == s)
        head = cfg["model.head"]
        folds = json.loads((root / "train" / "folds.json").read_text())
        variants = ([(2000.0, "Cx43 (common τ)"), (4000.0, "Cx43 (transit-matched τ)")]
                    if s == "cx43" else [(float(HELD[held_key][2]), PRETTY.get(s, s))])
        primary = variants[0][0]
        ev = json.loads((root / "train" / "evaluation.json").read_text())
        lab = pd.read_parquet(root / "labels" / "labels.parquet")
        ycol, vcol = f"y_tau{int(primary)}", f"valid_tau{int(primary)}"
        fitted = refused = 0
        chosen = []
        for fa in folds:
            if fa["tau"] != primary:
                continue
            a = fa.get("arms", {}).get(head, {})
            if "skipped" in a:
                refused += 1
            elif a:
                fitted += 1
                if "hyperparams" in a:
                    chosen.append(a["hyperparams"])
        pooled = ev["per_arm"].get(head, {}).get(str(primary), {}).get("pooled", {})
        entry = ev["per_arm"].get(head, {}).get(str(primary), {})
        acc_path = root / "report" / "acceptance.json"
        gates = "n/a"
        if acc_path.exists():
            # The gate document nests: some gates hold a verdict directly,
            # others hold one
            # per hyperparameter. Walking the whole tree for "ok" flags counts every
            # verdict actually recorded, so a gate added later is included
            # without this
            # summary having to know its name.
            def oks(node, out):
                if isinstance(node, dict):
                    for k, v in node.items():
                        if k == "ok" and isinstance(v, bool):
                            out.append(v)
                        else:
                            oks(v, out)
                elif isinstance(node, list):
                    for v in node:
                        oks(v, out)
                return out
            flags = oks(json.loads(acc_path.read_text()), [])
            gates = f"{sum(flags)}/{len(flags)} checks ok" if flags else "n/a"
        hp = "—"
        if chosen:
            c = chosen[0].get("chosen", chosen[0])
            hp = ", ".join(f"{k}={c[k]}" for k in ("depth", "k_events", "colsample")
                           if k in c)
            if chosen[0].get("M_rounds"):
                hp += f", M={chosen[0]['M_rounds']}"
        # Standard ML metrics on the OOF head scores and thresholded verdicts.
        oof = pd.read_parquet(root / "train" / "oof.parquet")
        sub = oof[(oof["arm"] == head) & (oof["tau"] == primary) & oof["valid"]]
        sc = sub["score_raw"].to_numpy(float)
        mfin = np.isfinite(sc)
        yv = sub["y"].to_numpy(int)[mfin]
        prauc = rocauc = mcc = f1 = bal = None
        if len(yv) and 0 < yv.sum() < len(yv):
            prauc = float(average_precision_score(yv, sc[mfin]))
            rocauc = float(roc_auc_score(yv, sc[mfin]))
            bench = pd.read_parquet(root / "benchmark" / "benchmark.parquet")
            if f"{head}_verdict" in bench.columns:
                bm = bench[bench["valid"]
                           & (bench[f"{head}_verdict"] != "N/A")]
                bt = bm["truth_crossing_within_tau"].to_numpy(int)
                bp = (bm[f"{head}_verdict"] == "READY").to_numpy(int)
                if len(bt) and 0 < bt.sum() < len(bt):
                    mcc = float(matthews_corrcoef(bt, bp))
                    f1 = float(f1_score(bt, bp))
                    bal = float(balanced_accuracy_score(bt, bp))
        fmt = lambda v, n=3: (round(v, n) if v is not None else "—")
        for primary, sys_label in variants:
            ycol, vcol = f"y_tau{int(primary)}", f"valid_tau{int(primary)}"
            rows.append({
                "System": sys_label,
                "Frames": len(lab),
                "Horizon τ (ps)": int(primary),
                "Positive frames": int(lab[ycol].sum()),
                "Positive fraction": round(float(lab[ycol][lab[vcol].astype(bool)].mean()), 4),
                "Folds fitted / refused": f"{fitted}/{refused}",
                "Chosen hyperparameters (first fold)": hp,
                "PR-AUC": fmt(prauc),
                "PR-AUC / chance": (f"x{pooled['ratio_to_chance']:.2f}"
                                    if pooled.get("defined") else "undefined"),
                "ROC-AUC": fmt(rocauc),
                "MCC": fmt(mcc),
                "F1": fmt(f1),
                "Balanced accuracy": fmt(bal),
                "Brier score": (round(entry["brier"], 4) if "brier" in entry else "—"),
                "Acceptance checks": gates,
            })
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "tableE_ml_quality.csv", index=False)
    (OUT / "tableE_ml_quality.md").write_text(df.to_markdown(index=False))


def main():
    """Write every block into runs/FINAL and list what was produced."""
    OUT.mkdir(exist_ok=True)
    block_a()
    block_b()
    block_c()
    block_d()
    block_e()
    block_f()
    (OUT / "README.md").write_text(FINAL_README)
    print(f"written: {OUT}")
    for p in sorted(OUT.iterdir()):
        print(" ", p.name)


if __name__ == "__main__":
    main()
