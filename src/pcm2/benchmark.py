"""Cross-method benchmark: every surveyed method judged on its own valid frames.

Output is one row per frame: the ground truth, and for every method the score,
the verdict in the method's own terms, and whether that verdict was correct,
alongside this system's answer (probability of readiness at the declared
horizon, the READY/NOT_READY verdict, and the diagnosed pore state).

Each method is scored against what it claims to predict. The published
Rao-2019 heuristic predicts local dewetting of the constriction rather than
conduction; its readiness proxy is "wetted => ready", so a weak proxy score is
a property of that transfer.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

from .baselines.published import (RAO_DEWETTING_THRESHOLD_KJMOL,
                                  RAO_SIGMA_D_THRESHOLD)
from .config import Config
from .evaluate import ap_with_base
from .runtime import StepLog, run_dir, step_output

# What each arm claims to predict, in its own terms.
METHOD_CLAIMS = {
    "trees": "this work (head model): P(crossing completes within tau), READY/NOT_READY "
             "verdict and mechanism axis (pore state) per frame",
    "linear": "this work (linear branch): same probabilistic claim, ridge-logistic on all "
              "descriptors",
    "plsfma_coords": "PLS-FMA (Krivobokova et al. 2012; the authors' g_fma PLS core, "
                     "Helland/Denham): collective mode maximally covarying with the "
                     "readiness label, regressed on superimposed C-alpha Cartesian "
                     "coordinates — the published method's machinery applied to a "
                     "binary functional quantity",
    "linear_control_structural": "fitted linear control on the published criterion's "
                                 "inputs (constriction radius + lining hydrophobicity)",
    "clock": "temporal null: event seriality only (time since last crossing, event rate, "
             "pore ion count); no structural input",
    "published_rao2019": "Rao et al. 2019 heuristic: local DEWETTING of the constriction "
                         "(frame-wise application); readiness proxy = wetted",
    "published_rao2019_sigma": "Rao et al. 2019 STRUCTURE classifier as published: "
                               "sum of contour distances (sigma_d) over all "
                               "pore-facing residue points past the 1RT line; "
                               "non-conductive when sigma_d > 0.55 — the rule the "
                               "paper's AUC 0.91 belongs to",
    "published_rao2019_win": "Rao et al. 2019 heuristic on window-averaged inputs "
                             "(closer to the authors' time-averaged definition)",
}


def _balanced_accuracy(y: np.ndarray, pred: np.ndarray) -> float:
    """Mean of sensitivity and specificity of a binary verdict; NaN if undefined.

    Plain verdict accuracy is uninformative at these base rates: on a system
    where one frame in a hundred is positive, answering NOT_READY everywhere
    scores 0.99. Averaging the two class-wise rates gives every method 0.5 for
    a constant answer, so the number can be compared across systems whose base
    rates differ by orders of magnitude. NaN when one class is missing or no
    frame carries a verdict, since the average is then not defined.
    """
    m = np.isfinite(pred.astype(float))
    y, pred = y[m], pred[m].astype(bool)
    if len(y) == 0 or y.sum() in (0, len(y)):
        return float("nan")
    tpr = float(np.mean(pred[y == 1]))
    tnr = float(np.mean(~pred[y == 0]))
    return 0.5 * (tpr + tnr)


def run_step(cfg: Config) -> None:
    """Write benchmark/benchmark.{parquet,csv}, BENCHMARK.md and summary.json.

    One row per frame at the primary horizon, with the ground truth once and,
    per method, its score, its verdict in its own vocabulary (READY/NOT_READY,
    WETTED/DEWETTED, CONDUCTIVE/NON_CONDUCTIVE) and whether that verdict was
    right. Nothing is fitted: the fitted arms are read from the train step's
    out-of-fold table, so no method ever sees a frame it was trained on, and the
    unfitted published rules are evaluated from the descriptor columns.
    """
    root = run_dir(cfg)
    primary = cfg["labels.primary_tau_ps"]
    head = cfg["model.head"]
    oof = pd.read_parquet(root / "train" / "oof.parquet")
    oof = oof[oof["tau"] == primary]
    # Per-frame answers exist only when the head model fitted at least one
    # fold; on a fully refused system the benchmark still compares the other
    # methods, just without the pcm_* columns.
    answers_path = root / "train" / "answers.parquet"
    answers = pd.read_parquet(answers_path) if answers_path.exists() else None
    folds_doc = json.loads((root / "train" / "folds.json").read_text())
    evaluation = json.loads((root / "train" / "evaluation.json").read_text())
    # Column list is explicit for speed, but the sigma column only exists in
    # tables built after the Rao structure classifier was added — include it
    # when the table carries it, or the sigma arm below silently never runs.
    feat_path = root / "features" / "features.parquet"
    available = pq.read_schema(feat_path).names
    feat_cols = ["condition", "replica", "time_ps",
                 "bl_rao_E_kJmol", "bl_rao_E_win_kJmol"]
    feat_cols += [c for c in ("bl_rao_sigma_d", "bl_rao_flag_n")
                  if c in available]
    features = pd.read_parquet(feat_path, columns=feat_cols)
    thresholds: dict[tuple[str, str], float] = {}
    for fa in folds_doc:
        if fa["tau"] != primary:
            continue
        for arm, a in fa.get("arms", {}).items():
            if "threshold" in a:
                thresholds[(arm, fa["fold"])] = a["threshold"]

    key = ["condition", "replica", "time_ps"]
    with step_output(cfg, "benchmark") as out:
        log = StepLog(out)
        # The frame base comes from the widest-coverage arm: on event-starved
        # systems the head may have refused folds (or every fold), and a
        # head-only base would silently shrink or empty the benchmark.
        prim = oof[oof["tau"] == primary] if "tau" in oof.columns else oof
        cov = prim.groupby("arm").size()
        base_arm = head if cov.get(head, 0) == cov.max() else cov.idxmax()
        if base_arm != head:
            log.say(f"frame base from '{base_arm}' (head '{head}' covers "
                    f"{cov.get(head, 0)} of {cov.max()} frames)")
        base = (oof[oof["arm"] == base_arm][key + ["y", "valid", "in_transit", "fold"]]
                .drop_duplicates(key).reset_index(drop=True))
        table = base.rename(columns={"y": "truth_crossing_within_tau"})
        table.insert(3, "tau_ps", primary)
        truth = table["truth_crossing_within_tau"].to_numpy()
        valid = table["valid"].to_numpy(bool)

        summary_rows = []
        arms = [a for a in oof["arm"].unique()]
        for arm in sorted(arms, key=lambda a: (a != head, a)):
            sub = oof[oof["arm"] == arm].drop_duplicates(key)
            merged = table[key + ["fold"]].merge(sub, on=key, how="left",
                                                validate="one_to_one",
                                                suffixes=("", "_oof"))
            score = merged["score_raw"].to_numpy(float)
            prob = merged["prob_cal"].to_numpy(float)
            if arm.startswith("published"):
                col = ("bl_rao_E_kJmol" if arm == "published_rao2019"
                       else "bl_rao_E_win_kJmol")
                e = table[key].merge(features, on=key, validate="one_to_one")[col]
                # Above the ~1 RT contour the constriction is predicted to dewet. The
                # rule's own output is a wetting call, not a readiness call; the
                # transfer
                # to this task is the declared proxy "wetted, therefore able to
                # conduct",
                # and a weak number here is a statement about that proxy rather than
                # about the authors' criterion.
                dewetted = e.to_numpy(float) > RAO_DEWETTING_THRESHOLD_KJMOL
                verdict = np.where(np.isfinite(e.to_numpy(float)),
                                   np.where(dewetted, "DEWETTED", "WETTED"), "N/A")
                ready_pred = np.where(verdict == "N/A", np.nan, ~dewetted)
                table[f"{arm}_energy_kJmol"] = e.to_numpy(float)
                table[f"{arm}_verdict"] = verdict
            else:
                thr = merged["fold"].map(
                    lambda f, a=arm: thresholds.get((a, f), np.nan)).to_numpy(float)
                ready_pred = np.where(np.isfinite(prob), prob >= thr, np.nan)
                table[f"{arm}_prob"] = prob
                table[f"{arm}_verdict"] = np.where(
                    np.isfinite(prob), np.where(prob >= thr, "READY", "NOT_READY"), "N/A")
            # Correctness is left missing, not False, on frames the method could
            # not judge
            # and on right-censored frames whose truth is unknown. nanmean below then
            # divides by the frames actually judged, so a method with narrow
            # coverage is
            # not rewarded and not punished for the frames it stayed silent on;
            # how many
            # those were is reported as frames_judged.
            correct = np.where(np.isfinite(ready_pred.astype(float)) & valid,
                               ready_pred.astype(float) == truth.astype(float), np.nan)
            table[f"{arm}_correct"] = correct
            pooled = (evaluation["per_arm"].get(arm, {})
                      .get(str(primary), {}).get("pooled", {}))
            summary_rows.append({
                "method": arm,
                "claims": METHOD_CLAIMS.get(arm, ""),
                "frames_judged": int(np.isfinite(correct.astype(float)).sum()),
                "verdict_accuracy": float(np.nanmean(correct.astype(float)))
                if np.isfinite(correct.astype(float)).any() else float("nan"),
                "balanced_accuracy": _balanced_accuracy(truth[valid],
                                                        ready_pred[valid]),
                "ap_over_chance": pooled.get("ratio_to_chance", float("nan"))
                if pooled.get("defined") else float("nan"),
            })

        # The published structure classifier, when this run's tables carry it.
        if "bl_rao_sigma_d" in features.columns:
            arm = "published_rao2019_sigma"
            sd = table[key].merge(features, on=key,
                                  validate="one_to_one")["bl_rao_sigma_d"].to_numpy(float)
            nonconductive = sd > RAO_SIGMA_D_THRESHOLD
            table[f"{arm}_score"] = sd
            table[f"{arm}_verdict"] = np.where(
                np.isfinite(sd),
                np.where(nonconductive, "NON_CONDUCTIVE", "CONDUCTIVE"), "N/A")
            ready_pred = np.where(np.isfinite(sd),
                                  (~nonconductive).astype(float), np.nan)
            correct = np.where(np.isfinite(ready_pred) & valid,
                               ready_pred == truth.astype(float), np.nan)
            table[f"{arm}_correct"] = correct
            m = valid & np.isfinite(sd)
            # sigma_d grows with the total contour distance past the 1 RT line,
            # so it is
            # large on non-conductive structures: its negative is the ranking
            # score in
            # the same orientation as every other arm, high meaning ready.
            e = ap_with_base(truth[m].astype(int), -sd[m])
            summary_rows.append({
                "method": arm,
                "claims": METHOD_CLAIMS.get(arm, ""),
                "frames_judged": int(np.isfinite(correct.astype(float)).sum()),
                "verdict_accuracy": float(np.nanmean(correct.astype(float)))
                if np.isfinite(correct.astype(float)).any() else float("nan"),
                "balanced_accuracy": _balanced_accuracy(truth[valid],
                                                        ready_pred[valid]),
                "ap_over_chance": e.get("ratio_to_chance", float("nan"))
                if e.get("defined") else float("nan"),
            })

        if answers is not None:
            ours = answers[answers["tau_ps"] == primary][
                key + ["prob_ready", "ready", "mechanism_axis", "reason"]]
            table = table.merge(ours, on=key, how="left", validate="one_to_one")
            table = table.rename(columns={
                "prob_ready": "pcm_prob_ready", "ready": "pcm_verdict",
                "mechanism_axis": "pcm_pore_state", "reason": "pcm_reason"})
        table.to_parquet(out / "benchmark.parquet")
        table.to_csv(out / "benchmark.csv", index=False)

        alarm = json.loads((root / "train" / "alarm.json").read_text())
        head_ev = alarm.get("headline_event_level", {})
        lines = [f"# Method benchmark — {cfg['system.id']} "
                 f"(horizon tau = {int(primary)} ps)", "",
                 "Each method is judged against what it claims to predict, on the "
                 "frames where both its own inputs and the label are defined — "
                 "which is why the frames column is not the same for every row. "
                 "`verdict_accuracy` is the fraction of valid frames "
                 "where the method's READY/NOT_READY reading (or its stated proxy) "
                 "matched whether a crossing actually completed within tau; "
                 "`ap_over_chance` is ranking quality relative to the base rate.", ""]
        lines.append("| method | claims | frames | verdict acc. | balanced acc. | AP / chance |")
        lines.append("|---|---|---|---|---|---|")
        for r in summary_rows:
            lines.append(
                f"| {r['method']} | {r['claims']} | {r['frames_judged']} | "
                f"{r['verdict_accuracy']:.3f} | {r['balanced_accuracy']:.3f} | "
                f"{r['ap_over_chance']:.2f} |")
        lines += ["",
                  "This system additionally reports, per frame: the calibrated "
                  "probability of readiness, the READY/NOT_READY verdict, and the "
                  "diagnosed pore state (mechanism axis with the measured value it "
                  "rests on) — columns `pcm_prob_ready`, `pcm_verdict`, "
                  "`pcm_pore_state`, `pcm_reason` in `benchmark.parquet`/`.csv`."]
        if head_ev:
            lines += ["",
                      f"Event-level monitor (this work only): "
                      f"{head_ev.get('warned_frac_mean', float('nan')):.2f} of crossings "
                      f"warned in advance vs "
                      f"{head_ev.get('random_null_warned_frac_mean', float('nan')):.2f} "
                      f"for a structure-preserving random alarm."]
        (out / "BENCHMARK.md").write_text("\n".join(lines))
        (out / "summary.json").write_text(json.dumps(summary_rows, indent=1,
                                                     ensure_ascii=False))
        log.say(f"benchmark: {len(table)} frames x {len(arms)} methods -> "
                f"benchmark.parquet / BENCHMARK.md")
        log.close()
