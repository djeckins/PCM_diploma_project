"""Run provenance and the human-readable report.

The report is generated as one block; every number in it is computed from the
run's artifacts at generation time.
"""

from __future__ import annotations

import hashlib
import json

from .config import Config
from .runtime import VENDORED_DIR, StepLog, run_dir, step_output


def _j(path):
    """Parse a JSON artifact, or None when the step that writes it has not run.

    A gate about an absent step must be able to report "not present" instead of
    crashing the report: partial runs are a normal state during development, and
    a report is most useful exactly then.
    """
    return json.loads(path.read_text()) if path.exists() else None


def acceptance(cfg: Config, root) -> dict:
    """Acceptance gates that are computable from this run's artifacts.

    Returns a dict keyed by gate name; a gate that has a verdict carries an
    "ok" boolean plus the evidence it was derived from, so a FAIL can be traced
    without rerunning anything. Gates read artifacts only — nothing is refitted,
    and the single numeric limit applied, the admissible fraction of grid-edge
    picks, is taken from the config rather than written in here.
    """
    gates: dict = {}
    folds = _j(root / "train" / "folds.json") or []
    # No tunable hyperparameter may land on a grid boundary more often than the
    # threshold.
    # A boundary pick means the optimum is outside the searched grid, so the
    # value that
    # was selected is a property of the grid, not of the data; the honest
    # response is to
    # widen the grid. Counted per parameter over all folds and horizons.
    edge_counts: dict[str, list[bool]] = {}
    ceiling_hits = []
    exception2_cases = []
    for fa in folds:
        for arm, a in fa.get("arms", {}).items():
            hp = a.get("hyperparams", {})
            if "edges" in hp:
                for k, v in hp["edges"].items():
                    edge_counts.setdefault(f"trees.{k}", []).append(bool(v))
                ceiling_hits.append(bool(hp.get("hit_ceiling")))
            if "edge" in hp:
                on_edge = hp["edge"].get("C_on_edge", hp["edge"].get("on_edge", False))
                param = "linear.C" if "C_on_edge" in hp["edge"] else f"{arm}.selection"
                exc = hp["edge"].get("exception2")
                if on_edge and exc:
                    # Edge of a physically meaningful range: the grid cannot be
                    # widened further, so the hit is recorded and not counted.
                    exception2_cases.append(
                        {"fold": fa["fold"], "tau": fa["tau"], "arm": arm, "note": exc})
                    edge_counts.setdefault(param, []).append(False)
                else:
                    edge_counts.setdefault(param, []).append(bool(on_edge))
    limit = cfg["accept.edge_fraction_max"]
    gates["hyperparam_edges"] = {
        k: {"edge_frac": sum(v) / len(v), "ok": sum(v) / len(v) <= limit}
        for k, v in edge_counts.items() if v}
    gates["edge_exception2_recorded"] = {"ok": True, "cases": exception2_cases}
    # A run that stopped at the round ceiling never triggered early stopping, so the
    # number of boosting rounds is the ceiling rather than a selected quantity
    # and the
    # model may be undertrained. Any hit fails the gate.
    gates["no_ceiling_hits"] = {"ok": not any(ceiling_hits),
                                "hits": int(sum(ceiling_hits))}
    # Regularization penalty expressed in comparable units (effective degrees of
    # freedom).
    # C is an arbitrary scale that says nothing across folds or feature sets;
    # the trace of
    # the hat matrix says how many parameters the penalized fit actually spends.
    dfs = [a["hyperparams"]["df_effective"] for fa in folds
           for arm, a in fa.get("arms", {}).items()
           if "df_effective" in a.get("hyperparams", {})]
    gates["penalty_in_df_units"] = {"ok": len(dfs) > 0, "df_by_fold": dfs}
    # Events per predictor, reported for each branch separately. The effective sample
    # size for a rare, temporally clustered label is the number of crossings, not the
    # number of frames: consecutive positive frames of one crossing carry one
    # piece of
    # information. Recording it beside the column count seen by the branch is
    # what makes
    # an over-parameterized fold visible; no pass/fail threshold is asserted here.
    epp = {}
    for fa in folds:
        for arm, a in fa.get("arms", {}).items():
            lm = a.get("hyperparams", {}).get("leaf_math")
            if lm:
                epp.setdefault(arm, []).append(
                    {"n_events": lm["n_events"], "branch": a.get("branch"),
                     "n_features": a.get("n_features_seen")})
    gates["events_and_features_per_branch"] = epp
    # Calibration must be recorded per fold. Pooling calibration over folds hides the
    # case where each fold is miscalibrated in its own direction and the errors
    # cancel;
    # the gate asks only that the per-fold record exists, since the numbers
    # themselves
    # are reported (slope and intercept) rather than thresholded.
    ev = _j(root / "train" / "evaluation.json") or {}
    head = cfg["model.head"]
    cal = (ev.get("per_arm", {}).get(head, {})
           .get(str(cfg["labels.primary_tau_ps"]), {}).get("calibration_per_fold"))
    gates["calibration_recorded_per_fold"] = {"ok": bool(cal)}
    # The label source must be declared. Everything downstream is a statement about a
    # particular definition of "a crossing happened", so which counter produced the
    # events — this project's state machine or an annotation shipped with the
    # dataset —
    # is part of the result. "auto" is not an admissible answer.
    ev_sum = _j(root / "events" / "summary.json") or {}
    gates["events_source_declared"] = {
        "ok": ev_sum.get("source_declared") in ("own", "provided"),
        "declared": ev_sum.get("source_declared")}
    # Checksums of the vendored reference files. The comparison arm is only a
    # comparison
    # with the published method as long as the tables are the authors' own bytes; the
    # expected digests below are pinned so an edited or re-downloaded file fails
    # the gate
    # instead of quietly changing what the arm computes.
    sums = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(VENDORED_DIR.glob("*.json"))}
    expected = {
        "heuristic_grid.json":
            "bb29c21e1339e71e7a877d5095a4d72c4ac3ab06ce04c863fbcfb407e567c4f0",
        "wimley_white_1996.json":
            "d898884281bf666074ea8c2a5629896b31b9bea6bed4d4389a736ef355542b5c"}
    gates["vendored_checksums"] = {"ok": sums == expected, "actual": sums}
    # The set of produced figures must match the declared catalog.
    cat = _j(root / "figures" / "catalog.json")
    if cat:
        gates["figure_catalog"] = {
            "ok": sorted(cat["produced"]) == sorted(f + ".png" for f in cat["catalog"])}
    # Answers must carry both a probability and a mechanism; specificity is
    # published.
    # A verdict without a named mechanism is not usable for the diagnosis
    # question, and a
    # mechanism axis without its specificity cannot be judged: an axis named on
    # nearly
    # every frame discriminates nothing, which is only visible from the odds ratio.
    gates["answers_with_mechanism"] = {
        "ok": (root / "train" / "answers.parquet").exists()
              and (root / "train" / "mechanism_specificity.json").exists()}
    # Anchor ablation: the conclusion must hold on resting frames alone. On a transit
    # frame the permeant is already between an observed entrance and its exit, so the
    # descriptors can see the crossing ion itself and a model can score well by
    # reading
    # the passage instead of forecasting it. Resting frames are the honest test; the
    # transit number is reported beside it for contrast, not as a result.
    anchor = ev.get("anchor_ablation", {})
    gates["resting_frames_defined"] = {
        "ok": bool(anchor.get("resting", {}).get("defined")),
        "resting": anchor.get("resting"), "transit": anchor.get("transit")}
    # Mask arithmetic: feature assembly crashes on any mismatch, so existence
    # suffices.
    gates["applicability_mask_written"] = {
        "ok": (root / "features" / "applicability.json").exists()}
    gates["note"] = ("linter/test and notebook gates are checked by make test and a "
                     "notebook run, outside this run's artifacts")
    return gates


def run_step(cfg: Config) -> None:
    """Write report/acceptance.json and report/report.md for one run.

    Both come from the artifacts already on disk, so the report can be
    regenerated at any time and will always quote what the train, events and
    figures steps actually wrote. The prose in report.md is fixed; every number
    in it is read, none is typed.
    """
    root = run_dir(cfg)
    with step_output(cfg, "report") as out:
        log = StepLog(out)
        gates = acceptance(cfg, root)
        (out / "acceptance.json").write_text(json.dumps(gates, indent=1,
                                                        ensure_ascii=False, default=str))
        lines = [f"# Run report {cfg['system.id']} (stride "
                 f"{cfg['data.stride']})", ""]
        ev = _j(root / "train" / "evaluation.json") or {}
        primary = str(cfg["labels.primary_tau_ps"])
        lines.append("## Per-arm metrics (pooled, primary horizon)")
        lines.append("")
        lines.append("| arm | AP | base rate | ×chance | n / n+ |")
        lines.append("|---|---|---|---|---|")
        for arm, taus in sorted((ev.get("per_arm") or {}).items()):
            p = taus.get(primary, {}).get("pooled", {})
            if p.get("defined"):
                lines.append(f"| {arm} | {p['ap']:.4f} | {p['base_rate']:.4f} | "
                             f"{p['ratio_to_chance']:.2f} | {p['n']}/{p['n_pos']} |")
            else:
                # An arm whose frames carry one class only has no average
                # precision at
                # all. Printing a dash keeps that distinct from a measured zero,
                # which
                # would read as an arm that ranks and ranks badly.
                lines.append(f"| {arm} | — | — | — | not defined |")
        lines.append("")
        lines.append("## Paired comparison with the head arm (AP difference, paired "
                     "bootstrap over trajectories)")
        for arm, d in (ev.get("paired_vs_head") or {}).items():
            if d.get("defined"):
                lines.append(f"- vs {arm}: Δ={d['mean_diff']:.4f}, "
                             f"CI {d['ci']}, groups {d['n_groups']} — {d['coverage_note']}")
        alarm = _j(root / "train" / "alarm.json") or {}
        h = alarm.get("headline_event_level")
        if h:
            lines.append("")
            lines.append("## Monitor (event-level headline quantity)")
            lines.append(f"- crossings warned: {h['warned_frac_mean']:.3f} "
                         f"(structure-matched random null: "
                         f"{h['random_null_warned_frac_mean']:.3f})")
            lines.append(f"- empty alarm episodes: {h['empty_episodes_total']} of "
                         f"{h['episodes_total']}")
        anchor = ev.get("anchor_ablation", {})
        if anchor:
            lines.append("")
            lines.append("## Anchor ablation")
            for k in ("resting", "transit"):
                a = anchor.get(k, {})
                if a.get("defined"):
                    lines.append(f"- {k}: AP={a['ap']:.4f} at base rate {a['base_rate']:.4f} "
                                 f"(×{a['ratio_to_chance']:.2f})")
        lines.append("")
        lines.append("## Acceptance (computable gates)")
        for k, v in gates.items():
            if isinstance(v, dict) and "ok" in v:
                lines.append(f"- {k}: {'PASS' if v['ok'] else 'FAIL'}")
        (out / "report.md").write_text("\n".join(lines))
        log.say(f"report: {out / 'report.md'}")
        log.close()
