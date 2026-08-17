"""Assemble the thesis results tree: one folder per narrative role.

The raw pipeline artifacts stay untouched under runs/<system>-stride<N>/ —
they are the provenance every number can be traced back to. This script
builds runs/RESULTS/ as the curated view the thesis is written from, with the
same folder skeleton everywhere:

    RESULTS/<role>/
        README.md      the role, its headline numbers, where they come from
        in_system/     the protein's own pipeline run (training proteins and
                       cx43); absent for negative controls, which have no
                       events to train on
        zero_shot/     this role's units predicted by a monitor trained
                       without them: timelines and full unit reports
        benchmark/     the surveyed methods on this role's frames

    RESULTS/cross_protein/   the experiments that span proteins: LOPO at two
                             horizons, the pooled-model horizon sweep, the
                             PLS-FMA self-check

Roles: three training proteins (gramicidin, MthK, KcsA-E71A K+), two
zero-shot tests (G77A-E71A Na+ — a related family; Cx43 — an unrelated one),
two negative controls (non-conducting K+ mutants; the non-conducting Na+
arm). Rerun after any pipeline refresh; missing artifacts are skipped and
listed, so the tree can be rebuilt at any stage of the project.

    python tools/collect_results.py        # writes runs/RESULTS/

Requirements: the committed run directories, nothing else. No trajectory is
opened and no model is refitted, so this runs on the published tree with no
trajectory data present. RESULTS/ is removed and rebuilt from scratch on every
invocation — an artifact that has since disappeared from runs/ disappears here
too, rather than surviving as a stale copy.
"""

from __future__ import annotations

import json
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from pcm2 import config  # noqa: E402
from pcm2.runtime import PROJECT_ROOT, run_dir  # noqa: E402

RESULTS = PROJECT_ROOT / "runs" / "RESULTS"
# Canonical three-pool runs only: the -cx43pool variants are the
# pool-inclusion experiment's own artifact and must not blend into the
# role folders' numbers or timelines.
GEN_DIRS = sorted(p for p in (PROJECT_ROOT / "runs").glob("monitor-generalisation-tau*")
                  if not p.name.endswith(("-cx43pool", "-pinner")))

# role -> (title, system id, condition ids or None for all, narrative)
ROLES = {
    "training_gramicidin": (
        "Gramicidin A", "gramicidin", None,
        "Training protein. A narrow single-file channel; both salt/temperature "
        "conditions conduct, and all ten trajectory units are in the pool."),
    "training_mthk": (
        "MthK", "mthk", None,
        "Training protein. A wide-vestibule K+ channel at both field signs; "
        "all four units are in the pool."),
    "training_kcsa_E71A": (
        "KcsA E71A (K+)", "kcsa", ["E71A"],
        "Training protein: the conducting K+ unit of the KcsA filter-mutant "
        "family. Its sibling mutants stay out of training and serve as "
        "negative controls."),
    "test_zero_shot_G77A_Na": (
        "G77A-E71A with Na+ (zero-shot)", "kcsa_na", ["G77AE71A"],
        "Zero-shot test within a known family: the rare Na+ leak of the G77A "
        "double mutant. The monitor never saw this unit, its protein, or any "
        "Na+ system in training."),
    "test_cx43": (
        "Connexin-43 (zero-shot)", "cx43", None,
        "Zero-shot test on an unrelated architecture: a wide aqueous "
        "gap-junction pore, double membrane, no selectivity filter; the most "
        "distant transfer in this study."),
    "control_negative_K_mutants": (
        "Non-conducting K+ mutants", "kcsa", ["G77AE71A", "T75AE71A"],
        "Negative control: filter mutants that do not conduct K+. A monitor "
        "that has learned physics must stay silent here; READY fractions "
        "above zero would be false alarms."),
    "control_negative_Na_arm": (
        "Non-conducting Na+ arms", "kcsa_na", ["E71A", "T75AE71A"],
        "Negative control: the same filter variants facing Na+, where the "
        "selective filter blocks conduction. The monitor must stay silent."),
}


def _copy(src: Path, dst_dir: Path) -> bool:
    """Copy one artifact if it exists; True when something was copied.

    The return value is what lets a role folder distinguish "the step has not
    run yet" from "the step ran and produced nothing", which is reported in the
    role README rather than left as an empty directory.
    """
    if src.exists():
        dst_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst_dir / src.name)
        return True
    return False


def _section(report_md: str, title: str) -> str:
    """Body of one '## title' section of a run report, without its heading."""
    m = re.search(rf"## {re.escape(title)}\n(.*?)(?=\n## |\Z)", report_md, re.S)
    return m.group(1).strip() if m else ""


def _headline_table(report_md: str) -> str:
    """The per-arm metric table of a run report, quoted verbatim.

    Quoting the report's own markdown instead of recomputing from
    evaluation.json guarantees the role README and the run report carry
    identical numbers even if one of them is regenerated alone.
    """
    m = re.search(r"\| arm \|.*?(?=\n\n)", report_md, re.S)
    return m.group(0) if m else ""


def _events_count(src: Path, conds: list[str] | None) -> int | None:
    """Complete crossings recorded for these conditions, or None if not collected.

    Counts n_own, this project's own state machine, so the number is on one
    definition of a crossing across all roles; conds=None means every condition
    of the system.
    """
    p = src / "events" / "summary.json"
    if not p.exists():
        return None
    ev = json.loads(p.read_text())
    total = 0
    for key, r in ev.get("replicas", {}).items():
        if conds is None or key.split("/")[0] in conds:
            total += r.get("n_own", 0)
    return total


def _gen_numbers(sysid: str, conds: list[str] | None) -> list[str]:
    """LOPO / negative-control readings for this role, across horizons.

    Returns ready-to-paste markdown bullets. Leave-one-protein-out units are
    quoted as AP over the unit's own base rate ("x2.3"), which is the only form
    comparable between units whose base rates differ by orders of magnitude.
    Negative-control units are quoted as a READY fraction and a mean
    probability instead: with no positive frames there is nothing to rank, and
    silence is the quantity of interest there.
    """
    lines = []
    for gd in GEN_DIRS:
        rp = gd / "results.json"
        if not rp.exists():
            continue
        tau = gd.name.split("tau")[1]
        res = json.loads(rp.read_text())
        for fam, blocks in res.get("lopo", {}).items():
            per = blocks.get("full", {}).get("per_unit_ratio", {})
            mine = {u: v for u, v in per.items()
                    if u.startswith(f"{sysid}:")
                    and (conds is None or u.split(":")[1].split("/")[0] in conds)}
            if mine:
                vals = ", ".join(f"`{u.split(':', 1)[1]}` x{v:.2f}"
                                 for u, v in sorted(mine.items()))
                lines.append(f"- τ = {tau} ps, leave-{fam}-out: {vals}")
        neg = res.get("negative_control", {}).get("full", {}).get("per_unit", {})
        mine = {u: v for u, v in neg.items()
                if u.startswith(f"{sysid}:")
                and (conds is None or u.split(":")[1].split("/")[0] in conds)}
        if mine:
            vals = ", ".join(
                f"`{u.split(':', 1)[1]}` READY {v['ready_fraction']:.2f} "
                f"(mean P {v['mean_prob']:.3f})"
                for u, v in sorted(mine.items()))
            lines.append(f"- τ = {tau} ps, negative control: {vals}")
    return lines


def collect_role(role: str, title: str, sysid: str,
                 conds: list[str] | None, narrative: str) -> str:
    """Build RESULTS/<role>/ and return a one-line status for the console.

    conds selects the conditions of the system that belong to this role, or None
    for all of them: one config can serve several roles, as the KcsA config does
    when its E71A condition is a training protein and its two other filter
    mutants are negative controls. Every number written into the role README is
    read from an artifact; the narrative text is the only prose.
    """
    cfg_path = PROJECT_ROOT / "configs" / f"{sysid}.yaml"
    out = RESULTS / role
    out.mkdir(parents=True)
    missing: list[str] = []
    if not cfg_path.exists():
        (out / "README.md").write_text(
            f"# {title}\n\n{narrative}\n\nNo config yet — the system has not "
            "been onboarded.\n")
        return f"{role}: config missing"
    cfg = config.load(cfg_path)
    src = run_dir(cfg)
    is_control = role.startswith("control_")

    # in_system: the protein's own run. Controls have no model by design.
    if not is_control:
        copied = False
        for f in sorted((src / "figures").glob("*.png")):
            copied |= _copy(f, out / "in_system")
        for f in ("figures/catalog.json", "report/report.md",
                  "report/acceptance.json", "features/applicability.json",
                  "train/final_model.json", "events/summary.json"):
            copied |= _copy(src / f, out / "in_system")
        if not copied:
            missing.append("in_system (pipeline run not finished yet)")

    # zero_shot: timelines for this role's conditions + unit reports.
    zs = out / "zero_shot"
    for gd in GEN_DIRS:
        for cond in cfg["data.conditions"]:
            if conds is not None and cond["id"] not in conds:
                continue
            for rep in cond["replicas"]:
                stem = f"{sysid}_{cond['id']}-{rep['id']}.png"
                for kind in ("lopo_gramicidin", "lopo_mthk", "lopo_kcsa_family",
                             "lopo_cx43", "negative_control"):
                    p = gd / f"timeline_{kind}_{stem}"
                    if p.exists():
                        zs.mkdir(parents=True, exist_ok=True)
                        shutil.copy2(p, zs / f"{gd.name.split('-')[-1]}_{p.name}")
    for d in sorted((PROJECT_ROOT / "runs").glob(f"unit-report-{sysid}-*")):
        cond_of_report = d.name.replace(f"unit-report-{sysid}-", "").rsplit("-", 1)[0]
        if conds is not None and cond_of_report not in conds:
            continue
        dest = zs / d.name
        dest.mkdir(parents=True, exist_ok=True)
        for f in sorted(d.glob("*")):
            if f.is_file():
                shutil.copy2(f, dest / f.name)
    if not zs.exists():
        missing.append("zero_shot (no cross-experiment artifacts yet)")

    # benchmark: the surveyed methods on this role's frames.
    bdir = src / "benchmark"
    if (bdir / "benchmark.parquet").exists():
        _copy(bdir / "BENCHMARK.md", out / "benchmark")
        b = pd.read_parquet(bdir / "benchmark.parquet")
        if conds is not None:
            b = b[b["condition"].isin(conds)]
            note = (f"Rows filtered to conditions {conds} of system "
                    f"`{sysid}`; the system-wide tables live in "
                    f"`runs/{src.name}/benchmark/`.\n")
            (out / "benchmark" / "NOTE.md").write_text(note)
        (out / "benchmark").mkdir(parents=True, exist_ok=True)
        b.to_csv(out / "benchmark" / "benchmark.csv", index=False)
    else:
        missing.append("benchmark (step not run yet)")

    # README with the numbers this role contributes.
    n_events = _events_count(src, conds)
    lines = [f"# {title}", "", narrative, ""]
    lines.append(f"System `{sysid}`" + ("" if conds is None else
                 f", conditions {', '.join(conds)}") +
                 (f"; complete crossings recorded: **{n_events}**." if
                  n_events is not None else "; events not collected yet."))
    lines.append("")
    rep_md = (src / "report" / "report.md")
    if not is_control and rep_md.exists():
        t = rep_md.read_text()
        tbl = _headline_table(t)
        if tbl:
            lines += ["## Own run (`in_system/`)", "",
                      "Trajectory-level splits; `x chance` = average precision "
                      "over the base rate.", "", tbl, ""]
        mon = _section(t, "Monitor (event-level headline quantity)")
        if mon:
            lines += ["Event-level monitor: "
                      + mon.replace("\n- ", "; ").lstrip("- "), ""]
    if is_control:
        lines += ["## No in_system folder", "",
                  "A negative control carries (nearly) no crossings, so no "
                  "in-system model can be trained on it. The reading is the one "
                  "below: a monitor trained on the conducting proteins must stay "
                  "silent here.", ""]
    gen = _gen_numbers(sysid, conds)
    if gen:
        lines += ["## Predicted without ever seeing these units (`zero_shot/`)",
                  ""] + gen + [""]
    units = sorted(p.name for p in zs.glob("unit-report-*")) if zs.exists() else []
    if units:
        lines += ["Full held-out unit reports with reality checks: "
                  + ", ".join(f"`{u}/`" for u in units) + ".", ""]
    if missing:
        lines += ["## Not present yet", ""] + [f"- {m}" for m in missing] + [""]
    (out / "README.md").write_text("\n".join(lines))
    n_png = len(list(out.rglob("*.png")))
    return f"{role}: {n_png} figures" + (f", missing: {missing}" if missing else "")


def collect_cross_protein() -> None:
    """Build RESULTS/cross_protein/ from the experiments that span proteins.

    These belong to no single role: leave-one-protein-out at each horizon, the
    negative-control application, and the PLS-FMA self-check. Only the readable
    summaries and the results JSON are copied; the per-frame score tables stay
    in the generating directories, which the role folders point back to.
    """
    out = RESULTS / "cross_protein"
    out.mkdir(parents=True)
    for gd in GEN_DIRS:
        dest = out / gd.name
        dest.mkdir()
        for f in ("GENERALISATION.md", "results.json"):
            _copy(gd / f, dest)
    sc = PROJECT_ROOT / "runs" / "plsfma-selfcheck"
    if sc.exists():
        dest = out / "plsfma-selfcheck"
        dest.mkdir()
        for f in sorted(sc.glob("*")):
            if f.is_file():
                shutil.copy2(f, dest / f.name)
    lines = ["# Cross-protein experiments", "",
             "The material of the generalisation chapter:", "",
             "- `monitor-generalisation-tau*/` — leave-one-protein-out and "
             "negative controls at each common horizon, in both model "
             "variants, with per-unit ratios and timelines.",
             "- `plsfma-selfcheck/` — the PLS-FMA method exercised in its "
             "published regime (a continuous structural quantity), showing "
             "the implementation delivers what the method promises.",
             "- the pooled-model horizon sweep lives in the notebook cache "
             "(`runs/pooled-predict/`); held-out lift runs from x5.95 at "
             "250 ps to x1.66 at 6 ns, physics-only matching or beating the "
             "full variant at every horizon.", ""]
    (out / "README.md").write_text("\n".join(lines))


def main() -> None:
    """Rebuild the whole RESULTS tree and print one status line per role."""
    if RESULTS.exists():
        shutil.rmtree(RESULTS)
    RESULTS.mkdir(parents=True)
    notes = [collect_role(role, *spec) for role, spec in ROLES.items()]
    collect_cross_protein()
    index = ["# Results by role", "",
             "One folder per narrative role. The five protein roles carry the "
             "same skeleton (`in_system/`, `zero_shot/`, `benchmark/`, "
             "`README.md`); the three exceptions are listed below. Raw "
             "pipeline artifacts stay in `runs/<system>-stride<N>/` — every "
             "number here traces back to them.", ""]
    for role, (title, _s, _c, _n) in ROLES.items():
        index.append(f"- [{title}]({role}/README.md)")
    index.append("- [Cross-protein experiments](cross_protein/README.md)")
    index += ["",
              "The two negative controls, `control_negative_K_mutants/` and "
              "`control_negative_Na_arm/`, have no `in_system/`: a "
              "non-conducting arm has no crossings to train on, so there is no "
              "own model to report. `cross_protein/` has none of the three, "
              "because it belongs to no single protein; under its own "
              "`README.md` it holds the readable summaries of the "
              "cross-protein runs "
              "(`monitor-generalisation-tau{1000,2000,4000}/`) and of the "
              "PLS-FMA self-check (`plsfma-selfcheck/`)."]
    (RESULTS / "README.md").write_text("\n".join(index) + "\n")
    for n in notes:
        print(n)
    print(f"index: {RESULTS / 'README.md'}")


if __name__ == "__main__":
    main()
