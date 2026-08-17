"""Generate numeric documentation blocks from run artifacts.

Any number that can be computed from an artifact is never typed into a
document. Blocks are delimited by BEGIN/END GENERATED markers; running with
--check compares and fails on divergence — the block-corruption test relies
on the same path.

    python tools/genblocks.py              # rewrite the blocks in docs/*.md
    python tools/genblocks.py --check      # exit 1 if a document has drifted

Requirements: the committed run directories only. Nothing is opened but JSON
artifacts, so both modes work on the published tree with no trajectory data
present, and --check is the quickest way to confirm that the numbers printed in
docs/ are the numbers in runs/. Run it last, after the artifacts are final:
rewriting the documents from a half-finished run would publish intermediate
numbers.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUNS = ROOT / "runs"
DOCS = ROOT / "docs"

MARK = re.compile(r"(<!-- BEGIN GENERATED: (?P<name>[\w-]+) -->)(?P<body>.*?)(<!-- END GENERATED -->)",
                  re.S)


def _j(p: Path):
    """Parse a JSON artifact, or None when the step that writes it has not run."""
    return json.loads(p.read_text()) if p.exists() else None


def _fmt(x, nd=4):
    """Fixed-point number, or an em dash for anything that is not a number.

    A missing or undefined quantity has to print as a dash rather than as 0.0000:
    in these tables a zero would read as a measured value.
    """
    return f"{x:.{nd}f}" if isinstance(x, (int, float)) else "—"


def block_results() -> str:
    """Per-arm metric tables, one section per run, for docs/RESULTS.md.

    The primary horizon is taken from each run's own PROVENANCE.json rather than
    from the current config file, so the table states the horizon the numbers
    were actually computed at even if the config has since been changed.
    """
    lines = [""]
    for run in sorted(RUNS.glob("*-stride*")):
        ev = _j(run / "train" / "evaluation.json")
        if not ev:
            continue
        prov = _j(run / "train" / "PROVENANCE.json")
        primary = str(prov["config"]["labels"]["primary_tau_ps"])
        lines.append(f"### {run.name}")
        lines.append("")
        lines.append("| arm | AP (pooled) | base rate | ×chance |")
        lines.append("|---|---|---|---|")
        for arm, taus in sorted(ev["per_arm"].items()):
            p = taus.get(primary, {}).get("pooled", {})
            if p.get("defined"):
                lines.append(f"| {arm} | {_fmt(p['ap'])} | {_fmt(p['base_rate'])} | "
                             f"{_fmt(p['ratio_to_chance'], 2)} |")
            else:
                lines.append(f"| {arm} | undefined | — | — |")
        alarm = _j(run / "train" / "alarm.json") or {}
        h = alarm.get("headline_event_level")
        if h:
            lines.append("")
            lines.append(f"Monitor: warned {_fmt(h['warned_frac_mean'], 3)} "
                         f"of crossings vs {_fmt(h['random_null_warned_frac_mean'], 3)} "
                         f"for a random alarm of the same structure; empty episodes "
                         f"{h['empty_episodes_total']} of {h['episodes_total']}.")
        anchor = ev.get("anchor_ablation", {})
        rest, tran = anchor.get("resting", {}), anchor.get("transit", {})
        if rest.get("defined"):
            lines.append(f"Anchor ablation: on resting frames AP={_fmt(rest['ap'])} "
                         f"(×{_fmt(rest['ratio_to_chance'], 2)}); on transit frames — "
                         + (f"AP={_fmt(tran['ap'])}" if tran.get("defined") else "undefined")
                         + ".")
        lines.append("")
    return "\n".join(lines) if len(lines) > 1 else "\nno runs yet\n"


def block_acceptance() -> str:
    """PASS/FAIL list of every recorded acceptance gate, per run.

    Only nodes that actually carry an "ok" flag are listed; gates that merely
    record evidence have no verdict to print. A FAIL is published as it stands —
    the point of the block is that a failing gate cannot be left out of the
    documents by hand.
    """
    lines = [""]
    for run in sorted(RUNS.glob("*-stride*")):
        acc = _j(run / "report" / "acceptance.json")
        if not acc:
            continue
        lines.append(f"### {run.name}")
        for k, v in acc.items():
            if isinstance(v, dict) and "ok" in v:
                lines.append(f"- {k}: {'PASS' if v['ok'] else 'FAIL'}")
        lines.append("")
    return "\n".join(lines) if len(lines) > 1 else "\nno runs yet\n"


def block_schema() -> str:
    """The full column table of every run, for docs/SCHEMA.md.

    Column count, units, the meaning of a missing value, the companion indicator
    and the applicability verdict all come from the run's own schema.json and
    applicability.json: the schema is computed from the config, so it differs
    between systems and must never be transcribed by hand.
    """
    lines = [""]
    for run in sorted(RUNS.glob("*-stride*")):
        doc = _j(run / "features" / "schema.json")
        appl = _j(run / "features" / "applicability.json")
        if not doc:
            continue
        lines.append(f"### {run.name} — {len(doc['columns'])} columns (count derived)")
        lines.append("")
        lines.append("| column | block | units | missing means | indicator | verdict |")
        lines.append("|---|---|---|---|---|---|")
        for c in doc["columns"]:
            lines.append(f"| `{c['name']}` | {c['block']} | {c['unit']} | "
                         f"{c['missing_means']} | {c['indicator'] or '—'} | "
                         f"{appl.get(c['name'], '?') if appl else '?'} |")
        lines.append("")
    return "\n".join(lines) if len(lines) > 1 else "\nno runs yet\n"


BLOCKS = {"results": block_results, "acceptance": block_acceptance,
          "schema": block_schema}


def process(check: bool) -> int:
    """Regenerate (or verify) every marked block in docs/; returns the exit code.

    An unknown block name aborts rather than being left untouched: a typo in a
    marker would otherwise turn a generated block into hand-written text that
    nothing keeps in step with the artifacts.
    """
    rc = 0
    for doc in sorted(DOCS.glob("*.md")):
        text = doc.read_text()

        # Replacement keeps the two marker comments (groups 1 and 4) and
        # substitutes only
        # what lies between them, so the markers survive every regeneration and
        # the block
        # stays addressable.
        def repl(m):
            name = m.group("name")
            if name not in BLOCKS:
                raise SystemExit(f"{doc}: unknown generated block {name}")
            return m.group(1) + BLOCKS[name]() + m.group(4)

        new = MARK.sub(repl, text)
        if new != text:
            if check:
                print(f"MISMATCH: {doc} diverges from its source artifacts")
                rc = 1
            else:
                doc.write_text(new)
                print(f"updated {doc}")
    return rc


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    args = ap.parse_args()
    sys.exit(process(args.check))
