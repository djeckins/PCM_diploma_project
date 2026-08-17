"""Figure 2.1 of the thesis: the workflow of the permeation compatibility monitor.

    python tools/workflow_figure.py            # writes runs/FINAL/figure_2_1_workflow.png

Four panels: (a) the input data and the role each system plays, (b) the nine
pipeline steps with what each one computes and writes, (c) what a finished run is
used for, (d) the products. The schematic is drawn from a script rather than in a
drawing program so that renaming a step or changing an output file is a source
edit, and the figure can be regenerated instead of retouched.

Requirements: matplotlib only. Nothing is read from runs/ — the content is the
literal text below, which means it does NOT update itself when the pipeline
changes and has to be edited by hand when a step or an artifact name moves.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "runs" / "FINAL" / "figure_2_1_workflow.png"

INK = "#1b1f24"
MUTED = "#5b6570"
FAINT = "#8b949e"
LINE = "#aab3bd"
POOL = "#2f6fb5"
POOL_FILL = "#eaf1f9"
NEG = "#b4453c"
NEG_FILL = "#f9eceb"
ZERO = "#a67c17"
ZERO_FILL = "#faf3e2"
STEP_FILL = "#f4f6f8"
PANEL_FILL = "#fbfcfd"

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 8.5,
    "text.color": INK,
})

# --------------------------------------------------------------------- panel a
SYSTEMS = [
    ("Gramicidin A", "single-file β-helical dimer, K⁺",
     ["+500 mV · 298 K / 2 M KCl (q 0.70)", "330 K / 1 M KCl (q 0.75) · POPC"],
     "Training pool", POOL, POOL_FILL),
    ("MthK", "tetrameric K⁺, TVGYG filter",
     ["±200 mV · 323 K · CHARMM36 + ECC", "q 0.78 · POPE"],
     "Training pool", POOL, POOL_FILL),
    ("KcsA E71A (K⁺)", "tetrameric K⁺, non-inactivating",
     ["300 mV · 290 K · DOPE/DOPG", "CHARMM36m · q 1.0 (no ECC)"],
     "Training pool", POOL, POOL_FILL),
    ("KcsA arms × 5", "filter mutants, K⁺ or Na⁺",
     ["E71A/Na⁺ · G77A/E71A · T75A/E71A", "with K⁺ or Na⁺ · same protocol"],
     "Negative controls", NEG, NEG_FILL),
    ("Connexin-43", "wide-pore connexin fold, 12 chains",
     ["200 mV · 303 K · POPC · K⁺", "2 replicas used (R3 = duplicate)"],
     "Held-out architecture control", ZERO, ZERO_FILL),
]

# --------------------------------------------------------------------- panel b
STAGES = [
    ("Characterise the system", "topology and a sample of frames", [
        ("1 · autodetect",
         "atom selections from composition and charge, pore axis and\n"
         "extent, coordination cutoff, thermostat temperature, box\n"
         "vectors; refuses rather than guess what it cannot read",
         "autodetect/report.json + config keys"),
    ]),
    ("Measure every frame", "streams the trajectory at the declared stride", [
        ("2 · events",
         "per-ion state machine over five axial zones, inside the\nconfinement cylinder",
         "events/events.parquet"),
        ("3 · features",
         "104 descriptors in ten physical blocks; every window\nlooks strictly backwards",
         "features/features.parquet, schema.json"),
        ("4 · coords",
         "superimposed channel Cα positions, for the PLS-FMA arm",
         "coords/coords.parquet"),
        ("5 · labels",
         "y(t) = 1 if a crossing completes in (t, t+τ];\nincomplete windows censored, not zeroed",
         "labels/labels.parquet"),
    ]),
    ("Learn and report", "splits over independent assemblies, never frames", [
        ("6 · train",
         "leave-one-assembly-out folds with an inner rotation;\nridge and boosted arms, calibration, TreeSHAP",
         "train/oof.parquet, final_model.pkl"),
        ("7 · figures",
         "fold curves, reliability, readiness timelines",
         "figures/fig_*.png"),
        ("8 · report",
         "every declared acceptance check: passed, failed or refused",
         "report/report.md, acceptance.json"),
        ("9 · benchmark",
         "the monitor and the four comparison arms scored on\nidentical frames",
         "benchmark/benchmark.parquet"),
    ]),
]

# --------------------------------------------------------------------- panel c
APPLICATIONS = [
    (ZERO, ZERO_FILL, "Connexin-43 — control on an unrelated architecture",
     "never trained on, scored by the pooled model, so any skill it shows must\n"
     "come from physics that transfers between architectures"),
    (NEG, NEG_FILL, "Five KcsA arms — negative controls",
     "four never conduct and must read not-ready; the G77A/E71A leak is the\n"
     "one rare-positive case where ranking stays definable"),
    (LINE, PANEL_FILL, "Pool ablation — three proteins against four",
     "the pool is refitted with connexin-43 added and both pools scored on\n"
     "identical held-out units: the recorded basis for keeping it outside"),
]

# --------------------------------------------------------------------- panel d
PRODUCTS = [
    ("Cross-protein experiments",
     "the runs behind panel (c), repeated at each of\nthe three horizons",
     "PCM2_GEN_TAU=τ tools/monitor_generalisation.py"),
    ("Tables and figures",
     "every reported number regenerated from the\nrun artefacts",
     "tools/final_results.py  →  runs/FINAL/"),
    ("Interactive monitor",
     "measured state at t, the verdict for (t, t+τ] or a\nwhole range, for any protein with a "
     "config and events",
     "notebooks/permeation_predict.ipynb"),
]


def card(ax, x, y, w, h, *, accent=None, fill="white", lw=0.9, edge=LINE):
    """Rounded box with its lower-left corner at (x, y), in axis data units.

    accent draws a coloured strip inside the top edge, which is how a card
    carries its role colour without a legend.
    """
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=0.09",
        linewidth=lw, edgecolor=edge, facecolor=fill, zorder=1))
    if accent:
        ax.add_patch(FancyBboxPatch(
            (x + 0.06, y + h - 0.20), w - 0.12, 0.13,
            boxstyle="round,pad=0,rounding_size=0.06",
            linewidth=0, facecolor=accent, zorder=2))


def chip(ax, xc, y, text, *, colour, fill):
    """Small filled label centred at (xc, y): the role a card plays."""
    ax.text(xc, y, text, ha="center", va="center", fontsize=7.2, weight="bold",
            color=colour, zorder=4,
            bbox=dict(boxstyle="round,pad=0.40", facecolor=fill, edgecolor="none"))


def arrow(ax, xy_from, xy_to, *, colour=LINE, lw=1.1, rad=0.0):
    """Arrow between two points; rad bows the line so parallel arrows separate.

    The ends are shrunk by a couple of points so the head does not sit on the card
    border it points at.
    """
    ax.add_patch(FancyArrowPatch(
        xy_from, xy_to, arrowstyle="-|>", mutation_scale=9, linewidth=lw,
        color=colour, connectionstyle=f"arc3,rad={rad}", zorder=3,
        shrinkA=2, shrinkB=3))


def panel_letter(ax, x, y, letter, title):
    """Panel label and its heading, on one baseline at (x, y)."""
    ax.text(x, y, letter, fontsize=13, weight="bold", va="center", ha="left")
    ax.text(x + 0.46, y, title, fontsize=10.5, weight="bold", va="center", ha="left")


def draw() -> Path:
    """Lay out all four panels and write the PNG; returns the path written.

    Everything is positioned by hand on a single axis with x in [0, 20] and y in
    [0.55, 18.5], panels stacked downwards in y. Those limits are the layout grid:
    the numbers in the body are positions on it, so moving a panel means changing
    the y band it occupies and the bands below it, and the aspect ratio of the
    figure has to stay as it is for the spacing to hold.
    """
    fig, ax = plt.subplots(figsize=(11.6, 10.6))
    ax.set_xlim(0, 20)
    ax.set_ylim(0.55, 18.5)
    ax.axis("off")

    # ---- (a) input data -------------------------------------------------
    panel_letter(ax, 0.2, 18.15, "a", "Input data: four channel systems, three pore architectures")
    w, gap = 3.72, 0.25
    for i, (name, arch, lines, role, colour, fill) in enumerate(SYSTEMS):
        x = 0.2 + i * (w + gap)
        card(ax, x, 15.55, w, 2.20, accent=colour)
        xc = x + w / 2
        ax.text(xc, 17.28, name, ha="center", fontsize=9.4, weight="bold")
        ax.text(xc, 16.95, arch, ha="center", fontsize=7.3, color=MUTED)
        for j, ln in enumerate(lines):
            ax.text(xc, 16.64 - j * 0.26, ln, ha="center", fontsize=6.7, color=MUTED)
        chip(ax, xc, 15.86, role, colour=colour, fill=fill)

    # ---- (b) the pipeline ----------------------------------------------
    panel_letter(ax, 0.2, 14.95, "b", "One pipeline run per system: input, what it computes, what it writes")
    sw, sgap = 6.3, 0.45
    top, bottom = 14.45, 8.25
    for si, (stage, feeds, steps) in enumerate(STAGES):
        x = 0.2 + si * (sw + sgap)
        card(ax, x, bottom, sw, top - bottom, fill=PANEL_FILL, lw=1.0)
        ax.text(x + 0.25, top - 0.32, stage, fontsize=9.2, weight="bold", va="center")
        ax.text(x + 0.25, top - 0.66, feeds, fontsize=7.2, style="italic",
                color=FAINT, va="center")
        y = top - 1.05
        for name, what, files in steps:
            nlines = what.count("\n") + 1
            h = 0.78 + 0.24 * nlines
            card(ax, x + 0.22, y - h, sw - 0.44, h, fill="white", lw=0.8)
            ax.text(x + 0.40, y - 0.24, name, fontsize=8.4, weight="bold", va="center")
            ax.text(x + 0.40, y - 0.42, what, fontsize=6.8, color=MUTED,
                    va="top", linespacing=1.38)
            ax.text(x + 0.40, y - h + 0.17, files, fontsize=6.5, family="monospace",
                    color=POOL, va="center")
            y -= h + 0.10
        if si == 0:
            ax.text(x + sw / 2, 10.15,
                    "the configuration written here\nfixes every later step: nothing\n"
                    "downstream re-measures or\nre-guesses these values",
                    ha="center", va="center", fontsize=7.2, style="italic",
                    color=FAINT, linespacing=1.6)
        if si < len(STAGES) - 1:
            arrow(ax, (x + sw, (top + bottom) / 2), (x + sw + sgap, (top + bottom) / 2), lw=1.2)
    ax.text(0.2, 7.95, "every step replaces its own directory as a whole — no cache, no partial "
                       "state — and writes PROVENANCE.json beside it: the configuration with the "
                       "origin of every value, a SHA-256 of each source file, library versions",
            fontsize=7.4, style="italic", color=MUTED)

    # ---- (c) evaluation design -----------------------------------------
    panel_letter(ax, 0.2, 7.62, "c", "What a finished run is used for")
    card(ax, 0.2, 4.30, 7.30, 3.00, fill=PANEL_FILL, lw=1.0)
    ax.text(3.85, 7.00, "Conducting pool — three proteins", ha="center",
            fontsize=9.2, weight="bold")
    for j, nm in enumerate(["Gramicidin A", "MthK", "KcsA E71A (K⁺)"]):
        y = 6.40 - j * 0.54
        card(ax, 0.62, y, 6.46, 0.46, fill=POOL_FILL, lw=1.0)
        ax.text(3.85, y + 0.23, nm, ha="center", va="center", fontsize=8.6,
                weight="bold", color=POOL)
    ax.text(3.85, 4.72, "leave-one-protein-out: train on two, score the held-out third\n"
                        "three horizons × full and physics-only column sets",
            ha="center", va="center", fontsize=7.4, color=MUTED, linespacing=1.5)

    card(ax, 0.2, 3.00, 7.30, 1.05, fill="white")
    ax.text(3.85, 3.78, "Benchmarked against", ha="center", fontsize=8.8, weight="bold")
    ax.text(3.85, 3.34, "Rao-2019 dewetting criterion · PLS-FMA (g_fma) · structural\n"
                        "linear control · clock — same frames, labels and splits",
            ha="center", va="center", fontsize=7.4, color=MUTED, linespacing=1.5)
    arrow(ax, (3.85, 4.30), (3.85, 4.05), lw=1.0)

    for k, (colour, fill, title, body) in enumerate(APPLICATIONS):
        y = 5.92 - k * 1.48
        card(ax, 9.60, y, 10.20, 1.34, fill=fill, lw=1.1, edge=colour)
        ax.text(14.70, y + 1.02, title, ha="center", fontsize=9.0, weight="bold",
                color=colour if colour != LINE else INK)
        ax.text(14.70, y + 0.46, body, ha="center", va="center", fontsize=7.4,
                color=MUTED, linespacing=1.5)
    ax.text(8.55, 6.28, "pooled model:\ntrain on all three", ha="center", va="center",
            fontsize=7.2, style="italic", color=MUTED, linespacing=1.4)
    arrow(ax, (7.50, 6.15), (9.60, 6.59), rad=0.05, lw=1.2)
    arrow(ax, (7.50, 5.80), (9.60, 5.11), rad=-0.05, lw=1.2)
    arrow(ax, (7.50, 4.75), (9.60, 3.63), rad=-0.05, lw=1.0, colour="#c2c9d1")

    # ---- (d) products ---------------------------------------------------
    panel_letter(ax, 0.2, 2.42, "d", "Products of a run")
    pw, pgap = 6.30, 0.45
    for k, (title, body, cmd) in enumerate(PRODUCTS):
        x = 0.2 + k * (pw + pgap)
        card(ax, x, 0.72, pw, 1.32, fill="white")
        ax.text(x + pw / 2, 1.83, title, ha="center", fontsize=8.8, weight="bold")
        ax.text(x + pw / 2, 1.42, body, ha="center", va="center", fontsize=7.1,
                color=MUTED, linespacing=1.45)
        ax.text(x + pw / 2, 0.93, cmd, ha="center", va="center", fontsize=6.6,
                family="monospace", color=POOL)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return OUT


if __name__ == "__main__":
    print("written:", draw())
