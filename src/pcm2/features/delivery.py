"""Delivery block: permeant availability — the most frequent failure route.

The channel may be ready while no ion is at the mouth; without these columns
the model would attribute to structure what is explained by delivery.

The entry-mouth columns are not a flip of the structural columns but a separate
estimator that reads the conduction direction of the condition. Canonicalization
does not touch ordered structural blocks, which would destroy the laboratory
frame, so upper/lower stay laboratory columns and the "entry" quantity is
defined on its own.
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, FrameData, SchemaCtx, backward_window_slice

BLOCK = "delivery"


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: mouth counts, the approach, and the crossing history.

    Distances are in angstroms, times in picoseconds and dlv_rate_win in crossings per
    nanosecond. None of these columns is marked sign_flip: the entry side is chosen from
    the conduction direction when the quantity is measured, so flipping it afterwards
    would apply the direction twice.
    """
    return [
        ColSpec("dlv_n_upper_mouth", BLOCK, "count",
                "permeant ions in the upper mouth zone (radius mouth_radius, depth mouth_depth)",
                "never missing"),
        ColSpec("dlv_n_lower_mouth", BLOCK, "count",
                "permeant ions in the lower mouth zone", "never missing"),
        ColSpec("dlv_n_entry_mouth", BLOCK, "count",
                "permeant ions in the ENTRY mouth zone (the side set by the conduction "
                "direction)",
                "never missing"),
        ColSpec("dlv_entry_wide_n", BLOCK, "count",
                "permeant ions in the wide entry zone (3×mouth_radius) — indicator",
                "never missing"),
        ColSpec("dlv_dist_entry_A", BLOCK, "A",
                "distance from the entry-mouth center to the nearest permeant outside the pore",
                "no permeant in the wide entry zone", indicator="dlv_entry_wide_n",
                conditional=True),
        ColSpec("dlv_approach_pair_n", BLOCK, "count",
                "at how many window endpoints (0..2, end ≠ start) the entry distance "
                "is measured (companion indicator of the approach)", "never missing"),
        ColSpec("dlv_approach_A_win", BLOCK, "A",
                "approach of the nearest permeant toward the entry over the backward "
                "window: dist(t−window)−dist(t)",
                "distance not measured at one of the window endpoints",
                indicator="dlv_approach_pair_n", indicator_min=2, conditional=True),
        ColSpec("dlv_t_since_cross_ps", BLOCK, "ps",
                "time since the last completed crossing in this trajectory",
                "before the trajectory's first crossing the measured object does not exist",
                indicator="dlv_crossed_before", conditional=True),
        ColSpec("dlv_crossed_before", BLOCK, "0/1",
                "whether a completed crossing has already occurred in this trajectory "
                "(indicator)",
                "never missing"),
        ColSpec("dlv_rate_win", BLOCK, "1/ns",
                "rate of completed crossings in the backward window", "never missing"),
    ]


def compute(ctx: SchemaCtx, fd: FrameData, i: int, cols: dict[str, np.ndarray]) -> None:
    """Fill row i of the per-frame delivery columns from the mouth counts of one frame.

    upper and lower are laboratory sides and are always present; entry is written only
    for a condition that declares a conduction direction, and its absence then reads as a
    count of zero rather than as a gap, because the zone itself is undefined.
    """
    mc = fd.mouth_counts or {}
    cols["dlv_n_upper_mouth"][i] = float(mc.get("upper", 0))
    cols["dlv_n_lower_mouth"][i] = float(mc.get("lower", 0))
    cols["dlv_n_entry_mouth"][i] = float(mc.get("entry", 0))
    cols["dlv_entry_wide_n"][i] = float(fd.entry_wide_n)
    cols["dlv_dist_entry_A"][i] = fd.entry_dist_A


def post(ctx: SchemaCtx, times_ps: np.ndarray, cols: dict[str, np.ndarray]) -> None:
    """Fill the crossing-history and approach columns for a whole replica.

    exit_times_dir holds the completion times, in picoseconds, of the crossings of this
    replica that run in the conduction direction of the condition; the caller has already
    restricted them to that direction. Everything here looks backward only, which is what
    keeps these columns out of the horizon the label is defined on.
    """
    exits = np.sort(np.asarray(ctx.params.get("exit_times_dir", [])))
    T = len(times_ps)
    for i in range(T):
        # side="right" puts a crossing that completes exactly at times_ps[i] in
        # the past.
        # That matches the label, which counts crossings in the half-open interval
        # (t, t+tau], so no crossing is both a past event here and a future one
        # there,
        # and none falls between the two definitions.
        k = int(np.searchsorted(exits, times_ps[i], side="right"))
        cols["dlv_crossed_before"][i] = float(k > 0)
        if k > 0:
            cols["dlv_t_since_cross_ps"][i] = float(times_ps[i] - exits[k - 1])
        sl = backward_window_slice(times_ps, i, ctx.window_ps)
        t0 = times_ps[sl.start]
        # Crossings strictly after the window start and up to the current time,
        # so one
        # landing exactly on the window edge belongs to the windows that contain it
        # properly rather than to this one as well.
        n_in_win = int(np.sum((exits > t0) & (exits <= times_ps[i])))
        # The rate is per nanosecond over the span actually covered, not over
        # the nominal
        # window: near the start of a trajectory the window is short and
        # dividing by its
        # nominal length would report a rate the data cannot support. The floor only
        # keeps the division defined on the first frame, where the span is zero.
        span_ns = max((times_ps[i] - t0) / 1000.0, 1e-9)
        cols["dlv_rate_win"][i] = n_in_win / span_ns
        d_now = cols["dlv_dist_entry_A"][i]
        # The distance at the window START frame, not the smallest over the
        # window: the
        # column is meant to be a displacement over a known interval.
        d_then = cols["dlv_dist_entry_A"][sl.start]
        pair = 0
        if sl.start < i:
            pair = int(np.isfinite(d_now)) + int(np.isfinite(d_then))
        cols["dlv_approach_pair_n"][i] = float(pair)
        # Both endpoints have to be measured: a difference taken against a frame
        # with no
        # permeant in the entry zone is not an approach, and substituting a
        # large number
        # for the missing end would manufacture one. Sign convention: positive
        # means the
        # nearest permeant moved TOWARD the entry over the window.
        if pair == 2:
            cols["dlv_approach_A_win"][i] = float(d_then - d_now)
