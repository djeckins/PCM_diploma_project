"""Dynamics block: drifts and residence times in a state.

A stalled single-file column is a failure route in which geometry and occupancy
look fine while the system does not move; these columns measure it.
Drift is a sign-carrying quantity, canonicalized along the conduction direction.
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, SchemaCtx, backward_window_slice

BLOCK = "dynamics"


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: one displacement and two state ages.

    The drift is in angstroms over one backward window, the ages in picoseconds. There
    is no per-frame part: everything here needs two times, so the block only has a post
    step.
    """
    return [
        ColSpec("dyn_ion_persist", BLOCK, "0/1",
                "the drift's measured object exists: the current innermost ion was in "
                "the pore at the window start as well (companion indicator)",
                "never missing"),
        ColSpec("dyn_ion_drift_A_win", BLOCK, "A",
                "axial displacement of the current constriction-nearest ion over the "
                "backward window; requires the ion to have been in the pore at the "
                "window start too; sign-carrying → canonicalized",
                "ion not in the pore now, or was not at the window start",
                indicator="dyn_ion_persist", sign_flip=True, conditional=True),
        ColSpec("dyn_wet_state_age_ps", BLOCK, "ps",
                "time in the current wet/dry state (from hyd_min_bin_water_n>0)",
                "never missing"),
        ColSpec("dyn_occ_change_age_ps", BLOCK, "ps",
                "time since the last change of occ_n_ions_pore", "never missing"),
    ]


def post(ctx: SchemaCtx, times_ps: np.ndarray, cols: dict[str, np.ndarray]) -> None:
    """Fill the drift and the two state ages for a whole replica.

    The three auxiliary arrays are per-ion histories the replica computer collected while
    it walked the frames; the columns of ion_z_all and ion_in_pore_all are the permeant
    ions in a fixed order, so a column index identifies the SAME ion at every time and
    the difference below is that ion's own displacement.
    """
    ion_z_all = ctx.params["_aux_ion_z_all"]          # [T, N] axial ion coordinates
    ion_in_pore_all = ctx.params["_aux_ion_in_pore"]  # [T, N]
    innermost = ctx.params["_aux_innermost"]          # [T] index or -1
    T = len(times_ps)
    for i in range(T):
        j = int(innermost[i])
        persist = False
        if j >= 0:
            sl = backward_window_slice(times_ps, i, ctx.window_ps)
            k0 = sl.start
            # Two conditions, both required. k0 < i means the window spans a real
            # interval, so at the very first frames there is no drift to report
            # rather
            # than a drift of zero. The pore membership at the window start rules out
            # measuring a displacement across an entry: an ion that arrived
            # mid-window
            # would report the length of its approach rather than motion inside
            # the pore.
            persist = bool(k0 < i and ion_in_pore_all[k0, j])
            if persist:
                # Anchor-relative minimum-image coordinates, so this difference
                # is not
                # disturbed by a rigid shift of the frame or a change of
                # periodic image.
                # The sign is still the laboratory one here; the
                # canonicalization to the
                # conduction direction happens once, after all blocks have run.
                cols["dyn_ion_drift_A_win"][i] = float(ion_z_all[i, j] - ion_z_all[k0, j])
        cols["dyn_ion_persist"][i] = float(persist)
    # Ages are accumulated forward, adding the real frame spacing rather than
    # counting
    # frames, so they stay in picoseconds under any stride and reset the moment
    # the state
    # changes. Frame 0 starts at zero: the age before the trajectory began is
    # unknown,
    # and this column reports the age observed within it.
    wet = cols["hyd_min_bin_water_n"] > 0
    occ = cols["occ_n_ions_pore"]
    age_wet = np.zeros(T)
    age_occ = np.zeros(T)
    for i in range(1, T):
        dt = times_ps[i] - times_ps[i - 1]
        age_wet[i] = age_wet[i - 1] + dt if wet[i] == wet[i - 1] else 0.0
        age_occ[i] = age_occ[i - 1] + dt if occ[i] == occ[i - 1] else 0.0
    cols["dyn_wet_state_age_ps"][:] = age_wet
    cols["dyn_occ_change_age_ps"][:] = age_occ
