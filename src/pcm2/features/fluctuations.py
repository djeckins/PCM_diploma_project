"""Fluctuations block: variances of counters as a leading indicator.

Fluctuations of the number of water molecules in an observation volume are a
measure of surface hydrophobicity: Patel A.J., Varilly P., Chandler D. (2010)
J. Phys. Chem. B 114(4):1632–1637, doi:10.1021/jp909048f.
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, SchemaCtx, backward_window_slice

BLOCK = "fluctuations"


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: windowed variances of three counters already measured.

    Units are the square of the counted quantity, so squared counts for the water and
    ion columns and A^2 for the constriction radius. Nothing here is measured from
    coordinates: the block is a second moment of columns the other blocks wrote.
    """
    return [
        ColSpec("flu_water_pore_var_win", BLOCK, "count^2",
                "variance of hyd_water_pore_n over the strictly backward window",
                "never missing"),
        ColSpec("flu_ion_pore_var_win", BLOCK, "count^2",
                "variance of occ_n_ions_pore over the strictly backward window",
                "never missing"),
        ColSpec("flu_r_constr_var_win", BLOCK, "A^2",
                "variance of geo_r_constriction_A over the window (measured frames only)",
                "no measured radius values in the window",
                indicator="geo_n_search", conditional=True),
    ]


def post(ctx: SchemaCtx, times_ps: np.ndarray, cols: dict[str, np.ndarray]) -> None:
    """Fill the three variance columns for a whole replica.

    Requires the hydration, occupancy and geometry columns to be complete, which the
    caller guarantees by running every per-frame block before any post step. The window
    is the same backward window as everywhere else, so a fluctuation and the level it is
    supposed to run ahead of are measured on one time scale. np.var is the population
    variance, so the first frames of a trajectory, whose window holds one frame, get an
    exact zero rather than a missing value.
    """
    for i in range(len(times_ps)):
        sl = backward_window_slice(times_ps, i, ctx.window_ps)
        # The two counters are never missing, so their variance needs no mask.
        cols["flu_water_pore_var_win"][i] = float(np.var(cols["hyd_water_pore_n"][sl]))
        cols["flu_ion_pore_var_win"][i] = float(np.var(cols["occ_n_ions_pore"][sl]))
        # The radius can be absent on frames where the profile estimator refused
        # every
        # slice. Those frames are dropped instead of being read as a radius, so the
        # variance is over the measured frames of the window only, and the
        # column keeps
        # geo_n_search as its indicator.
        r = cols["geo_r_constriction_A"][sl]
        r = r[np.isfinite(r)]
        if len(r):
            cols["flu_r_constr_var_win"][i] = float(np.var(r))
