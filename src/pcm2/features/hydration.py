""""hydration" block: wetting — the dehydration route of a hydrophobic constriction.

In a narrow apolar pore the barrier is thermodynamic rather than steric, so the
radius must be accompanied by water counts, windowed wetting log-odds and the
hydrophobicity of the lining.
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, FrameData, SchemaCtx, backward_window_slice

BLOCK = "hydration"


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: counts and densities of pore water, plus the lining.

    Counts are numbers of water oxygens (one per molecule), densities are per cubic
    angstrom of measured lumen, and the hydrophobicity is on the vendored Wimley–White
    scale in [-1, 1] with +1 the most hydrophobic residue.
    """
    cols = [
        ColSpec("hyd_water_pore_n", BLOCK, "count",
                "number of water oxygens in the pore cylinder within the bin range",
                "never missing"),
        ColSpec("hyd_min_bin_water_n", BLOCK, "count",
                "minimum water count over bins — dry-spot detector", "never missing"),
        ColSpec("hyd_wet_frac_win", BLOCK, "frac",
                "fraction of window frames (strictly backward) with min_bin_water_n>0",
                "never missing"),
        ColSpec("hyd_wet_logodds_win", BLOCK, "logit",
                "log((n_wet+1)/(n_dry+1)) over the backward window — windowed estimate "
                "of the wetting free energy in ln-odds units", "never missing"),
        ColSpec("hyd_lining_hydrophobicity", BLOCK, "WW[-1,1]",
                "mean hydrophobicity (Wimley–White scale, vendored file) of residues "
                "facing the pore in the constriction band",
                "no pore-facing residues in the band",
                indicator="hyd_lining_n_facing"),
        ColSpec("hyd_lining_n_facing", BLOCK, "count",
                "number of pore-facing residues in the constriction band (indicator)",
                "never missing"),
        # System-relative forms: absolute water counts and residue counts scale
        # with architecture (a 12-chain junction holds hundreds of waters, a
        # single-file channel a handful) and get screened out as system labels
        # in cross-protein experiments; densities and per-subunit counts are
        # the same physics on a shared scale.
        ColSpec("hyd_water_pore_density", BLOCK, "A^-3",
                "water oxygens per searched pore volume (sum of pi*R^2*dz over "
                "search slices in the bin range) — system-relative form of "
                "hyd_water_pore_n",
                "no searched profile volume in the bin range",
                indicator="geo_n_search"),
        ColSpec("hyd_lining_facing_per_sub", BLOCK, "count/subunit",
                "pore-facing residues in the constriction band per subunit — "
                "system-relative form of hyd_lining_n_facing",
                "subunit count unknown"),
    ]
    for k in range(ctx.n_bins):
        cols.append(ColSpec(f"hyd_water_bin{k}_n", BLOCK, "count",
                            "number of water oxygens in the axial bin inside the cylinder",
                            "never missing"))
    for k in range(ctx.n_bins):
        cols.append(ColSpec(f"hyd_water_bin{k}_density", BLOCK, "A^-3",
                            "water oxygens per searched volume of the axial bin — "
                            "system-relative form of the bin count",
                            "no search slices in the bin",
                            indicator=f"geo_bin{k}_nsearch"))
    return cols


def compute(ctx: SchemaCtx, fd: FrameData, i: int, cols: dict[str, np.ndarray]) -> None:
    """Fill row i of the hydration columns from one frame; windowed columns come later.

    fd.wat_z_rel already holds only the oxygens inside the pore cylinder, so the masks
    here select along the axis alone.
    """
    edges = ctx.params["bin_edges"]
    wz = fd.wat_z_rel
    in_range = (wz >= edges[0]) & (wz <= edges[-1])
    cols["hyd_water_pore_n"][i] = float(np.sum(in_range))
    # np.histogram closes the last bin on the right, matching in_range, so the bin
    # counts sum to hyd_water_pore_n exactly and no water is counted twice or lost.
    per_bin = np.histogram(wz[in_range], bins=edges)[0]
    for k in range(ctx.n_bins):
        cols[f"hyd_water_bin{k}_n"][i] = float(per_bin[k])
    # The minimum over bins, not the total, is what a dewetted state shows: one empty
    # bin anywhere along the range breaks the water column even while the pore as a
    # whole still looks well filled. It is also the quantity the wet/dry state
    # is read
    # from, in the windowed columns below and in the age column of the dynamics
    # block.
    cols["hyd_min_bin_water_n"][i] = float(per_bin.min()) if len(per_bin) else np.nan
    ww = ctx.params["ww_scale"]
    # Residues the scale does not cover (ligands, ions modelled as residues) are left
    # out of both the mean and the indicator, so hyd_lining_n_facing is exactly the
    # number of residues the reported hydrophobicity averages over.
    vals = [ww[rn] for rn in fd.lining_resnames if rn in ww]
    cols["hyd_lining_n_facing"][i] = float(len(vals))
    if vals:
        cols["hyd_lining_hydrophobicity"][i] = float(np.mean(vals))
    # Searched pore volume of this frame: pi*R^2*dz summed over search slices.
    # The denominator therefore covers only slices the profiler measured, while the
    # numerator counts every water in the range; on a frame where part of the
    # range was
    # refused the density is an upper bound, and geo_n_search is its indicator.
    z, R, search = fd.prof_z_offsets, fd.prof_R, fd.prof_search
    s_int = search & (z >= edges[0]) & (z <= edges[-1])
    if np.any(s_int):
        # The slice grid is uniform by construction, so the median spacing IS
        # the step;
        # the fallback matters only for a degenerate one-slice profile.
        dz = float(np.median(np.diff(z))) if len(z) > 1 else 1.0
        vol = float(np.sum(np.pi * R[s_int] ** 2) * dz)
        if vol > 0:
            cols["hyd_water_pore_density"][i] = cols["hyd_water_pore_n"][i] / vol
        for k in range(ctx.n_bins):
            in_bin = s_int & (z >= edges[k]) & (z < edges[k + 1])
            if np.any(in_bin):
                v = float(np.sum(np.pi * R[in_bin] ** 2) * dz)
                if v > 0:
                    cols[f"hyd_water_bin{k}_density"][i] = float(per_bin[k]) / v
    if ctx.n_subunits:
        cols["hyd_lining_facing_per_sub"][i] = (
            cols["hyd_lining_n_facing"][i] / ctx.n_subunits)


def post(ctx: SchemaCtx, times_ps: np.ndarray, cols: dict[str, np.ndarray]) -> None:
    """Add the two windowed wetting columns once the per-frame columns are complete.

    Runs per replica, over the whole trajectory at once, and only reads columns this
    block wrote itself; times_ps is in picoseconds and ascending.
    """
    # A frame is wet when no bin of the range is empty: the state is the
    # continuity of
    # the water column, not its total occupancy.
    wet = cols["hyd_min_bin_water_n"] > 0
    for i in range(len(times_ps)):
        sl = backward_window_slice(times_ps, i, ctx.window_ps)
        w = wet[sl]
        n_wet, n_dry = int(np.sum(w)), int(np.sum(~w))
        # Every window contains frame i itself, so the slice is empty only if the
        # trajectory is; the guard keeps the fraction defined in that case.
        cols["hyd_wet_frac_win"][i] = n_wet / max(len(w), 1)
        # One pseudo-count on each side keeps the estimate finite for a window
        # that is
        # entirely wet or entirely dry, where the plain odds would diverge; the
        # price is
        # that a saturated window is shrunk toward zero. The value is in natural-log
        # odds, which for a two-state equilibrium is the wetting free energy in units
        # of kT up to the sign convention (SCHEMA.md records the unit as logit).
        cols["hyd_wet_logodds_win"][i] = float(np.log((n_wet + 1) / (n_dry + 1)))
