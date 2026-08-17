""""water_wire" block: continuity and orientation of the water column.

For single-file channels a broken wire is a failure route in its own right;
dipole orientation is a sign-bearing quantity and is canonized along the
conduction direction (the field is flipped, the protein is not).
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, FrameData, SchemaCtx

BLOCK = "water_wire"


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: how continuous the water column is and how it points.

    Gaps are in angstroms along the axis; the dipole columns are cosines against the
    axis, so +1 is a dipole pointing along it and -1 against it, and the mean is the
    only column here that carries a sign.
    """
    return [
        ColSpec("ww_n_waters", BLOCK, "count",
                "number of water oxygens in the pore cylinder within the bin range (indicator)",
                "never missing"),
        ColSpec("ww_max_gap_A", BLOCK, "A",
                "maximum axial gap between adjacent oxygens of the column",
                "fewer than two waters in the pore", indicator="ww_n_waters", indicator_min=2,
                conditional=True),
        ColSpec("ww_continuous", BLOCK, "0/1",
                "1 if max_gap < coordination cutoff (first minimum of the ion–water g(r))",
                "fewer than two waters in the pore", indicator="ww_n_waters", indicator_min=2,
                conditional=True),
        ColSpec("ww_dipole_cos_mean", BLOCK, "cos",
                "mean cos of the water dipole to the axis; sign-bearing, hence canonized",
                "no water in the pore", indicator="ww_n_waters", sign_flip=True,
                conditional=True),
        ColSpec("ww_dipole_cos_std", BLOCK, "cos",
                "spread of pore-water dipole cos — disorder of the column",
                "fewer than two waters in the pore", indicator="ww_n_waters", indicator_min=2,
                conditional=True),
        # System-relative forms: the absolute water count and gap length scale
        # with pore length and get screened out as system labels cross-protein.
        ColSpec("ww_water_linear_density", BLOCK, "A^-1",
                "water oxygens per angstrom of the bin range — system-relative "
                "form of ww_n_waters", "never missing"),
        ColSpec("ww_max_gap_rel", BLOCK, "frac",
                "maximum axial gap as a fraction of the bin range — "
                "system-relative form of ww_max_gap_A",
                "fewer than two waters in the pore", indicator="ww_n_waters",
                indicator_min=2, conditional=True),
    ]


def compute(ctx: SchemaCtx, fd: FrameData, i: int, cols: dict[str, np.ndarray]) -> None:
    """Fill row i of the water-wire columns from the pore water of one frame.

    The column is treated as an ordered chain along the axis, which is why the oxygens
    are sorted first; radial confinement was already applied when fd.wat_z_rel was
    built, so a "gap" here is a gap inside the cylinder and not a detour around it.
    """
    edges = ctx.params["bin_edges"]
    m = (fd.wat_z_rel >= edges[0]) & (fd.wat_z_rel <= edges[-1])
    wz = np.sort(fd.wat_z_rel[m])
    span = float(edges[-1] - edges[0])
    cols["ww_n_waters"][i] = float(len(wz))
    cols["ww_water_linear_density"][i] = len(wz) / span if span > 0 else np.nan
    if len(wz) >= 2:
        gaps = np.diff(wz)
        # Gaps to the mouths count too: the column must reach the edges of the range.
        # Without the two end terms a wire that stops short of a mouth would look
        # continuous, and a pore closed at one end would score like an open one.
        full = np.concatenate([[wz[0] - edges[0]], gaps, [edges[-1] - wz[-1]]])
        mg = float(np.max(full))
        cols["ww_max_gap_A"][i] = mg
        cols["ww_max_gap_rel"][i] = mg / span if span > 0 else np.nan
        # Continuity is judged against the same length that defines a first
        # coordination shell in this system (first minimum of the ion-water g(r)):
        # neighbours closer than that can still hand a hydrogen bond along, a wider
        # break cannot. Reusing the detected cutoff keeps the threshold out of
        # the code.
        cols["ww_continuous"][i] = float(mg < ctx.params["coordination_cutoff_A"])
    if len(fd.wat_dipole_cos):
        # The same axial mask as for the oxygens: fd.wat_dipole_cos is stored in the
        # order of fd.wat_z_rel, one entry per pore-cylinder water.
        cos = fd.wat_dipole_cos[m]
        if len(cos):
            # Mean orientation of the column against the axis. A single water
            # still has
            # an orientation, so this column needs only one; the spread below
            # needs two.
            cols["ww_dipole_cos_mean"][i] = float(np.mean(cos))
        if len(cos) >= 2:
            cols["ww_dipole_cos_std"][i] = float(np.std(cos))
