"""Symmetry block: symmetry breaking between equivalent subunits.

All quantities are built from distances and are therefore invariant to global
rotation about the pore axis; they are structural and not canonicalized.
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, FrameData, SchemaCtx

BLOCK = "symmetry"


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: spreads of the subunit centres of mass.

    Spreads are population standard deviations in angstroms, over the subunits of the
    frame; sym_nn_dist_cv is their dimensionless coefficient of variation. A perfectly
    symmetric oligomer gives zero, so these columns measure the departure from symmetry
    rather than symmetry itself.
    """
    return [
        ColSpec("sym_com_radial_spread_A", BLOCK, "A",
                "std of radial offsets of subunit centers of mass from the axis",
                "fewer than two subunits (structural mask)"),
        ColSpec("sym_com_z_spread_A", BLOCK, "A",
                "std of axial coordinates of subunit centers of mass; for a dimer, ×2 "
                "gives the axial interface gap (gating by dissociation)",
                "fewer than two subunits (structural mask)"),
        ColSpec("sym_nn_dist_cv", BLOCK, "frac",
                "coefficient of variation of distances between circularly adjacent subunits",
                "fewer than three subunits (structural mask)"),
        ColSpec("sym_com_z_spread_rel", BLOCK, "frac",
                "subunit COM z-spread as a fraction of the bin range — "
                "system-relative form of sym_com_z_spread_A (absolute spreads "
                "scale with pore length and screen out as system labels)",
                "fewer than two subunits (structural mask)"),
    ]


def compute(ctx: SchemaCtx, fd: FrameData, i: int, cols: dict[str, np.ndarray]) -> None:
    """Fill row i of the symmetry columns from the subunit centres of mass of one frame.

    Returns without writing when the assembly has fewer than two subunits: for a single
    chain the quantity does not exist, and applicability marks it inapplicable rather
    than the build treating the gap as a failed measurement.
    """
    if fd.subunit_com is None or len(fd.subunit_com) < 2:
        return
    # Components are laboratory x, y, z with the pore axis taken as z, the same
    # convention the profile estimator uses; every shipped config sets
    # data.pore.axis: z.
    com = fd.subunit_com
    r = np.hypot(com[:, 0], com[:, 1])
    cols["sym_com_radial_spread_A"][i] = float(np.std(r))
    # For a dimer the standard deviation of two numbers is half their
    # separation, which
    # is why the schema notes that twice this column is the axial interface gap.
    cols["sym_com_z_spread_A"][i] = float(np.std(com[:, 2]))
    edges = ctx.params["bin_edges"]
    span = float(edges[-1] - edges[0])
    if span > 0:
        cols["sym_com_z_spread_rel"][i] = float(np.std(com[:, 2])) / span
    if len(com) >= 3:
        # Sorting the subunits by azimuth about the axis puts them in ring
        # order, so the
        # distances below are between CIRCULAR neighbours. Taking a nearest
        # neighbour by
        # distance instead would pair two subunits that had collapsed together
        # and hide
        # exactly the asymmetry this column exists to see.
        ang = np.argsort(np.arctan2(com[:, 1], com[:, 0]))
        ring = com[ang]
        d = np.linalg.norm(np.roll(ring, -1, axis=0) - ring, axis=1)
        # The small term only guards against a division by zero for coincident
        # centres;
        # it is far below any real inter-subunit distance in angstroms.
        cols["sym_nn_dist_cv"][i] = float(np.std(d) / (np.mean(d) + 1e-12))
