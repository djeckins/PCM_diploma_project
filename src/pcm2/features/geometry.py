""""geometry" block: pore shape + lining orientation.

Lining carbonyl orientation is the column for the "conformational rearrangement
of the lining" failure route: a C=O flip changes the coordination capacity of a
site without touching the lumen; without this column that diagnostic axis could
never fire.
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, FrameData, SchemaCtx

BLOCK = "geometry"


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: seven whole-pore columns plus three per axial bin.

    The list is a function of ctx.n_bins alone, so it can be produced from the config
    before any trajectory is opened. Radii and axial offsets are in angstroms.
    """
    cols = [
        ColSpec("geo_r_constriction_A", BLOCK, "A",
                "min R(z) over search slices inside the bin range; inscribed-sphere estimator",
                "no search slice inside the pore", indicator="geo_n_search"),
        ColSpec("geo_z_constriction_A", BLOCK, "A",
                "axial anchor offset of the min-R slice; structural quantity, never flipped",
                "no search slice inside the pore", indicator="geo_n_search"),
        ColSpec("geo_n_search", BLOCK, "count",
                "number of search slices inside the bin range (companion indicator)",
                "never missing: the counter is always defined"),
        ColSpec("geo_search_frac", BLOCK, "frac",
                "fraction of search slices per profile frame — estimator-mixing diagnostic",
                "never missing"),
        ColSpec("geo_boundary_frac", BLOCK, "frac",
                "fraction of interior slices with the optimum on the search-disc edge — "
                "checks whether probe_search_radius_A is too small",
                "never missing"),
        ColSpec("geo_lining_carbonyl_cos", BLOCK, "cos",
                "mean cos of the lining backbone C=O angle to the axis in the constriction band; "
                "structural orientation, never flipped (the protein is not flipped)",
                "no lining carbonyls in the band", indicator="geo_lining_n"),
        ColSpec("geo_lining_n", BLOCK, "count",
                "number of lining backbone carbonyls in the constriction band (indicator)",
                "never missing"),
    ]
    for k in range(ctx.n_bins):
        cols.append(ColSpec(f"geo_r_bin{k}_A", BLOCK, "A",
                            "min R(z) over the bin's search slices; bins are fixed anchor offsets",
                            "no search slices in the bin", indicator=f"geo_bin{k}_nsearch"))
    for k in range(ctx.n_bins):
        cols.append(ColSpec(f"geo_bin{k}_nsearch", BLOCK, "count",
                            "number of search slices in the bin (indicator)", "never missing"))
    # System-relative form: absolute bin radii separate architectures by scale
    # (2 A single-file vs 20 A junction) and are screened out as system labels
    # in cross-protein experiments; the ratio to the frame's own constriction
    # is the profile SHAPE, comparable across architectures.
    for k in range(ctx.n_bins):
        cols.append(ColSpec(f"geo_r_bin{k}_rel", BLOCK, "ratio",
                            "bin min-radius over the frame's constriction radius — "
                            "profile shape, system-relative form of geo_r_bin{k}_A",
                            "no search slices in the bin (or no constriction)",
                            indicator=f"geo_bin{k}_nsearch"))
    return cols


def compute(ctx: SchemaCtx, fd: FrameData, i: int, cols: dict[str, np.ndarray]) -> None:
    """Fill row i of the geometry columns from the pore profile of one frame.

    Writes into cols in place. Anything not written stays NaN, which the schema reads
    as "nothing to measure" — a refused slice must never appear as a radius of zero.
    """
    edges = ctx.params["bin_edges"]
    z, R = fd.prof_z_offsets, fd.prof_R
    search = fd.prof_search
    # interior is the bin range the descriptors cover; s_int narrows it to the
    # slices on
    # which the inscribed-sphere estimator actually returned a radius. A refused
    # slice is
    # exported as a gap and is never interpolated over by its neighbours, so
    # every radius
    # below rests on one estimator alone.
    interior = (z >= edges[0]) & (z <= edges[-1])
    s_int = search & interior
    # Taken over the whole slice grid, not the bin range: it diagnoses the estimator
    # (how much of the profile it could measure), not the region the model reads.
    cols["geo_search_frac"][i] = float(np.mean(search)) if len(search) else np.nan
    # A probe optimum sitting on the rim of the search disc means the true
    # maximum may
    # lie outside it, so a large fraction says probe_search_radius_A is set too
    # small.
    cols["geo_boundary_frac"][i] = (float(np.mean(fd.prof_boundary[interior]))
                                    if np.any(interior) else np.nan)
    cols["geo_n_search"][i] = float(np.sum(s_int))
    if np.any(s_int):
        # The narrowest measured slice is the gate of the frame; the shared frame
        # measurement locates the same slice, and the blocks that need only its
        # position
        # read it from there. argmin keeps the first minimum, so ties go to the
        # lowest z
        # and the choice does not depend on the order of the grid.
        j = np.flatnonzero(s_int)[np.argmin(R[s_int])]
        cols["geo_r_constriction_A"][i] = R[j]
        cols["geo_z_constriction_A"][i] = z[j]
    # Bins are half-open, [edges[k], edges[k+1]), so a slice enters at most one
    # of them;
    # a slice landing exactly on the top edge counts in geo_n_search and in no bin.
    for k in range(ctx.n_bins):
        in_bin = s_int & (z >= edges[k]) & (z < edges[k + 1])
        cols[f"geo_bin{k}_nsearch"][i] = float(np.sum(in_bin))
        if np.any(in_bin):
            cols[f"geo_r_bin{k}_A"][i] = float(np.min(R[in_bin]))
    # The shape columns are meaningful only against a measured, non-zero
    # constriction:
    # a division by a missing or vanishing radius would leave an infinity in the
    # table,
    # which the model would read as an ordinary large value.
    rc = cols["geo_r_constriction_A"][i]
    if np.isfinite(rc) and rc > 0:
        for k in range(ctx.n_bins):
            rb = cols[f"geo_r_bin{k}_A"][i]
            if np.isfinite(rb):
                cols[f"geo_r_bin{k}_rel"][i] = float(rb / rc)
    cols["geo_lining_carbonyl_cos"][i] = fd.lining_carbonyl_cos
    cols["geo_lining_n"][i] = float(fd.lining_n_carbonyls)
