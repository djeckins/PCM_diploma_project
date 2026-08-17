"""Named-sites block: occupancies of selectivity-filter sites, where a filter exists.

The schema is the union across systems: there are always NS_UNION_SITES columns
(S0–S4, defined in _common); for systems without a filter they are structurally
inapplicable and are masked rather than removed. For a filter system with fewer
rings than the union count, the extra sites are inapplicable in the same way.
"""

from __future__ import annotations

import numpy as np

from ._common import NS_UNION_SITES, ColSpec, FrameData, SchemaCtx

BLOCK = "named_sites"


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: the filter ion count and one soft occupancy per site.

    The count of site columns is NS_UNION_SITES for every system, so tables from a
    filter system and from a filter-less one stay column-compatible; the difference
    between them is carried by applicability, not by a different set of columns.
    """
    cols = [ColSpec("ns_filter_n_ions", BLOCK, "count",
                    "number of permeant ions within the axial range of the filter rings",
                    "no filter (structural mask)")]
    for s in range(NS_UNION_SITES):
        cols.append(ColSpec(f"ns_S{s}_occ", BLOCK, "count",
                            "soft occupancy of a named filter site: "
                            "Σ exp(−(z−z_s)²/2σ²); the center is the midpoint between "
                            "adjacent oxygen rings of the motif, per frame",
                            "no filter, or site beyond the ring count (structural mask)"))
    return cols


def compute(ctx: SchemaCtx, fd: FrameData, i: int, cols: dict[str, np.ndarray]) -> None:
    """Fill row i of the named-site columns; leaves them missing where no filter exists.

    The site centres arrive per frame from fd, as midpoints between adjacent oxygen
    rings of the motif, so a site follows the rings as they move instead of sitting at a
    fixed offset. Returning early leaves the row NaN, which applicability then reports
    as structurally inapplicable rather than as a failed measurement.
    """
    if not ctx.filter_present or fd.site_centers_z is None:
        return
    centers = fd.site_centers_z
    sigma = ctx.params["site_sigma_A"]
    # The filter is read only over ions inside the pore cylinder, as everywhere else.
    zin = fd.ion_z_rel[fd.ion_in_pore]
    # Axial extent of the filter: the outermost site centres widened by one kernel
    # width, so an ion sitting in an end site falls inside the counted range while an
    # ion a full site spacing beyond the filter does not.
    lo, hi = centers.min() - sigma, centers.max() + sigma
    cols["ns_filter_n_ions"][i] = float(np.sum((zin >= lo) & (zin <= hi)))
    # A filter with fewer rings than the union count fills fewer columns; the
    # remaining
    # ones stay NaN and are masked, never zero-filled, since S4 of a four-site filter
    # does not exist rather than standing empty.
    for s in range(min(len(centers), NS_UNION_SITES)):
        cols[f"ns_S{s}_occ"][i] = float(
            np.sum(np.exp(-((zin - centers[s]) ** 2) / (2 * sigma ** 2))))
