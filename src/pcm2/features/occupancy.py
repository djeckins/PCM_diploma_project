"""Occupancy block: where the ions are — mutual arrangement and desolvation.

Effective conduction proceeds by concerted knock-on, so what matters is not
occupancy per se but the mutual arrangement of the ions and their degree of
desolvation. Soft per-site columns live here only for systems without a
selectivity filter (sites_basis=density_peaks): for filter systems the same
quantities are reported by the named_sites block, and keeping two
implementations of one quantity is forbidden.
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, FrameData, SchemaCtx

BLOCK = "occupancy"


def _density_sites(ctx: SchemaCtx) -> int:
    """How many soft-site columns this block owns; zero unless sites are density peaks.

    A system with a selectivity filter has physical sites and gets them from the
    named_sites block, so returning zero here is what keeps one quantity to one
    estimator instead of two site definitions in one table.
    """
    if ctx.sites_basis == "density_peaks" and ctx.n_sites:
        return ctx.n_sites
    return 0


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: the innermost ion, its solvation, and soft sites.

    Axial quantities are in angstroms, coordination numbers are counts of water
    oxygens, and occ_desolvation is a dimensionless fraction of the same frame's bulk
    coordination.
    """
    cols = [
        ColSpec("occ_n_ions_pore", BLOCK, "count",
                "number of permeant ions in the pore cylinder within the bin range "
                "(indicator for ion-conditional columns)", "never missing"),
        ColSpec("occ_has_innermost", BLOCK, "0/1",
                "the measured object (ion nearest to the constriction) exists: an ion is "
                "in the pore AND the constriction is measured (companion indicator of a "
                "composite object)",
                "never missing"),
        ColSpec("occ_innermost_z_A", BLOCK, "A",
                "axial offset of the constriction-nearest ion from the anchor; "
                "structural coordinate",
                "no innermost ion (no ion in pore, or constriction not measured)",
                indicator="occ_has_innermost", conditional=True),
        ColSpec("occ_innermost_dz_constr_A", BLOCK, "A",
                "|ion z − constriction z| for the innermost ion",
                "no innermost ion", indicator="occ_has_innermost",
                conditional=True),
        ColSpec("occ_innermost_coord_n", BLOCK, "count",
                "water oxygens in the first shell of the innermost ion (cutoff at the "
                "first minimum of g(r))",
                "no innermost ion", indicator="occ_has_innermost", conditional=True),
        ColSpec("occ_desolvation", BLOCK, "frac",
                "1 − coord_n/bulk_coord_n: desolvation degree of the innermost ion; "
                "bulk is the median coordination of ions outside the pore in the same frame",
                "no innermost ion, or no ions outside the pore for the bulk reference",
                indicator="occ_has_innermost", conditional=True),
    ]
    for s in range(_density_sites(ctx)):
        cols.append(ColSpec(f"occ_site{s}_soft", BLOCK, "count",
                            "soft site occupancy: Σ exp(−(z−z_s)²/2σ²) over ions in the "
                            "cylinder; centers are density peaks (a discretization, not "
                            "physical sites)",
                            "never missing: the sum is defined even for an empty pore"))
    if _density_sites(ctx) >= 2:
        cols.append(ColSpec("occ_adjacent_pair_n", BLOCK, "count",
                            "number of adjacent site pairs with both soft occupancies > 1/2",
                            "never missing"))
    return cols


def compute(ctx: SchemaCtx, fd: FrameData, i: int, cols: dict[str, np.ndarray]) -> None:
    """Fill row i of the occupancy columns from one frame.

    The constriction position and the frame's bulk coordination number are read from
    ctx.params, where the replica computer left them: both are measured once per frame
    by their own estimator and are not re-derived here.
    """
    n_in = int(np.sum(fd.ion_in_pore))
    cols["occ_n_ions_pore"][i] = float(n_in)
    cols["occ_has_innermost"][i] = float(fd.innermost is not None)
    if fd.innermost is not None:
        zi = fd.ion_z_rel[fd.innermost]
        cols["occ_innermost_z_A"][i] = float(zi)
        # First-shell count of the ion nearest the gate. The cutoff is the first
        # minimum
        # of this system's own ion-water g(r), detected rather than assumed and
        # recorded
        # as features.coordination_cutoff_A.
        cols["occ_innermost_coord_n"][i] = fd.innermost_coord_n
        zc = ctx.params.get("_z_constr_now", np.nan)
        if np.isfinite(zc):
            cols["occ_innermost_dz_constr_A"][i] = float(abs(zi - zc))
        # Bulk reference from the SAME frame (median coordination of the ions well
        # outside the pore), so the fraction carries no force-field or
        # thermostat offset
        # that a tabulated bulk coordination number would bring in. Not clipped: a
        # negative value means the innermost ion is better coordinated than
        # bulk, which
        # happens and should stay visible.
        bulk = ctx.params.get("_bulk_coord_now", np.nan)
        if np.isfinite(bulk) and bulk > 0 and np.isfinite(fd.innermost_coord_n):
            cols["occ_desolvation"][i] = float(1.0 - fd.innermost_coord_n / bulk)
    ns = _density_sites(ctx)
    if ns:
        # Gaussian kernel of width sigma = features.site_sigma_A, whose rule is
        # half the
        # site spacing; at that width an ion sitting on a neighbouring centre
        # still adds
        # exp(-2) ~ 0.14 here, which is the leakage this discretization accepts.
        sigma = ctx.params["site_sigma_A"]
        centers = fd.site_centers_z
        occ = np.zeros(ns)
        # Only ions inside the pore cylinder enter the sum: an ion at the same axial
        # position but out in the bulk is not in the site.
        zin = fd.ion_z_rel[fd.ion_in_pore]
        for s in range(ns):
            occ[s] = float(np.sum(np.exp(-((zin - centers[s]) ** 2) / (2 * sigma ** 2))))
            cols[f"occ_site{s}_soft"][i] = occ[s]
        if ns >= 2:
            # Half a unit of occupancy is reached either by one ion within
            # sqrt(2 ln 2) sigma ~ 1.18 sigma of the centre or by several further out
            # adding up to it. Two adjacent sites loaded at once is the arrangement a
            # concerted passage proceeds from, which is why the count of such
            # pairs, and
            # not the total occupancy, is the column reported here.
            cols["occ_adjacent_pair_n"][i] = float(np.sum((occ[:-1] > 0.5) & (occ[1:] > 0.5)))
