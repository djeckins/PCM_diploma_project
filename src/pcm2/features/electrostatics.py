"""Electrostatics block: partner contributions to the potential at the permeant ion.

Coulomb sums within the declared cutoff, in units of e/Å, with topology charges
taken as-is (any ECC charge scaling is part of the replica protocol). No
conversion to kT is applied: a reference temperature is not measurable for
every system (no external MD engine is involved), and mixing normalized and
raw scales in one column is forbidden. The axial field drop is a sign-carrying
quantity and is canonicalized.
"""

from __future__ import annotations

import numpy as np

from ._common import ColSpec, FrameData, SchemaCtx

BLOCK = "electrostatics"
PARTNERS = ("protein", "water", "ions", "lipid")


def schema(ctx: SchemaCtx) -> list[ColSpec]:
    """Column list of the block: one partner contribution per class, plus the asymmetry.

    The rounded cutoff in angstroms goes into the column name, so a table computed with
    another cutoff cannot be mistaken for this one or concatenated with it by accident.
    Partner columns are sums of q/r in e/A, the asymmetry is a net charge in e.
    """
    cut = int(round(ctx.partner_cutoff_A))
    cols = []
    for p in PARTNERS:
        cols.append(ColSpec(f"ele_pot_{p}_cut{cut}_e_per_A", BLOCK, "e/A",
                            f"Σ q_j/r_ij over atoms of class '{p}' within {cut} A of the "
                            "constriction-nearest ion; the cutoff is part of the name",
                            "no innermost ion", indicator="occ_has_innermost",
                            conditional=True))
    cols.append(ColSpec("ele_axial_charge_asym_e", BLOCK, "e",
                        "total charge of mobile particles in the cylinder above the "
                        "constriction minus below; sign-carrying → canonicalized",
                        "constriction not measured", indicator="geo_n_search",
                        sign_flip=True, conditional=True))
    return cols


def compute(ctx: SchemaCtx, fd: FrameData, i: int, cols: dict[str, np.ndarray]) -> None:
    """Fill row i of the electrostatics columns from the sums measured for this frame.

    The sums themselves are computed once per frame by the replica computer, which owns
    the neighbour search; this function only routes them into their named columns. The
    rounding must match schema() exactly, or the write would land on a name that the
    schema never created.
    """
    cut = int(round(ctx.partner_cutoff_A))
    # A class with no atom inside the cutoff still contributes a measured zero,
    # whereas
    # a frame without an innermost ion has nothing to measure at all: the get()
    # default
    # leaves NaN for a partner class the estimator did not report.
    if fd.pot_by_partner is not None:
        for p in PARTNERS:
            cols[f"ele_pot_{p}_cut{cut}_e_per_A"][i] = fd.pot_by_partner.get(p, np.nan)
    cols["ele_axial_charge_asym_e"][i] = fd.axial_charge_asym_e
