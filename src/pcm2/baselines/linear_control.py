"""Linear control on descriptors.

Ridge-penalised logistic regression on the same structural inputs as the
published criterion (constriction radius and lining hydrophobicity). It is
named for that, and not for functional mode analysis, which regresses on
superimposed Cartesian coordinates.
"""

# Radius in A and hydrophobicity on the Wimley-White scale: the same two numbers
# the published criterion reads off its surface (which takes the radius in nm).
# Fitting them instead of looking them up separates two questions — whether the
# published surface is the right function of these inputs, and whether these
# inputs carry the signal at all.
CONTROL_COLS = ["geo_r_constriction_A", "hyd_lining_hydrophobicity"]


def control_columns(available: list[str]) -> list[str]:
    """The control columns present in `available`, in the fixed order declared above.

    Intersected for the same reason as the clock columns: a system without the
    geometry or hydration block has fewer of them, and an empty result makes the
    caller skip the arm by name.
    """
    return [c for c in CONTROL_COLS if c in available]
