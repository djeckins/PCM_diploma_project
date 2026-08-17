"""Van der Waals radii: a fallback ladder with provenance.

Matching order: (resname, atom name) → atom name → element. The first two
rungs hold sourced values only; unattributed constants are not admitted. An
unrecognized atom is a refusal rather than a carbon substitute, since a wrong
radius on one element shows up nowhere else in the run.
"""

from __future__ import annotations

import numpy as np


class VdwError(Exception):
    """Raised when an atom matches no rung of the ladder; it never resolves to a default."""
    pass


# Rung 3: by element.
# Source: Bondi A. (1964) "van der Waals Volumes and Radii",
# J. Phys. Chem. 68(3):441-451, doi:10.1021/j100785a001 — all elements except H.
# H: Rowland R.S., Taylor R. (1996) J. Phys. Chem. 100(18):7384-7391,
# doi:10.1021/jp953141+ — 1.10 A refined from crystallographic contacts.
ELEMENT_RADII_A: dict[str, float] = {
    "H": 1.10,
    "C": 1.70,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "CL": 1.75,
    "K": 2.75,
    "NA": 2.27,
    "MG": 1.73,
    "CA": 2.31,  # no Ca in Bondi; Alvarez S. (2013) Dalton Trans. 42:8617, doi:10.1039/c3dt50599e
}

# Rungs 1-2 are empty: no sourced per-atom tables are vendored. The rung on
# which a match fired is still recorded, and adding tables later leaves the
# interface unchanged.
NAME_RESNAME_RADII_A: dict[tuple[str, str], float] = {}
NAME_RADII_A: dict[str, float] = {}


def resolve_radii(resnames: np.ndarray, names: np.ndarray,
                  elements: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    """Per-atom radius + a counter of fired rungs (recorded in the artifact).

    The three arrays are per-atom and must be in one atom order, the order the caller then
    reads positions in: the profile estimator pairs radii with coordinates by index. → an
    array of radii in A of the same length, and the rung counts. Nothing is returned when
    any atom is unrecognized — the whole call raises, since a run with one wrong radius is
    not distinguishable afterwards from a correct one.
    """
    n = len(names)
    out = np.empty(n, dtype=np.float64)
    steps = {"resname+name": 0, "name": 0, "element": 0}
    unknown: set[str] = set()
    for i in range(n):
        # Keys are upper-cased on both sides: force fields differ in case (Cl,
        # CL, cl) and
        # the case of a name is not chemistry.
        key2 = (str(resnames[i]).upper(), str(names[i]).upper())
        if key2 in NAME_RESNAME_RADII_A:
            out[i] = NAME_RESNAME_RADII_A[key2]
            steps["resname+name"] += 1
            continue
        key1 = str(names[i]).upper()
        if key1 in NAME_RADII_A:
            out[i] = NAME_RADII_A[key1]
            steps["name"] += 1
            continue
        el = str(elements[i]).upper()
        if el in ELEMENT_RADII_A:
            out[i] = ELEMENT_RADII_A[el]
            steps["element"] += 1
            continue
        unknown.add(f"{resnames[i]}/{names[i]}/element={elements[i]!r}")
    # Every unrecognized atom is collected before raising, so one run names the
    # whole gap in
    # the tables; the message is capped at twenty entries to stay readable.
    if unknown:
        raise VdwError(
            "unrecognized atoms (refusal, not a carbon substitute): "
            + ", ".join(sorted(unknown)[:20]))
    return out, steps
