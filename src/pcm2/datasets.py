"""Model matrix: applying the applicability mask and structural exclusions.

The order is strict: the applicability mask decides what enters the matrix (a
verdict about the system, over the whole run); the activity rail decides what
the model sees within a fold. The reverse order would make the matrix
composition depend on the split.

Excluded structurally rather than via config: the baseline score (a competing
predictor) and exact affine functions of other columns. The affine identity is
re-verified on the data, because a rework of the estimator can silently make it
false.

Every column of the schema leaves this module in exactly one of two states:
trained on, or excluded with a recorded reason. The arithmetic is asserted, so a
column cannot disappear from both the matrix and the account of the matrix.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from .config import Config
from .features._common import ColSpec

# The identity of a frame: which trajectory it came from and when. Carried
# alongside the matrix so a score can be traced back to a frame and to a moment
# in the trajectory, and kept out of X, where the time column would be a clock
# the model could read instead of the structure. The clock is measured on its own,
# by the arm that has no structural columns at all.
META_COLS = ("condition", "replica", "time_ps")


def find_affine_duplicates(X: pd.DataFrame, tol: float = 1e-12) -> list[tuple[str, str]]:
    """Pairs (kept, excluded) where the second is an exact affine function of the first.

    This is a search for identities, not for correlation: with tol at the level of
    floating-point noise, only a column that is literally a rescaling of another
    is found. Columns that merely covary strongly are left in the matrix and only
    reported, since dropping one of a correlated pair would be a decision about
    which of two measured quantities matters.
    """
    cols = list(X.columns)
    out = []
    excluded: set[str] = set()
    # Left to right, so in a chain of mutually affine columns the first in schema
    # order is the one kept and the later ones are reported against it.
    for i, a in enumerate(cols):
        if a in excluded:
            continue
        va = X[a].to_numpy()
        for b in cols[i + 1:]:
            if b in excluded:
                continue
            vb = X[b].to_numpy()
            m = np.isfinite(va) & np.isfinite(vb)
            # Three shared points at least: a straight line fits any two points
            # exactly, so a pair that happens to be measured together on only two
            # frames would otherwise be declared an identity. A constant first
            # column has no slope to fit; constants are the activity rail's
            # business, not this one's.
            if m.sum() < 3 or np.nanstd(va[m]) == 0:
                continue
            # Fit the affine relation and check it holds on every row.
            # The residual bound is relative to the scale of the second column,
            # with an absolute floor of tol so that a column whose values are all
            # far below unity is not held to a stricter test than the others.
            k, c = np.polyfit(va[m], vb[m], 1)
            if np.max(np.abs(vb[m] - (k * va[m] + c))) <= tol * max(1.0, np.abs(vb[m]).max()):
                # Identical missingness is part of being a duplicate. Two columns
                # can lie on one line wherever both are measured and still differ
                # in where they are measurable at all; the pattern of missing
                # values is itself a descriptor here, so such a pair carries two
                # different pieces of information and both are kept.
                if np.array_equal(np.isnan(va), np.isnan(vb)):
                    out.append((a, b))
                    excluded.add(b)
    return out


def build_matrix(cfg: Config, table: pd.DataFrame, schema: list[ColSpec],
                 verdict: dict[str, str], log=print) -> dict:
    """The model matrix and the labels that define the split units.

    verdict is the applicability mask from the features step, one entry per schema
    column. Returns the matrix X, the names of its columns, the reason each
    excluded column was excluded, and, per row, the trajectory it belongs to and
    the assembly that trajectory was built from. The two label arrays are what the
    fold machinery groups on: rows are split over assemblies, and scored per
    trajectory.
    """
    reasons: dict[str, str] = {}
    candidates = []
    # One pass over the schema, in schema order, so the column order of X is the
    # documented order and not the order of the parquet file.
    for c in schema:
        # A published-criterion column is a competing predictor, not an input: it
        # stays in the table, where the comparison arm reads it, and is kept out
        # of the matrix so the model cannot learn from its rival's answer.
        if not c.to_model:
            reasons[c.name] = "structural exclusion: baseline score / outside the model"
        # The subject of the measurement does not exist in this system (no
        # selectivity filter, no subunits): the column is missing everywhere, and
        # missing here means "nothing to measure", never zero.
        elif verdict[c.name] == "inapplicable_structural":
            reasons[c.name] = "mask: structurally inapplicable"
        # Measured, present, and the same in every frame of the run. A column with
        # no variation cannot separate ready frames from the rest, and its
        # constancy is a fact about the system that belongs in the record.
        elif verdict[c.name] == "constant":
            reasons[c.name] = "mask: measured and constant"
        else:
            candidates.append(c.name)
    # An explicit copy rather than a view: the caller goes on reading the same
    # feature table -- the comparison arm takes its baseline columns from it -- and
    # nothing done to the matrix should be able to reach back into it.
    X = table[candidates].copy()
    for keep, drop in find_affine_duplicates(X):
        reasons[drop] = f"structural exclusion: exact affine function of {keep}"
        X = X.drop(columns=[drop])
    used = list(X.columns)
    total = len(schema)
    # Used plus excluded must be the whole schema. If a column were ever dropped
    # without a reason being recorded, the published count of modelled columns
    # would stop matching the column dictionary and nothing else would notice.
    assert len(used) + len(reasons) == total, "excluded-column arithmetic does not add up"
    log(f"model matrix: training on {len(used)} of {total} columns; {len(reasons)} excluded:")
    by_reason: dict[str, int] = {}
    for r in reasons.values():
        by_reason[r.split(":")[0]] = by_reason.get(r.split(":")[0], 0) + 1
    for r, n in by_reason.items():
        log(f"  {r}: {n}")
    # The trajectory label: the unit within which a score is computed, and the
    # unit the paired bootstrap resamples. condition/replica is unique by the
    # cross-check in config.py, so the label cannot fuse two trajectories.
    groups = np.array([f"{c}/{r}" for c, r in zip(table["condition"], table["replica"])])
    # The assembly label: the unit of splitting. Several trajectories can share
    # one lineage, and then they must not be separated by a fold boundary. The map
    # is built from the config rather than from the table, so a trajectory that is
    # in the table but not in the config raises here instead of being scored
    # against an assembly nobody declared.
    lineage = {}
    for cond in cfg["data.conditions"]:
        for rep in cond["replicas"]:
            lineage[f"{cond['id']}/{rep['id']}"] = rep["lineage"]
    assemblies = np.array([lineage[g] for g in groups])
    return {"X": X, "columns": used, "excluded": reasons,
            "trajectory": groups, "assembly": assemblies,
            "meta": table[list(META_COLS)].copy()}
