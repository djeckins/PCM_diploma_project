"""Splits by trajectory and by time; NEVER by frames.

The unit of generalization is an independent assembly (lineage): replicas from
one equilibrated ancestor form one group. Adjacent frames are near-duplicates,
and a random frame-level split silently inflates every metric.
"""

from __future__ import annotations

import numpy as np


def outer_folds(groups: np.ndarray) -> list[tuple[np.ndarray, np.ndarray, str]]:
    """Leave-one-assembly-out: (train_idx, test_idx, held_out_group).

    groups holds the assembly label of every ROW of the model matrix, so one
    label covers all frames of all replicas descended from the same equilibrated
    ancestor. Held out one at a time rather than in k folds: the number of
    independent assemblies is small, and one fold per assembly is the finest
    split that still keeps a whole assembly on one side.

    The groups are visited in sorted order, which fixes both the sequence of the
    folds and the fold identifier. Per-fold seeds are derived from the single
    recorded seed and that identifier, so the order has to be reproducible.
    """
    out = []
    for g in sorted(set(groups)):
        test = np.flatnonzero(groups == g)
        train = np.flatnonzero(groups != g)
        out.append((train, test, g))
    return out


def inner_rotation(groups: np.ndarray) -> list[tuple[np.ndarray, np.ndarray, str]]:
    """Rotation over the training groups: each in turn serves as validation.

    Hyperparameters are chosen on assemblies the fit did not see, for the same
    reason the outer split exists: a parameter picked on frames of the assembly
    it was fitted on is a parameter tuned to that assembly. Called with the
    training rows of one outer fold, so the held-out assembly never takes part in
    the selection. The construction is the outer one; it carries its own name so
    that the two levels of the nesting are readable at the call sites.
    """
    return outer_folds(groups)


def time_block_folds(times_ps: np.ndarray, n_blocks: int,
                     embargo_ps: float) -> list[tuple[np.ndarray, np.ndarray]]:
    """Time blocks within a single trajectory with a TWO-SIDED embargo.

    embargo_ps = longest feature window + longest label horizon: the feature
    window looks backward (a test frame right after the boundary is built from
    training frames), and the label looks forward (training right before the
    boundary knows about an event inside the test block). An embargo on the
    window alone closes only half of the hole.

    times_ps and embargo_ps are in ps. The returned index arrays are positions in
    the arrays passed in, ordered by time, so the caller can map them straight
    back onto its own rows. This is the fallback used when the assembly rotation
    cannot validate at all, and it never replaces the outer lineage folds.
    """
    order = np.argsort(times_ps)
    # Blocks of equal DURATION, not of equal frame count: the embargo is a time,
    # so a block has to be measured in the same units for the two to be
    # comparable. The blocks are half-open, which leaves the single frame at the
    # upper edge outside every test block.
    edges = np.linspace(times_ps.min(), times_ps.max(), n_blocks + 1)
    folds = []
    for b in range(n_blocks):
        lo, hi = edges[b], edges[b + 1]
        test = order[(times_ps[order] >= lo) & (times_ps[order] < hi)]
        train = order[(times_ps[order] < lo - embargo_ps)
                      | (times_ps[order] >= hi + embargo_ps)]
        # A trajectory short against its embargo can leave one side empty; such a
        # block yields no fold rather than a fold trained on nothing.
        if len(test) and len(train):
            folds.append((train, test))
    return folds


def check_no_frame_leak(train_idx: np.ndarray, test_idx: np.ndarray,
                        groups: np.ndarray) -> None:
    """No trajectory may appear in both train and test at once.

    Called with the per-row TRAJECTORY label, one level finer than the assembly
    the split was built on: disjoint assemblies already imply disjoint
    trajectories, so a violation here means the index arrays, not the grouping,
    went wrong. It raises rather than warns, because the symptom of frame-level
    leakage is a metric that looks good.
    """
    inter = set(groups[train_idx]) & set(groups[test_idx])
    if inter:
        raise RuntimeError(f"groups {inter} on both sides of the split — data leakage")
