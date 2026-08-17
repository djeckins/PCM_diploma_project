"""Temporal ("clock") arm: uses no structure at all.

Predictors are the time since the last completed crossing, a sliding event
rate, and the number of permeant ions in the cylinder. Its score measures how
much of the signal is carried by the burstiness of event arrivals alone; the
gap to the structural arms bounds what conformation contributes.
"""

# Waiting time in ps, event rate in 1/ns, ion count in the pore cylinder.
# None of the three describes the conformation: they are read off the event
# record and the permeant positions, so an arm restricted to them cannot know
# what the channel looks like.
CLOCK_COLS = ["dlv_t_since_cross_ps", "dlv_rate_win", "occ_n_ions_pore"]


def clock_columns(available: list[str]) -> list[str]:
    """The clock columns present in `available`, in the fixed order declared above.

    Intersected rather than assumed: a system whose config leaves the delivery
    block out has no waiting-time column, and the arm then runs on what remains.
    An empty result makes the caller skip the arm by name instead of training on
    an empty matrix.
    """
    return [c for c in CLOCK_COLS if c in available]
