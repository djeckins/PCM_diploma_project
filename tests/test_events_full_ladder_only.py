"""Defect classes: a crossing counted without a complete ladder pass, a bilayer
bypass past the pore counted, a left-censored ion given a fake entry time."""

import numpy as np

from pcm2.events import ladder_pass, zone_planes


def _run(z_series, r_val=1.0):
    T = len(z_series)
    times = np.arange(T, dtype=float) * 10
    z = np.asarray(z_series, float)[:, None]
    r = np.full((T, 1), r_val)
    p_low = np.full(T, -20.0)
    p_high = np.full(T, 20.0)
    pl, cl, ch, ph = zone_planes(p_low, p_high)
    return ladder_pass(times, z, r, pl, ph, cl, ch, cylinder_A=8.0)


def test_full_upward_pass_counts_once():
    ev, n, left = _run([-25, -15, 0, 15, 25])
    assert n[0] == 1 and len(ev) == 1
    assert ev[0]["direction"] == 1
    assert ev[0]["t_exit_ps"] == 40.0
    assert ev[0]["entrance_observed"]


def test_partial_pass_does_not_count():
    ev, n, _ = _run([-25, -15, 0, -15, -25])
    assert n[0] == 0 and len(ev) == 0


def test_bilayer_bypass_outside_cylinder_does_not_count():
    ev, n, _ = _run([-25, -15, 0, 15, 25], r_val=9.5)  # r > cylinder in the channel slab
    assert n[0] == 0 and len(ev) == 0


def test_left_censored_start_flagged_and_no_fake_entry():
    ev, n, left = _run([0, 15, 25])  # particle starts inside
    assert left[0]
    assert len(ev) == 0  # origin unknown: no complete pass was observed


def test_downward_pass_negative_direction():
    ev, n, _ = _run([25, 15, 0, -15, -25])
    assert len(ev) == 1 and ev[0]["direction"] == -1


def test_zone_planes_only_from_headgroups():
    pl, cl, ch, ph = zone_planes(np.array([-20.0]), np.array([20.0]))
    assert pl[0] == -20 and ph[0] == 20 and cl[0] == -10 and ch[0] == 10


def test_pbc_wrap_is_not_a_crossing():
    # A bulk ion jumps across the periodic boundary: 4 → 0 without passing the membrane.
    zs = [25.0, 28.0, -28.0, -25.0]
    T = len(zs)
    times = np.arange(T, dtype=float) * 10
    z = np.asarray(zs, float)[:, None]
    r = np.full((T, 1), 1.0)
    p_low, p_high = np.full(T, -20.0), np.full(T, 20.0)
    pl, cl, ch, ph = zone_planes(p_low, p_high)
    box = np.full(T, 60.0)
    ev, n, _ = ladder_pass(times, z, r, pl, ph, cl, ch, cylinder_A=8.0, box_z=box)
    assert n[0] == 0 and len(ev) == 0
    # A jump smaller than half the box is not an event either: a complete pass
    # requires observing the particle in the channel zone between extreme zones.
    ev2, n2, _ = ladder_pass(times, z, r, pl, ph, cl, ch, cylinder_A=8.0, box_z=None)
    assert n2[0] == 0


def test_rigid_frame_shift_cancels_in_relative_axial():
    """Defect class: a whole-frame periodic shift (unstable unwrap image choice)
    teleports absolute coordinates and resets particle origins mid-pass.
    Anchor-relative minimum-image coordinates must cancel the shift exactly."""
    from pcm2.events import relative_axial
    box = 90.0
    z = np.array([10.0, -30.0, 44.0])
    anchor = 5.0
    base = relative_axial(z, anchor, box)
    shifted = relative_axial(z + box, anchor + box, box)  # rigid frame shift
    assert np.allclose(base, shifted)
    # A particle wrapped to the far image resolves to the near one.
    assert np.isclose(relative_axial(np.array([anchor + box - 1.0]), anchor, box)[0], -1.0)


def test_direct_extreme_jump_without_channel_observation_is_not_a_crossing():
    # 0 → 4 in a single sampling step, channel zone never observed: no event counted.
    zs = [-25.0, -22.0, 22.0, 25.0]
    T = len(zs)
    times = np.arange(T, dtype=float) * 10
    z = np.asarray(zs, float)[:, None]
    r = np.full((T, 1), 1.0)
    pl, cl, ch, ph = zone_planes(np.full(T, -20.0), np.full(T, 20.0))
    ev, n, _ = ladder_pass(times, z, r, pl, ph, cl, ch, cylinder_A=8.0)
    assert n[0] == 0 and len(ev) == 0
