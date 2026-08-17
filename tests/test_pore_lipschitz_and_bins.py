"""Defect classes: a slice below the Lipschitz envelope is provably not a radius;
bins recomputed per frame denote different places on different frames."""

import numpy as np

from pcm2.features._common import bin_edges
from pcm2.pore import PoreProfiler, lipschitz_repair


def _ring_atoms(z: float, ring_r: float, n: int = 24, atom_r: float = 1.5):
    th = np.linspace(0, 2 * np.pi, n, endpoint=False)
    xyz = np.stack([(ring_r + atom_r) * np.cos(th),
                    (ring_r + atom_r) * np.sin(th),
                    np.full(n, z)], axis=1)
    return xyz, np.full(n, atom_r)


def _cylinder(z_lo=-10, z_hi=10, ring_r=3.0, step=1.0):
    xs, rs = [], []
    for z in np.arange(z_lo, z_hi + step, step):
        a, r = _ring_atoms(z, ring_r)
        xs.append(a)
        rs.append(r)
    return np.vstack(xs), np.concatenate(rs)


def _profiler(z_offsets):
    return PoreProfiler(z_offsets=z_offsets, search_radius_A=4.0, dr_A=0.5,
                        n_theta=24, fence_max_gap_deg=120.0, lipschitz_tol_A=0.05,
                        slab_pad_A=8.0)


def test_inscribed_sphere_matches_analytic_cylinder():
    atoms, radii = _cylinder(ring_r=3.0)
    prof = _profiler(np.arange(-6.0, 6.1, 1.0))
    res = prof.profile(atoms, radii, np.zeros(3))
    got = res["R"][res["search"]]
    assert len(got) > 5
    # Analytic: a ring at (ring_r + atom_r) with atoms of radius atom_r gives an
    # on-axis clearance of ring_r, but the 3D sphere also sees the neighbouring
    # rings: R is slightly smaller on the inter-ring planes.
    assert np.all(got <= 3.0 + 0.05)
    assert np.all(got >= 2.5)


def test_lipschitz_repair_removes_provable_non_radius():
    z = np.arange(0.0, 10.0, 1.0)
    R = np.full_like(z, 5.0)
    R[4] = 1.0  # a 4 A dip at 1 A spacing violates 1-Lipschitz continuity
    search = np.ones(len(z), dtype=bool)
    fixed, n = lipschitz_repair(R.copy(), z, search, tol_A=0.05)
    assert n == 1 and not fixed[4] and fixed.sum() == len(z) - 1


def test_lipschitz_envelope_uses_all_slices_not_neighbours():
    # Consecutive dips mask each other under a pairwise neighbour check.
    z = np.arange(0.0, 10.0, 1.0)
    R = np.full_like(z, 5.0)
    R[4] = 1.0
    R[5] = 1.2
    search = np.ones(len(z), dtype=bool)
    fixed, n = lipschitz_repair(R.copy(), z, search, tol_A=0.05)
    assert n == 2 and not fixed[4] and not fixed[5]


def test_pocket_behind_lining_is_not_pore():
    # A pocket unreachable from the axis: a wide cavity behind the lining must
    # not stand in for the pore radius (connectivity is required from the axis).
    atoms, radii = _cylinder(ring_r=1.2)
    prof = _profiler(np.array([0.0]))
    res = prof.profile(atoms, radii, np.zeros(3))
    assert res["search"][0]
    assert res["R"][0] < 1.3  # the probe did not leak past the lining atoms


def test_bin_edges_are_fixed_offsets():
    e1 = bin_edges(-14.0, 3.5, 8)
    e2 = bin_edges(-14.0, 3.5, 8)
    assert np.array_equal(e1, e2)
    # Guarded defect: edges derived from per-frame extent — here they CANNOT
    # depend on the frame, since the function takes no frame input; bin width is
    # constant across frames by construction.
    assert np.allclose(np.diff(e1), 3.5)


def test_refused_slice_is_nan_not_fallback():
    # A slice with no atoms at all is non-search; the value is NaN, never a
    # substituted estimate from another method.
    atoms, radii = _cylinder(z_lo=-2, z_hi=2, ring_r=3.0)
    prof = _profiler(np.array([0.0, 30.0]))
    res = prof.profile(atoms, radii, np.zeros(3))
    assert res["search"][0] and not res["search"][1]
    assert np.isnan(res["R"][1])
