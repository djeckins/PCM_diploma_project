"""Pore profile R(z): largest inscribed sphere in full 3D.

clearance(c) = min_i(|c − a_i| − r_i) over ALL atoms, with their radii;
R(z) = max_c clearance(c), the probe centre restricted to a search disk around
the axis.

Requirements implemented here:
  * the slice grid is anchored to the protein (offsets from the anchor), not to
    the box;
  * admissible centres form the connectivity component reachable from the axis
    (otherwise the probe leaks through the lining into a pocket);
  * "fence": atoms whose spheres reach the slice plane must surround the found
    centre — otherwise a slice at the pore mouth is spuriously wide;
  * a slice on which the estimator refuses is "non-search": it is exported as a
    gap, never filled in by a different estimator;
  * Lipschitz repair: the envelope is built over ALL search slices at once; a
    slice below the envelope is provably not a radius and is dropped as a
    failed measurement.

Atoms: everything that physically occupies volume in the pore cylinder — the
channel AND lipids (acyl-chain entry into the pore is a documented closure
mechanism); water and ions are excluded: the measured quantity is the lumen
available to the mobile phase, and this decision is recorded in the schema.
Hydrogens are included with their radii (accounted for explicitly, and the
convention is recorded).

PBC: the frame is prepared before computation (io.attach_prep) — the channel is
centred in the box and the rest is wrapped by residue; every atom able to
constrain a sphere inside the pore lies in the primary cell near the channel,
so the nearest image is the atom itself.
"""

from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _clearance(cands: np.ndarray, atoms: np.ndarray, radii: np.ndarray,
               z_plane: float) -> np.ndarray:
    """Signed clearance in A of every candidate probe centre on one slice plane.

    cands: [n_cand, 2], in-plane (x, y) positions; atoms: [n_atom, 3]; radii: [n_atom]
    in the atom order of atoms. A candidate centre is the point (x, y, z_plane), so the
    axial coordinate enters only through z_plane. The returned value is
    min_i(|c - a_i| - r_i), the radius of the largest sphere centred at c that touches
    no van der Waals sphere; it is negative when the centre lies inside an atom.

    Distances are taken as raw differences. After frame preparation (io.attach_prep) the
    channel sits at the cell centre, so for every atom able to constrain a probe inside
    the pore the nearest image is the atom itself. The loop is compiled: it runs
    n_cand x n_atom times for every slice of every frame.
    """
    out = np.empty(cands.shape[0])
    for c in range(cands.shape[0]):
        best = 1e9
        cx, cy = cands[c, 0], cands[c, 1]
        for a in range(atoms.shape[0]):
            dx = cx - atoms[a, 0]
            dy = cy - atoms[a, 1]
            dz = z_plane - atoms[a, 2]
            d = (dx * dx + dy * dy + dz * dz) ** 0.5 - radii[a]
            if d < best:
                best = d
        out[c] = best
    return out


class ProbeGrid:
    """Polar grid of candidate centres in the search disk; connectivity counted from the axis."""

    def __init__(self, search_radius_A: float, dr_A: float, n_theta: int):
        """Candidate offsets from the axis: rings dr_A apart, n_theta candidates per ring.

        The grid holds offsets in A, not positions: it is built once per system and every
        slice adds its own axis position to it. max_r_A is the requested search radius
        rounded UP to a whole number of rings, so the disk actually searched is never
        smaller than the one asked for. At least one ring is always built — with the axis
        point alone the grid is a single candidate and the flood fill has nothing to do.

        Angular spacing on ring k is 2*pi*k*dr_A/n_theta, so the grid is finest near the
        axis. That is where the constriction sits, and the constriction is the quantity
        the profile is read for.
        """
        rings = max(1, int(np.ceil(search_radius_A / dr_A)))
        self.max_r_A = rings * dr_A
        pts = [(0.0, 0.0)]
        ring_of = [0]
        theta_of = [0]
        for k in range(1, rings + 1):
            r = k * dr_A
            for j in range(n_theta):
                th = 2 * np.pi * j / n_theta
                pts.append((r * np.cos(th), r * np.sin(th)))
                ring_of.append(k)
                theta_of.append(j)
        self.xy = np.asarray(pts)
        self.ring = np.asarray(ring_of)
        self.theta = np.asarray(theta_of)
        self.n_theta = n_theta
        self.n_rings = rings
        # Adjacency: same ring, ±1 in angle; ring ±1 at the same angle; the
        # centre adjoins ring 1.
        idx = {}
        for i, (k, j) in enumerate(zip(self.ring, self.theta)):
            idx[(k, j)] = i
        adj: list[list[int]] = [[] for _ in range(len(pts))]
        for i, (k, j) in enumerate(zip(self.ring, self.theta)):
            if k == 0:
                adj[i] = [idx[(1, jj)] for jj in range(n_theta)] if rings >= 1 else []
                continue
            adj[i].append(idx[(k, (j + 1) % n_theta)])
            adj[i].append(idx[(k, (j - 1) % n_theta)])
            if k > 1:
                adj[i].append(idx[(k - 1, j)])
            else:
                adj[i].append(0)
            if k < rings:
                adj[i].append(idx[(k + 1, j)])
        self.adj = adj


def slice_radius(grid: ProbeGrid, atoms_xyz: np.ndarray, radii: np.ndarray,
                 axis_xy: np.ndarray, z_plane: float,
                 fence_max_gap_deg: float) -> tuple[float, bool, bool]:
    """→ (R, is_search, on_boundary). is_search=False: the estimator refused this slice.

    R is in A. atoms_xyz [n, 3] and radii [n] are in A and in the same atom order; the
    caller has already restricted them to the slab around this plane. axis_xy is the
    in-plane position of the axis, z_plane the absolute axial coordinate of the plane,
    both in A.

    on_boundary is a diagnostic rather than a failure: the optimum sits on the outermost
    ring, so a wider disk might have found a larger sphere and R is a lower bound. The
    share of such slices is exported as geo_boundary_frac.

    Three conditions make the slice non-search, and all three return NaN so that it is
    exported as a gap and never filled in: the axis is inside an atom, no atom reaches
    the plane, and the fence around the optimum is not closed.
    """
    cands = grid.xy + axis_xy[None, :]
    clr = _clearance(cands, atoms_xyz, radii, z_plane)
    # Candidate 0 is the axis itself. With the axis occluded no admissible centre is
    # connected to it, and the widest sphere elsewhere in the disk would describe a
    # side pocket rather than the lumen.
    if clr[0] <= 0.0:
        return np.nan, False, False
    # Reachability from the axis: flood fill over cells with positive clearance.
    # Only positive-clearance cells are traversed, so the component cannot pass
    # through
    # the lining; a pocket separated from the lumen by even one occluded cell is
    # excluded.
    reach = np.zeros(len(clr), dtype=np.bool_)
    stack = [0]
    reach[0] = True
    while stack:
        i = stack.pop()
        for nb in grid.adj[i]:
            if not reach[nb] and clr[nb] > 0.0:
                reach[nb] = True
                stack.append(nb)
    # argmax within the reachable subset, mapped back to the index in the full grid.
    best = int(np.flatnonzero(reach)[np.argmax(clr[reach])])
    r_best = float(clr[best])
    on_boundary = grid.ring[best] == grid.n_rings
    # Fence: atoms whose spheres reach the slice plane and lie near the found
    # centre must surround it without a large angular gap.
    c_xy = cands[best]
    dz = np.abs(atoms_xyz[:, 2] - z_plane)
    # Elementwise against each atom's own radius: the sphere of atom i cuts the plane
    # when its centre lies closer to the plane than r_i.
    reach_plane = dz < radii
    if not np.any(reach_plane):
        return np.nan, False, on_boundary
    d_xy = atoms_xyz[reach_plane, :2] - c_xy[None, :]
    rad_d = np.hypot(d_xy[:, 0], d_xy[:, 1])
    # How far out a wall atom is looked for. Contact with the probe takes only
    # r_best + r_i,
    # so a ring this wide also admits the atoms just behind the contacting
    # shell, and a
    # single atom short of contact cannot open a gap on its own.
    near = rad_d < (r_best + 2.0 * np.max(radii))
    if not np.any(near):
        return np.nan, False, on_boundary
    # Angles are measured around the found centre and sorted; the last gap closes the
    # circle through ang[0] + 2*pi. A gap wider than fence_max_gap_deg (default
    # 120 deg)
    # means the wall is open in that direction: the plane lies at a pore mouth or the
    # probe has escaped sideways, and either way the value is not a pore radius.
    ang = np.sort(np.arctan2(d_xy[near, 1], d_xy[near, 0]))
    gaps = np.diff(np.concatenate([ang, [ang[0] + 2 * np.pi]]))
    if np.degrees(np.max(gaps)) > fence_max_gap_deg:
        return np.nan, False, on_boundary
    return r_best, True, bool(on_boundary)


def lipschitz_repair(R: np.ndarray, z: np.ndarray, search: np.ndarray,
                     tol_A: float) -> tuple[np.ndarray, int]:
    """R(z) is 1-Lipschitz in z (this holds only across search slices that share
    the same set of admissible centres). The envelope is built over all slices
    at once: consecutive dips must not justify each other.

    Where the bound comes from: a sphere of radius R centred at z, shrunk to
    R - |z' - z| and moved to z', is contained in the original sphere and therefore
    still touches no atom. Hence R(z') >= R(z) - |z' - z| for any two slices, and a
    slice more than tol_A below max_j(R_j - |z - z_j|) cannot be the radius of an
    inscribed sphere. Such a slice is a failed measurement: its search flag is cleared
    and the value stays out of the table. It is never replaced by the envelope value —
    the envelope is an upper bound on what the estimator should have found, not a
    measurement. tol_A (default 0.05 A) absorbs the discretisation of the probe grid.

    R, z, search: [n_slices], z in A and in the order of R. → (search, n_removed)."""
    search = search.copy()
    idx = np.flatnonzero(search)
    if len(idx) < 2:
        return search, 0
    removed = 0
    # Iterate: removing a dip can expose another dip that it was masking.
    changed = True
    while changed:
        changed = False
        idx = np.flatnonzero(search)
        if len(idx) < 2:
            break
        zz, rr = z[idx], R[idx]
        # The maximum includes the slice itself (its own term is rr), so the
        # envelope is
        # never below the measured value and only a deficit against ANOTHER slice can
        # mark a slice bad.
        envelope = np.max(rr[None, :] - np.abs(zz[:, None] - zz[None, :]), axis=1)
        bad = rr < envelope - tol_A
        if np.any(bad):
            search[idx[bad]] = False
            removed += int(np.sum(bad))
            changed = True
    return search, removed


class PoreProfiler:
    """Single profile estimator per system; both descriptors and diagnostics call it."""

    def __init__(self, z_offsets: np.ndarray, search_radius_A: float, dr_A: float,
                 n_theta: int, fence_max_gap_deg: float, lipschitz_tol_A: float,
                 slab_pad_A: float):
        """Fix the slice positions and the three guards for one system.

        z_offsets: slice positions in A as offsets from the anchor, so the grid follows
        the protein and a slice means the same thing in every frame. search_radius_A,
        dr_A and n_theta define the candidate grid (see ProbeGrid); fence_max_gap_deg,
        lipschitz_tol_A and slab_pad_A are the guards described in the module docstring.
        All parameters come from the data.pore section of the config.
        """
        self.z_offsets = z_offsets
        self.grid = ProbeGrid(search_radius_A, dr_A, n_theta)
        self.fence = fence_max_gap_deg
        self.tol = lipschitz_tol_A
        self.slab_pad = slab_pad_A

    def profile(self, atoms_xyz: np.ndarray, radii: np.ndarray,
                anchor: np.ndarray) -> dict:
        """R(z) for one frame → z_offsets, R, search, boundary, n_lipschitz_removed.

        atoms_xyz [n, 3] and radii [n] are in A and in the same atom order; together they
        must cover every volume-occupying atom near the pore (channel and lipids, see the
        module docstring). anchor is a position in A, in practice the centre of mass of
        the anchor selection in this frame.

        R is in A and is NaN wherever search is False. The caller must read those slices
        as absent and put no number in their place; boundary marks slices whose optimum
        sat on the edge of the search disk, and n_lipschitz_removed counts the slices the
        envelope repair discarded.

        The pore axis is the third Cartesian column: this estimator works in the
        laboratory frame with the membrane normal along z, which is what data.pore.axis
        reports for every system in this repository.
        """
        z_planes = anchor[2] + self.z_offsets
        R = np.full(len(z_planes), np.nan)
        search = np.zeros(len(z_planes), dtype=bool)
        boundary = np.zeros(len(z_planes), dtype=bool)
        axis_xy = anchor[:2]
        # Prefilter: measurable clearance is bounded by slab_pad (slices with a
        # larger clearance are pore mouths, rejected by the fence); an atom farther
        # than search disk + slab_pad + its radius cannot constrain a sphere
        # within that bound.
        xy_lim = self.grid.max_r_A + self.slab_pad + float(np.max(radii))
        d_xy = np.hypot(atoms_xyz[:, 0] - axis_xy[0], atoms_xyz[:, 1] - axis_xy[1])
        near = d_xy < xy_lim
        sub_xyz, sub_r = atoms_xyz[near], radii[near]
        for k, zp in enumerate(z_planes):
            dz = np.abs(sub_xyz[:, 2] - zp)
            # Per-atom slab: an atom whose own sphere reaches the plane is kept
            # whatever
            # slab_pad is, and a discarded atom can only be the nearest one when the
            # clearance itself exceeds slab_pad — a mouth slice, which the fence
            # rejects.
            slab = dz < (self.slab_pad + sub_r)
            # No atom near the plane at all: R stays NaN and search stays False.
            # This is
            # how the padding beyond the ends of the channel is reported.
            if not np.any(slab):
                continue
            R[k], search[k], boundary[k] = slice_radius(
                self.grid, sub_xyz[slab], sub_r[slab], axis_xy, zp, self.fence)
        # The repair runs on absolute plane positions; only differences enter the
        # envelope, so the anchor offset cancels and the result is the same as
        # on offsets.
        search, n_removed = lipschitz_repair(R, z_planes, search, self.tol)
        return {"z_offsets": self.z_offsets, "R": R, "search": search,
                "boundary": boundary, "n_lipschitz_removed": n_removed}
