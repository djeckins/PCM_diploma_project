"""Trajectory opening, atom-group resolution, frame preparation.

All per-frame reading is done by MDAnalysis: streaming (one frame in memory),
an offset index for random access, PBC distances in a triclinic cell.

Frame preparation is an explicit, recorded decision:
  1) the channel is made whole as a molecule (unwrap along bonds);
  2) the system is shifted so the channel centre of mass sits at the cell centre;
  3) the remaining molecules are wrapped into the cell WHOLE (by residue).
After this the membrane and the pore lie far from the cell faces, and the
nearest image for all near-pore distances is the atom itself; distances where
this does not hold are computed with MDAnalysis functions given box (minimum
image in three dimensions; orthogonality is not assumed — the gramicidin cell
is hexagonal). No rotational fit is applied: all exported quantities are
computed relative to the per-frame anchor and the membrane axis and do not
depend on rotation about that axis.
The operation is applied once, before descriptor computation.

What this layer is allowed to do: open files, resolve selections into atom
groups, prepare the frame, and refuse input it cannot read. It computes no
descriptor and holds no state between frames, so there is exactly one place
where a change to frame handling can alter a measured value. Lengths are
angstroms and times picoseconds throughout, as MDAnalysis reports them.
"""

from __future__ import annotations

from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis import transformations as trf

from .config import Config


class IoError(Exception):
    """Input that cannot be read or cannot be trusted. Always a refusal, never a fallback."""
    pass


# Config axis name -> column of a position array. The pore axis is the membrane
# normal; everything axial in the feature layer is a signed offset along this
# column, measured from the per-frame anchor.
AXIS_INDEX = {"x": 0, "y": 1, "z": 2}

# Standard atomic masses as force fields tabulate them. Elements in a GROMACS
# .top come out of a name-based guesser that misreads CHARMM type names
# (HAL2 -> "AL"); the mass column is exact, so it is the measurement to trust.
_MASS_ELEMENTS = (
    (1.008, "H"), (12.011, "C"), (14.007, "N"), (15.999, "O"), (15.9994, "O"),
    (30.974, "P"), (32.06, "S"), (35.45, "Cl"), (22.99, "Na"), (39.0983, "K"),
    (24.305, "Mg"), (40.08, "Ca"),
)


def _elements_from_masses(u: mda.Universe) -> None:
    """Set u.atoms.elements from the mass column; refuse if any atom is left unassigned.

    Elements are needed because the van der Waals radius of an atom is resolved
    by element on the last rung of the ladder in vdw.py, and the pore profile is
    an inscribed sphere over those radii: one element read wrong shifts a radius
    and nothing else in the run reveals it.
    """
    masses = u.atoms.masses
    elements = np.full(len(masses), "", dtype=object)
    # 0.05 amu: two orders of magnitude wider than the disagreement between
    # force-field roundings of one element (15.999 against 15.9994) and far
    # narrower than the closest pair of different elements in the table, K and Ca
    # at 0.98 amu apart, so no two windows can overlap. An atom whose mass has
    # been altered on purpose -- a united-atom bead, a deuterium -- falls outside
    # every window and is refused below instead of taking the nearest element.
    for ref, el in _MASS_ELEMENTS:
        elements[np.abs(masses - ref) < 0.05] = el
    if (elements == "").any():
        bad = u.atoms[elements == ""]
        seen = sorted({f"{a.resname}/{a.name}/mass={a.mass:.4f}" for a in bad[:50]})
        raise IoError("cannot assign elements from masses (refusal, not a "
                      f"guess): {', '.join(seen)}")
    u.atoms.elements = elements.astype(str)


def open_replica(rep: dict) -> mda.Universe:
    """Universe for one replica record of the config; a missing file is named and refused.

    Frames are not prepared here: a caller that measures anything geometric must
    also call attach_prep. Frames stay streamed -- only the current one is in
    memory -- so a replica of any length can be read on a workstation.
    """
    top, traj = Path(rep["topology"]), Path(rep["trajectory"])
    if not top.exists():
        raise IoError(f"missing topology {top}")
    if not traj.exists():
        raise IoError(f"missing trajectory {traj}")
    if top.suffix == ".top":
        # GROMACS system topology: the .itp parser handles it, but include
        # paths inside the file are relative to the topology's own directory
        # (resolved through symlinks).
        u = mda.Universe(str(top), str(traj), topology_format="ITP",
                         include_dir=str(top.resolve().parent))
        # The text topology is the case where the guessed elements are wrong, so
        # they are replaced by the mass-based assignment (see _MASS_ELEMENTS).
        _elements_from_masses(u)
        return u
    return mda.Universe(str(top), str(traj))


def resolve_groups(u: mda.Universe, cfg: Config) -> dict[str, mda.AtomGroup]:
    """Groups from config selections; an undefined or empty selection raises.

    Both failures are refusals on purpose. An undefined selection means
    autodetect has not run on this config, and an empty one means the selection
    string does not fit this topology; either way the features that depend on the
    group would come out missing for every frame, the missing-value budget would
    absorb it, and the run would still finish.
    """
    groups: dict[str, mda.AtomGroup] = {}
    for name in ("ion", "ion_negative", "water_oxygen", "water_hydrogen",
                 "phosphate", "lipid", "channel"):
        sel = cfg[f"data.selections.{name}"]
        if sel is None:
            raise IoError(f"selection {name} is undefined: run the autodetect step first")
        ag = u.select_atoms(sel)
        if len(ag) == 0:
            raise IoError(f"selection {name} = {sel!r} is empty on this topology")
        groups[name] = ag
    return groups


class _FragmentImageAligner:
    """Pin every channel fragment to the periodic image nearest the first one.

    Unwrapping makes each fragment whole but chooses its periodic image
    independently, and for a multi-chain channel that choice can alternate
    between frames: the channel centre of mass then jumps by one cell, the
    whole recentred frame shifts with it, and every axial series develops
    box-sized teleports. Aligning fragment images to the first fragment before
    recentring removes the instability at its source; when fragments are
    already together this is an exact no-op.
    """

    def __init__(self, channel: mda.AtomGroup):
        # Fragments are bonded components, so for a multi-chain channel this is
        # the list of chains. Resolved once: the bond graph does not change.
        self.frags = list(channel.fragments)

    def __call__(self, ts):
        # A channel that is one bonded molecule has nothing to align, and the loop
        # below would be an identity on it.
        if len(self.frags) < 2:
            return ts
        from MDAnalysis.lib.mdamath import triclinic_vectors
        # Rows of L are the cell vectors in angstroms, so x @ Linv expresses a
        # displacement in fractional cell coordinates. Working in fractions is
        # what makes this correct for a non-orthogonal cell.
        L = triclinic_vectors(ts.dimensions)
        Linv = np.linalg.inv(L)
        ref = self.frags[0].atoms.center_of_mass()
        for f in self.frags[1:]:
            delta = f.atoms.center_of_mass() - ref
            # Rounding the fractional separation gives the integer image offset;
            # it is zero whenever the chains already sit in the same image, and
            # the subunits of an assembled channel are much closer to each other
            # than half a cell, so a nonzero n means a wrapped chain and not a
            # genuinely distant one.
            n = np.round(delta @ Linv)
            if np.any(n != 0):
                f.atoms.translate(-(n @ L))
        return ts


def attach_prep(u: mda.Universe, channel: mda.AtomGroup) -> None:
    """Frame preparation as on-the-fly transformations: once, before descriptors.

    Attached to the trajectory rather than applied by hand, so every frame a
    caller touches -- including a random-access jump to frame i -- is prepared
    the same way. MDAnalysis refuses a second attachment on the same trajectory,
    which is the guard against preparing a frame twice.
    """
    others = u.atoms - channel
    # The order is strict. unwrap first, because a channel split across the
    # boundary has no meaningful centre of mass; then pin the chains to one image,
    # because unwrap chooses each chain's image on its own; then centre the
    # channel, which is what puts the pore far from every cell face; and only then
    # wrap the rest back in. Wrapping by residue rather than by atom keeps each
    # water molecule whole, so a dipole is never split across the boundary and the
    # dipole cosine along the wire stays a physical quantity.
    workflow = [
        trf.unwrap(channel),
        _FragmentImageAligner(channel),
        trf.center_in_box(channel, center="mass"),
        trf.wrap(others, compound="residues"),
    ]
    u.trajectory.add_transformations(*workflow)


def frame_times_ps(u: mda.Universe) -> tuple[float, float]:
    """Frame step from the first two frame times in the file; missing timing is a refusal.

    Returns (t0, dt) in picoseconds. Every horizon, every backward window and
    every rate in the pipeline is expressed in time units, so a fabricated unit
    step would rescale all of them at once and the numbers would still look
    plausible. Hence the refusal: a file without usable timing is not read.
    """
    if u.trajectory.n_frames < 2:
        raise IoError("fewer than two frames — the time step cannot be measured")
    t0 = float(u.trajectory[0].time)
    t1 = float(u.trajectory[1].time)
    if not np.isfinite(t0) or not np.isfinite(t1) or t1 <= t0:
        raise IoError("frame timing is missing or non-monotonic — "
                      "refuse rather than assume a unit step")
    return t0, t1 - t0


def check_replica_sanity(u: mda.Universe, groups: dict[str, mda.AtomGroup],
                         horizon_ps: float, n_probe: int = 25) -> dict:
    """Cheap fitness checks, run before the main computation.

    Returns the timing of the replica: t0_ps, dt_ps, n_frames, span_ps. Each
    check refuses a file whose defect would otherwise be absorbed by the
    pipeline and reappear as a number.
    """
    t0, dt = frame_times_ps(u)
    n = u.trajectory.n_frames
    t_last = float(u.trajectory[-1].time)
    span = t_last - t0
    # The step measured on the first pair must agree with the average over the
    # whole file to one part in a thousand. A file assembled from segments
    # written at different intervals would put the frames on an uneven grid, and
    # the backward windows -- which are counted in time, not in frames -- would
    # cover a different number of frames in different parts of the trajectory.
    implied = span / (n - 1)
    if abs(implied - dt) > 1e-3 * dt:
        raise IoError(f"time step varies within the file: first {dt}, mean {implied}")
    # Frames whose forward window runs past the end of the file are censored, not
    # negative. Below a span of a few horizons the censored tail would be a large
    # fraction of the trajectory, so the demand is that the span clearly exceed
    # the horizon rather than merely exceed it.
    if span <= 3 * horizon_ps:
        raise IoError(f"trajectory span {span} ps does not clearly exceed the "
                      f"horizon {horizon_ps} ps: censoring would consume the sample")
    # Probe frames spread over the whole file. Each cell length is compared with
    # the first probe, so both a discontinuity (a restart against another box)
    # and a slow drift are caught. A cell that changes by that much invalidates
    # the fixed offsets the axial bins are defined on.
    idx = np.linspace(0, n - 1, n_probe, dtype=int)
    box0 = None
    for i in idx:
        ts = u.trajectory[int(i)]
        box = ts.dimensions[:3].copy()
        if box0 is None:
            box0 = box
        elif np.any(np.abs(box - box0) / box0 > 0.1):
            raise IoError(f"box changes abruptly at frame {i}: {box0} -> {box}")
        for gname, ag in groups.items():
            if len(ag) == 0:
                raise IoError(f"group {gname} is empty at frame {i}")
    return {"t0_ps": t0, "dt_ps": dt, "n_frames": n, "span_ps": span}


def triclinic_flag(u: mda.Universe) -> bool:
    """True when the cell is not rectangular; recorded by autodetect as a measured fact.

    Angles are in degrees, and the tolerance is only there to absorb the rounding
    of a nominally right angle in the file. The flag matters because a distance
    in a non-orthogonal cell cannot be minimum-imaged component by component --
    the gramicidin box is hexagonal -- so the artifact states which case the run
    was in rather than leaving it to be inferred.
    """
    ang = u.trajectory[0].dimensions[3:]
    return bool(np.any(np.abs(np.asarray(ang) - 90.0) > 1e-3))
