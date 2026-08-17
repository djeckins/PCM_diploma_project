"""Feature-table assembly, output writing, and the missing-value budget.

The schema is COMPUTED: the column list is a function of the enabled blocks ×
the detected number of sites × the number of bins × flags; no literal column
count exists anywhere. All descriptors are per-frame; downstream of this table
the trajectory no longer exists in the system.

Work is split in two: ReplicaComputer measures each frame once and the blocks read
those measurements, so a quantity has one estimator no matter how many columns
report it. run_step then fans the replicas out over processes, concatenates their
tables and writes features.parquet next to schema.json, applicability.json,
resolution.json and diagnostics.json. Every check that can fail the build lives in
run_step, so a table that exists on disk has already passed all of them.
"""

from __future__ import annotations

import os

import json
from pathlib import Path

import numpy as np
import pandas as pd
from MDAnalysis.lib.distances import capped_distance

from .. import io as pio
from ..autodetect import detect_filter, filter_oxygen_layers
from ..events import relative_axial
from ..baselines.published import rao_energy, rao_sigma_d
from ..config import Config
from ..config import load as load_config
from ..pore import PoreProfiler
from ..runtime import StepLog, run_dir, step_output
from ..vdw import resolve_radii
from . import (
    delivery,
    dynamics,
    electrostatics,
    fluctuations,
    geometry,
    hydration,
    named_sites,
    occupancy,
    symmetry,
    water_wire,
)
from ._common import ColSpec, FrameData, SchemaCtx, bin_edges, load_ww_scale

BLOCK_MODULES = {
    "geometry": geometry, "hydration": hydration, "water_wire": water_wire,
    "occupancy": occupancy, "named_sites": named_sites,
    "electrostatics": electrostatics, "symmetry": symmetry,
    "fluctuations": fluctuations, "dynamics": dynamics, "delivery": delivery,
}

# The published-criterion arm. These four columns carry to_model=False: they are the
# comparison arm of the evaluation, not inputs, and feeding them to the model would
# make the arm part of what it is supposed to be compared against. They live in the
# same table so that the arm and the model are scored on identical rows.
BASELINE_COLS = [
    ColSpec("bl_rao_E_kJmol", "baseline", "kJ/mol",
            "dewetting energy from the vendored Rao-2019 surface as a function "
            "of (lining hydrophobicity, constriction radius) — PER FRAME; a "
            "declared deviation from the paper (the authors define the inputs "
            "on time averages)",
            "inputs not measured on this frame", indicator="geo_n_search",
            conditional=True, to_model=False),
    ColSpec("bl_rao_E_win_kJmol", "baseline", "kJ/mol",
            "the same surface on inputs averaged over a backward window — the "
            "version comparable to the authors' averaging; its difference from "
            "the per-frame value measures how justified the transfer is",
            "no measured inputs in the window", indicator="geo_n_search",
            conditional=True, to_model=False),
    ColSpec("bl_rao_sigma_d", "baseline", "contour distance sum",
            "the published Rao-2019 STRUCTURE score: every pore-facing residue "
            "(CHAP rule, COG within local radius + 7.5 A margin) contributes a "
            "(kernel-smoothed hydrophobicity, local radius) point; points past "
            "the 1RT contour by the authors' nearest-node rule contribute their "
            "distance to it; the paper calls sigma_d > 0.55 non-conductive",
            "no pore-facing residue points measured on this frame",
            indicator="bl_rao_flag_n", conditional=True, to_model=False),
    ColSpec("bl_rao_flag_n", "baseline", "count",
            "number of residue points past the 1RT contour (indicator)",
            "no pore-facing residue points measured on this frame",
            conditional=True, to_model=False),
]

RAO_HYDROPHOBICITY_BANDWIDTH_A = 4.5  # CHAP's kernel bandwidth (0.45 nm)
RAO_NARROW_FILTER_NM = 0.7  # the authors keep residues with local radius < 0.7 nm


def _rao_sigma_frame(ctx: SchemaCtx, fd: FrameData) -> tuple[float, float] | None:
    """One frame's Rao structure score from the per-residue lining points.

    Returns (sigma_d, number of flagged points), or None when the frame offers no
    lining point the scale covers, in which case the two baseline columns stay missing.
    The score is a sum of distances in the authors' mixed (hydrophobicity, nm) metric,
    so it has no physical unit; SCHEMA.md calls it a contour distance sum.
    """
    if fd.rao_res_s_A is None or len(fd.rao_res_s_A) == 0:
        return None
    ww = ctx.params["ww_scale"]
    known = np.array([rn in ww for rn in fd.rao_res_resnames])
    if not known.any():
        return None
    s = fd.rao_res_s_A[known]
    # The vendored surface is tabulated in nanometres; our profile radii are
    # angstroms.
    # Getting this conversion wrong would move every point to a different part of the
    # surface and still return a plausible-looking score.
    r_nm = fd.rao_res_R_A[known] / 10.0
    raw_h = np.array([ww[rn] for rn, k in zip(fd.rao_res_resnames, known) if k])
    # CHAP's hydrophobicity profile: Nadaraya-Watson Gaussian smoothing along
    # the pore axis, evaluated at each residue's own position.
    d2 = (s[:, None] - s[None, :]) ** 2
    K = np.exp(-d2 / (2.0 * RAO_HYDROPHOBICITY_BANDWIDTH_A ** 2))
    h = (K @ raw_h) / K.sum(axis=1)
    # Only the narrow part of the pathway can gate, which is the authors' own filter;
    # keeping wide-lumen residues would add points that their contour was never drawn
    # for and inflate the score of an open channel.
    narrow = r_nm < RAO_NARROW_FILTER_NM
    if not narrow.any():
        return 0.0, 0.0  # nothing narrow enough to gate: wet by construction
    sigma, n_flag = rao_sigma_d(h[narrow], r_nm[narrow])
    return sigma, float(n_flag)


def build_ctx(cfg: Config) -> SchemaCtx:
    """Schema context for one config: the sizes that shape the column list, plus params.

    Reads only the config, so it can be called before a trajectory is opened. The params
    dict carries the estimators' numeric settings in the units their names state —
    angstroms for the pore geometry, picoseconds for the window — and the vendored
    hydrophobicity scale, read from its checksummed file instead of being written into
    the code.
    """
    edges = bin_edges(cfg["data.pore.bin_low_offset_A"], cfg["data.pore.bin_width_A"],
                      cfg["data.pore.n_bins"])
    return SchemaCtx(
        blocks=list(cfg["features.blocks"]),
        n_bins=cfg["data.pore.n_bins"],
        n_sites=cfg["system.arch_profile.n_sites"],
        sites_basis=cfg["system.arch_profile.sites_basis"],
        filter_present=bool(cfg["system.arch_profile.filter_present"]),
        n_subunits=cfg["system.arch_profile.n_subunits"],
        partner_cutoff_A=cfg["features.partner_cutoff_A"],
        window_ps=cfg["features.window_ps"],
        params={
            "bin_edges": edges,
            "ww_scale": load_ww_scale(),
            "coordination_cutoff_A": cfg["features.coordination_cutoff_A"],
            "site_sigma_A": cfg["features.site_sigma_A"],
            "site_centers_offset_A": cfg["features.site_centers_offset_A"],
            "mouth_radius_A": cfg["features.mouth_radius_A"],
            "mouth_depth_A": cfg["features.mouth_depth_A"],
            "lining_band_A": cfg["features.lining_band_A"],
        })


def build_schema(cfg: Config) -> list[ColSpec]:
    """The ordered column list for a config, with the three build-time invariants checked.

    The order is the config's block order followed by the baseline columns, and it is the
    order the table's columns are written in, so schema and table can be compared name by
    name later. The checks are on the schema alone and cost nothing, but each of them
    catches a defect that would otherwise surface as a silently wrong table: a duplicated
    name means one quantity has two estimators, a column with no block cannot be
    aggregated to a mechanism axis, and an indicator that names a column outside the
    schema leaves a gap that nothing can explain.
    """
    ctx = build_ctx(cfg)
    cols: list[ColSpec] = []
    for b in ctx.blocks:
        cols.extend(BLOCK_MODULES[b].schema(ctx))
    cols.extend(BASELINE_COLS)
    names = [c.name for c in cols]
    if len(names) != len(set(names)):
        dup = sorted({n for n in names if names.count(n) > 1})
        raise RuntimeError(f"two columns share one name: {dup} — two "
                           f"implementations of one quantity are forbidden")
    # Block classification must be complete: a column without a block fails the
    # build.
    for c in cols:
        if not c.block:
            raise RuntimeError(f"column {c.name} has no block")
    # The column → indicator map must close within the schema.
    for c in cols:
        if c.indicator is not None and c.indicator not in names:
            raise RuntimeError(f"{c.name}: indicator {c.indicator} does not exist in schema")
    return cols


class ReplicaComputer:
    """Per-frame shared measurements of one replica — one estimator per quantity.

    Holds one open trajectory and measures a frame into a FrameData, which the blocks
    then read. Everything that does not change from frame to frame is resolved once in
    the constructor: selections, atom radii, the water topology, the subunit list, the
    backbone carbonyl pairs. A frame measurement therefore does no selection parsing and
    no topology work, which is what makes a full trajectory affordable.
    """

    def __init__(self, cfg: Config, cond: dict, rep: dict):
        """Open the replica and resolve everything that is constant over its frames.

        cond and rep are the config's condition and replica records; cond may carry a
        conduction direction, and the entry-side columns exist only when it does.
        """
        self.cfg = cfg
        self.cond = cond
        self.rep = rep
        self.ctx = build_ctx(cfg)
        self.u = pio.open_replica(rep)
        self.groups = pio.resolve_groups(self.u, cfg)
        pio.attach_prep(self.u, self.groups["channel"])
        # ax is the pore axis, plane the two directions perpendicular to it, so that
        # radial distances are the same expression whichever axis a system declares.
        self.ax = pio.AXIS_INDEX[cfg["data.pore.axis"]]
        self.plane = [i for i in range(3) if i != self.ax]
        self.anchor_ag = self.u.select_atoms(cfg["data.selections.anchor"])
        # Two ion groups with different jobs: self.ions is the permeant species,
        # whose
        # positions define occupancy, sites and delivery; self.all_ions adds the
        # counter-ions, which carry charge and so enter the electrostatic sums even
        # though they never permeate.
        self.ions = self.u.select_atoms(f"resname {cfg['system.permeant']}")
        self.all_ions = self.groups["ion"] + self.groups["ion_negative"]
        self.cyl = cfg["data.pore.cylinder_radius_A"]
        self.edges = self.ctx.params["bin_edges"]
        # Coordinate-based applicability: the declaration is checked against the
        # structure.
        filt = detect_filter(self.u, cfg["data.selections.channel"], self.ax,
                     cfg["system.arch_profile.filter_motif_pattern"])
        declared = bool(cfg["system.arch_profile.filter_present"])
        if filt["filter_present"] != declared:
            raise RuntimeError(
                f"config declares filter_present={declared}, the structure says "
                f"{filt['filter_present']} ({filt['basis']}) — a config that silently "
                f"contradicts the coordinates fails the run")
        self.filter_layers = (filter_oxygen_layers(
            self.u, filt["motif_resids_by_copy"], cfg["data.selections.channel"], self.ax)
            if filt["filter_present"] else None)
        # Volume-occupying atoms for the profile: channel and lipids, with water
        # and ions
        # left out on purpose, so the measured radius is the lumen available to
        # the mobile
        # phase (see pore.py). Lipids belong in the list because an acyl chain
        # entering
        # the pore is a documented way for a channel to close.
        occl = self.groups["channel"] + self.groups["lipid"]
        self.occl = occl
        # Radii come from the sourced fallback ladder; the counter of fired
        # rungs travels
        # into the step's diagnostics, so the provenance of every radius stays
        # visible.
        self.occl_radii, self.vdw_steps = resolve_radii(
            occl.resnames, occl.names, occl.elements)
        z_lo = cfg["data.pore.z_low_offset_A"]
        z_hi = cfg["data.pore.z_high_offset_A"]
        self.profiler = PoreProfiler(
            z_offsets=np.arange(z_lo, z_hi + 1e-9, cfg["data.pore.z_step_A"]),
            search_radius_A=cfg["data.pore.probe_search_radius_A"],
            dr_A=cfg["data.pore.probe_grid_dr_A"],
            n_theta=cfg["data.pore.probe_grid_n_theta"],
            fence_max_gap_deg=cfg["data.pore.fence_max_gap_deg"],
            lipschitz_tol_A=cfg["data.pore.lipschitz_tol_A"],
            slab_pad_A=cfg["data.pore.slab_pad_A"])
        # Water: match each oxygen to its two hydrogens, once.
        # Both groups are sorted by residue index with a stable sort, so entry
        # 2k and 2k+1
        # of the hydrogens belong to oxygen k and the dipole of a molecule can
        # be formed
        # by reshaping rather than by searching. The count check is what makes that
        # pairing safe: a three-site water model with a dummy site, or a
        # topology where
        # the two groups disagree, would pair hydrogens with the wrong oxygens
        # silently.
        wox, wh = self.groups["water_oxygen"], self.groups["water_hydrogen"]
        o_order = np.argsort(wox.resindices, kind="stable")
        h_order = np.argsort(wh.resindices, kind="stable")
        if len(wh) != 2 * len(wox):
            raise RuntimeError("water is not O+2H — contradicts composition detection")
        self.wox_sorted = wox[o_order]
        self.wh_sorted = wh[h_order]
        # Subunits: connected fragments of the channel.
        # The size floor keeps the list to chains: a connected fragment of a few
        # atoms is
        # a ligand, an ion or a fragment of the topology, and counting it as a
        # subunit
        # would corrupt the symmetry columns.
        self.frags = [f for f in self.groups["channel"].fragments if len(f) >= 30]
        # Channel backbone carbonyls: per-residue C→O pairs, once.
        # A residue enters only if it holds exactly one atom named C and one
        # named O, so a
        # residue with a second carbonyl, or with none, drops out instead of
        # contributing
        # an ambiguous vector to the orientation average.
        ch = self.groups["channel"]
        c_idx, o_idx, res_of = [], [], []
        for res in ch.residues:
            cs = res.atoms.select_atoms("name C")
            os_ = res.atoms.select_atoms("name O")
            if len(cs) == 1 and len(os_) == 1:
                c_idx.append(cs[0].ix)
                o_idx.append(os_[0].ix)
                res_of.append(res.resindex)
        self.carbonyl_c = np.asarray(c_idx)
        self.carbonyl_o = np.asarray(o_idx)
        self.carbonyl_res = np.asarray(res_of)
        self.ch_resindices = ch.resindices
        self.ch_resnames_by_resindex = {r.resindex: r.resname for r in ch.residues}
        # Per-residue index compression for whole-lining COG computation.
        # The inverse map and the atom counts let all residue centres of geometry be
        # obtained from one bincount per Cartesian direction instead of a Python
        # loop over
        # residues, which matters because this runs on every frame of every
        # replica. The
        # sums are unweighted, so these are centres of geometry, not of mass.
        self._res_unique, self._res_inverse = np.unique(ch.resindices,
                                                        return_inverse=True)
        self._res_names = np.array([self.ch_resnames_by_resindex[r]
                                    for r in self._res_unique])
        self._res_counts = np.bincount(self._res_inverse).astype(float)
        # +1 or -1 when the condition declares a conduction direction, None when
        # it does
        # not; the entry-side columns and the sign canonicalization both depend
        # on it.
        self.direction = cond.get("direction")

    def measure_frame(self, ts) -> FrameData:
        """Everything measured from the coordinates of one frame, as a FrameData.

        ts is the current MDAnalysis timestep, already passed through the preparation
        transformations, so the channel is centred and the rest is wrapped by residue.
        The anchor's centre of mass is the origin of every axial coordinate written here,
        and the blocks read these arrays rather than the trajectory, which is why a
        neighbour search or a selection appears once per frame and not once per block.
        """
        cfg, ax, plane = self.cfg, self.ax, self.plane
        fd = FrameData()
        fd.time_ps = float(ts.time)
        anchor = self.anchor_ag.center_of_mass()
        fd.anchor = anchor
        prof = self.profiler.profile(self.occl.positions, self.occl_radii, anchor)
        fd.prof_R, fd.prof_search = prof["R"], prof["search"]
        fd.prof_boundary, fd.prof_z_offsets = prof["boundary"], prof["z_offsets"]
        # Axial position of the narrowest measured slice, published in
        # ctx.params for the
        # blocks that need the gate position without holding the profile. NaN
        # when the
        # profiler measured nothing inside the bin range, and the ion distance
        # and the
        # charge asymmetry then have nothing to be referenced against.
        interior = (fd.prof_z_offsets >= self.edges[0]) & (fd.prof_z_offsets <= self.edges[-1])
        s_int = fd.prof_search & interior
        z_constr = np.nan
        if np.any(s_int):
            z_constr = float(fd.prof_z_offsets[np.flatnonzero(s_int)[
                np.argmin(fd.prof_R[s_int])]])
        self.ctx.params["_z_constr_now"] = z_constr

        ipos = self.ions.positions
        # Minimum-image axial coordinates relative to the anchor: a change of
        # periodic
        # image of the channel shifts the whole frame, and only anchor-relative
        # coordinates are untouched by that shift.
        fd.ion_z_rel = relative_axial(ipos[:, ax], anchor[ax], ts.dimensions[ax])
        fd.ion_rxy = np.hypot(ipos[:, plane[0]] - anchor[plane[0]],
                              ipos[:, plane[1]] - anchor[plane[1]])
        # "In the pore" is radial confinement AND the axial range the
        # descriptors cover,
        # so an ion beside the protein at the right height does not count as inside.
        fd.ion_in_pore = ((fd.ion_rxy <= self.cyl)
                          & (fd.ion_z_rel >= self.edges[0])
                          & (fd.ion_z_rel <= self.edges[-1]))
        # Coordination of all permeant ions from a single distance map.
        # The cutoff is this system's first minimum of the ion-water g(r), and
        # the pair
        # search is periodic (box=), so an ion near a box face is not reported as
        # under-coordinated. Counting pair rows per ion gives the first-shell number.
        cut = self.ctx.params["coordination_cutoff_A"]
        pairs = capped_distance(ipos, self.wox_sorted.positions, max_cutoff=cut,
                                box=ts.dimensions, return_distances=False)
        coord_n = np.bincount(pairs[:, 0], minlength=len(self.ions)).astype(float)
        # Bulk reference for the desolvation column, taken in this same frame:
        # ions that
        # are out of the pore and well away from the axis, the radial margin
        # keeping ions
        # loitering at a mouth out of the reference. The median, not the mean,
        # so a single
        # oddly coordinated ion cannot move it.
        bulk_mask = (~fd.ion_in_pore) & (fd.ion_rxy > 2 * self.cyl)
        self.ctx.params["_bulk_coord_now"] = (float(np.median(coord_n[bulk_mask]))
                                              if np.any(bulk_mask) else np.nan)
        # The innermost ion is the one nearest the gate along the axis. It is a
        # composite
        # object: it needs both an ion in the pore and a measured constriction,
        # which is
        # why occ_has_innermost, not the ion count, is the indicator of its columns.
        fd.innermost = None
        if np.any(fd.ion_in_pore) and np.isfinite(z_constr):
            cand = np.flatnonzero(fd.ion_in_pore)
            fd.innermost = int(cand[np.argmin(np.abs(fd.ion_z_rel[cand] - z_constr))])
            fd.innermost_coord_n = coord_n[fd.innermost]

        # Water in the cylinder: axial coordinates and dipoles.
        # Radial confinement only, with no axial cut: the water blocks apply
        # their own
        # axial ranges, and both arrays below are stored in the same order so a mask
        # taken on one is valid on the other.
        wpos = self.wox_sorted.positions
        w_rxy = np.hypot(wpos[:, plane[0]] - anchor[plane[0]],
                         wpos[:, plane[1]] - anchor[plane[1]])
        w_z = relative_axial(wpos[:, ax], anchor[ax], ts.dimensions[ax])
        w_in = w_rxy <= self.cyl
        fd.wat_z_rel = w_z[w_in]
        # Dipole direction from the oxygen to the midpoint of its two hydrogens, i.e.
        # toward the positive end of the molecule, so cos = +1 means the dipole
        # moment
        # points along the axis. Only the direction is used; the magnitude,
        # which would
        # be the same for every rigid water anyway, divides out. The pairing works
        # because the hydrogens were sorted with their oxygens in the constructor.
        h = self.wh_sorted.positions.reshape(-1, 2, 3)
        dip = 0.5 * (h[:, 0] + h[:, 1]) - wpos
        # Guard against a zero-length vector; far below any real O-H geometry.
        norm = np.linalg.norm(dip, axis=1) + 1e-12
        fd.wat_dipole_cos = (dip[:, ax] / norm)[w_in]

        # Pore lining within the constriction band.
        fd.lining_resnames = []
        fd.lining_carbonyl_cos = np.nan
        fd.lining_n_carbonyls = 0
        if np.isfinite(z_constr):
            band = self.ctx.params["lining_band_A"]
            chpos = self.groups["channel"].positions
            ch_rxy = np.hypot(chpos[:, plane[0]] - anchor[plane[0]],
                              chpos[:, plane[1]] - anchor[plane[1]])
            # Plain difference, not minimum-image: the channel is centred in the
            # box by
            # the preparation, so its own atoms are all in the primary cell.
            ch_z = chpos[:, ax] - anchor[ax]
            # A residue counts as lining when it has an atom within the band
            # around the
            # constriction and close enough to the axis. The radial test is the
            # cylinder
            # widened by 2 A, because the atoms that FORM the wall sit just
            # outside the
            # cylinder that the mobile phase occupies; without the margin the
            # wall of a
            # narrow pore would be missed and the band would come out empty.
            facing = (ch_rxy <= self.cyl + 2.0) & (np.abs(ch_z - z_constr) <= band)
            # One entry per residue, so a large residue does not weight the
            # hydrophobicity
            # average by its atom count.
            res_idx = np.unique(self.ch_resindices[facing])
            fd.lining_resnames = [self.ch_resnames_by_resindex[r] for r in res_idx]

        # Rao heuristic: one point per pore-facing residue along the WHOLE pore.
        # Facing follows the CHAP rule — residue centre of geometry closer to
        # the axis than the local pathway radius plus the tool's default margin
        # (pm-pl-margin = 0.75 nm); the radius carried by the point is the
        # profile radius at the residue's axial position.
        s_ok = fd.prof_search
        if np.any(s_ok):
            chpos = self.groups["channel"].positions
            cogs = np.empty((len(self._res_unique), 3))
            for d in range(3):
                cogs[:, d] = (np.bincount(self._res_inverse, weights=chpos[:, d])
                              / self._res_counts)
            cog_s = cogs[:, ax] - anchor[ax]
            cog_rxy = np.hypot(cogs[:, plane[0]] - anchor[plane[0]],
                               cogs[:, plane[1]] - anchor[plane[1]])
            # Radius at a residue's own axial position, interpolated over the
            # measured
            # slices only. The interpolation bridges refused slices, which is
            # acceptable
            # here and only here: this radius feeds the published-criterion arm,
            # whose
            # rule needs a value at every lining residue, while the geometry columns
            # continue to export a refusal as a gap.
            zs, Rs = fd.prof_z_offsets[s_ok], fd.prof_R[s_ok]
            inside = (cog_s >= zs[0]) & (cog_s <= zs[-1])
            R_at = np.full(len(cog_s), np.nan)
            R_at[inside] = np.interp(cog_s[inside], zs, Rs)
            facing_pore = inside & (cog_rxy <= R_at + 7.5)
            fd.rao_res_resnames = list(self._res_names[facing_pore])
            fd.rao_res_s_A = cog_s[facing_pore]
            fd.rao_res_R_A = R_at[facing_pore]
            # Backbone C=O orientation of the residues in the constriction band,
            # as the
            # cosine of the C→O vector against the axis. This is the observable
            # of the
            # rearrangement route: a carbonyl that turns away from the lumen
            # removes the
            # coordination it was providing while leaving the radius as it was.
            if len(self.carbonyl_c) and len(res_idx):
                in_band = np.isin(self.carbonyl_res, res_idx)
                if np.any(in_band):
                    cpos = self.u.atoms.positions[self.carbonyl_c[in_band]]
                    opos = self.u.atoms.positions[self.carbonyl_o[in_band]]
                    v = opos - cpos
                    cos = v[:, ax] / (np.linalg.norm(v, axis=1) + 1e-12)
                    fd.lining_carbonyl_cos = float(np.mean(cos))
                    fd.lining_n_carbonyls = int(np.sum(in_band))

        # Sites: filter rings per frame, or static density peaks.
        # For a filter the ring positions are re-read every frame and a site is the
        # midpoint between adjacent rings, so the sites follow the filter as it
        # breathes
        # instead of standing at a fixed offset. n_rings rings give n_rings-1
        # sites, which
        # is where the site count of the architecture profile comes from.
        if self.filter_layers is not None:
            layer_z = np.array([float(g.positions[:, ax].mean()) for g in self.filter_layers])
            layer_z = np.sort(layer_z) - anchor[ax]
            fd.site_centers_z = 0.5 * (layer_z[1:] + layer_z[:-1])
        elif self.ctx.params["site_centers_offset_A"]:
            # Without a filter the centres are declared in the config as fixed
            # offsets
            # from the anchor: they discretize an ion-density profile and are
            # not sites in
            # the structural sense, which is why they only ever feed the soft
            # columns.
            fd.site_centers_z = np.asarray(self.ctx.params["site_centers_offset_A"], float)

        # Subunits.
        # Centres of mass relative to the anchor, so the symmetry columns see the
        # arrangement of the chains and not the position of the assembly in the box.
        if len(self.frags) >= 2:
            fd.subunit_com = np.array([f.center_of_mass() - anchor for f in self.frags])

        # Electrostatics around the innermost ion.
        # Sum of q_j / r_ij over each partner class within the cutoff, in e/A: the
        # potential at the ion up to the factor 1/(4 pi eps0), which is a
        # constant and is
        # deliberately not applied, since the columns are compared with each
        # other and
        # with themselves across frames rather than converted to an energy.
        if fd.innermost is not None:
            target = ipos[fd.innermost][None, :]
            pcut = self.ctx.params.get("partner_cutoff_A", self.ctx.partner_cutoff_A)
            fd.pot_by_partner = {}
            # Each class is summed on its own, so the four columns decompose one
            # quantity
            # by source. Water enters with both its oxygens and its hydrogens:
            # the molecule
            # is neutral overall, and only the separation of its charges
            # produces a field
            # at the ion, so dropping the hydrogens would leave a spurious net
            # charge.
            partners = {
                "protein": self.groups["channel"],
                "water": self.wox_sorted + self.wh_sorted,
                "ions": self.all_ions,
                "lipid": self.groups["lipid"],
            }
            for name, ag in partners.items():
                pr, dist = capped_distance(target, ag.positions, max_cutoff=pcut,
                                           box=ts.dimensions, return_distances=True)
                # The ion is a member of its own class: it must be removed by
                # global atom
                # index, or its own charge at r = 0 would dominate the sum.
                if name == "ions":
                    keep = ag.ix[pr[:, 1]] != self.ions.ix[fd.innermost]
                    pr, dist = pr[keep], dist[keep]
                q = ag.charges[pr[:, 1]]
                # Charges are taken from the topology as they are; any ECC-style
                # scaling
                # belongs to the replica's own protocol and is not undone here.
                # The floor
                # on the distance only protects the division.
                fd.pot_by_partner[name] = float(np.sum(q / np.maximum(dist, 1e-6)))
        # Axial charge asymmetry of mobile particles in the cylinder.
        # Net charge above the constriction minus below, in units of e: the sign
        # function
        # turns the sum into a difference between the two sides. Only mobile
        # particles
        # (ions and water) are counted, so this is the redistribution the field
        # acts on
        # rather than the fixed charge of the protein.
        if np.isfinite(z_constr):
            asym = 0.0
            for ag, rxy, zz in (
                    (self.all_ions, None, None),
                    (self.wox_sorted + self.wh_sorted, None, None)):
                pos = ag.positions
                r = np.hypot(pos[:, plane[0]] - anchor[plane[0]],
                             pos[:, plane[1]] - anchor[plane[1]])
                z = pos[:, ax] - anchor[ax]
                m = (r <= self.cyl) & (z >= self.edges[0]) & (z <= self.edges[-1])
                asym += float(np.sum(ag.charges[m] * np.sign(z[m] - z_constr)))
            fd.axial_charge_asym_e = asym

        # Delivery: pore mouths.
        # Two disc-shaped zones just outside the ends of the descriptor range,
        # of radius
        # mouth_radius_A and axial depth mouth_depth_A. They are wider than the pore
        # cylinder on purpose: an ion about to enter is still in the vestibule,
        # not yet in
        # the lumen, and counting it only once inside would make delivery invisible.
        mr = self.ctx.params["mouth_radius_A"]
        md = self.ctx.params["mouth_depth_A"]
        top, bot = self.edges[-1], self.edges[0]
        up_n = int(np.sum((fd.ion_rxy <= mr) & (fd.ion_z_rel > top)
                          & (fd.ion_z_rel <= top + md)))
        lo_n = int(np.sum((fd.ion_rxy <= mr) & (fd.ion_z_rel < bot)
                          & (fd.ion_z_rel >= bot - md)))
        # upper and lower stay laboratory sides; the entry zone below is the one
        # quantity
        # in this frame that reads the conduction direction.
        fd.mouth_counts = {"upper": up_n, "lower": lo_n}
        if self.direction is not None:
            # direction = +1 means conduction toward increasing anchor-relative axial
            # coordinate, the same convention the events step uses, so ions
            # enter at the
            # lower edge; direction = -1 reverses both.
            entry_edge = bot if self.direction > 0 else top
            outside = (fd.ion_z_rel < entry_edge) if self.direction > 0 else (fd.ion_z_rel > entry_edge)
            fd.mouth_counts["entry"] = lo_n if self.direction > 0 else up_n
            # Centre of the entry mouth on the axis, in laboratory coordinates,
            # since the
            # distance below is a full three-dimensional one: an ion off to the
            # side is
            # farther from the entry than an ion straight above it at the same
            # height.
            mouth_center = anchor.copy()
            mouth_center[ax] = anchor[ax] + entry_edge
            cand = np.flatnonzero(outside)
            if len(cand):
                d = np.linalg.norm(ipos[cand] - mouth_center[None, :], axis=1)
                # The wide zone, three mouth radii, exists to make the absence of a
                # candidate a measured fact: dlv_dist_entry_A is missing when
                # the zone is
                # empty, and dlv_entry_wide_n is the companion that says so.
                # Without it a
                # frame with no ion anywhere near the entry and a frame whose
                # distance
                # simply was not measured would look the same.
                wide = d <= 3 * mr
                fd.entry_wide_n = int(np.sum(wide))
                if np.any(wide):
                    fd.entry_dist_A = float(np.min(d[wide]))
        return fd


def compute_replica(cfg: Config, cond: dict, rep: dict,
                    exit_times_dir: np.ndarray) -> tuple[pd.DataFrame, dict]:
    """One replica's feature table plus its diagnostics.

    exit_times_dir holds the completion times in picoseconds of this replica's crossings
    in the condition's conduction direction; the delivery block turns them into the
    history columns. Returns the table with condition, replica and time_ps in front of the
    schema's columns, in schema order, and a diagnostics dict recording the van der Waals
    ladder, the columns that were sign-flipped and the direction that flipped them.
    """
    rc = ReplicaComputer(cfg, cond, rep)
    ctx = rc.ctx
    ctx.params["exit_times_dir"] = exit_times_dir
    ctx.params["direction"] = rc.direction
    schema = build_schema(cfg)
    stride = cfg["data.stride"]
    frames = list(range(0, rc.u.trajectory.n_frames, stride))
    T = len(frames)
    # Preallocated as NaN, never as zero: a column the blocks do not write on a
    # frame must
    # come out missing, because zero is a measurement and absence is not.
    cols = {c.name: np.full(T, np.nan) for c in schema}
    times = np.empty(T)
    # Per-ion histories for the dynamics post step. The second index is the
    # permeant ion,
    # in a fixed order for the whole replica, so a column follows one ion
    # through time.
    N = len(rc.ions)
    ion_z_all = np.full((T, N), np.nan)
    ion_in_pore = np.zeros((T, N), dtype=bool)
    innermost = np.full(T, -1, dtype=int)
    active_blocks = [BLOCK_MODULES[b] for b in ctx.blocks]
    for k, fi in enumerate(frames):
        ts = rc.u.trajectory[fi]
        fd = rc.measure_frame(ts)
        times[k] = fd.time_ps
        ion_z_all[k] = fd.ion_z_rel
        ion_in_pore[k] = fd.ion_in_pore
        innermost[k] = fd.innermost if fd.innermost is not None else -1
        for mod in active_blocks:
            if hasattr(mod, "compute"):
                mod.compute(ctx, fd, k, cols)
        # The published structure score is computed here rather than in a block:
        # it is a
        # comparison arm reading the same frame, not a descriptor family of this
        # work.
        sig = _rao_sigma_frame(ctx, fd)
        if sig is not None:
            cols["bl_rao_sigma_d"][k], cols["bl_rao_flag_n"][k] = sig
    # The post steps run only now, because a windowed column needs the whole
    # trajectory in
    # hand; they are also the only place the per-ion histories are read.
    ctx.params["_aux_ion_z_all"] = ion_z_all
    ctx.params["_aux_ion_in_pore"] = ion_in_pore
    ctx.params["_aux_innermost"] = innermost
    for mod in active_blocks:
        if hasattr(mod, "post"):
            mod.post(ctx, times, cols)
    # Baseline columns: excluded from the model, kept in the table.
    # The authors' surface takes the radius in nanometres, hence the division by 10.
    r_nm = cols["geo_r_constriction_A"] / 10.0
    h_in = cols["hyd_lining_hydrophobicity"]
    cols["bl_rao_E_kJmol"] = rao_energy(h_in, r_nm)
    from ._common import backward_window_slice
    win = ctx.window_ps
    for i in range(T):
        sl = backward_window_slice(times, i, win)
        rr, hh = r_nm[sl], h_in[sl]
        # The inputs are averaged and the surface is then read once, not the
        # other way
        # round: the paper's inputs are time averages, and averaging the
        # energies of a
        # non-linear surface would give a different number that no published
        # rule defines.
        # Frames where either input is missing drop out, so the average is over
        # measured
        # frames alone.
        m = np.isfinite(rr) & np.isfinite(hh)
        if np.any(m):
            cols["bl_rao_E_win_kJmol"][i] = float(
                rao_energy(np.array([hh[m].mean()]), np.array([rr[m].mean()]))[0])
    # Sign canonicalization: exactly once, before the table is written.
    # Only the columns the schema marks as direction-bearing, and only for a
    # condition
    # conducting toward decreasing axial coordinate, so that a drift toward the
    # exit is
    # positive in every condition. Applying this twice, or to a structural
    # column, would
    # silently mirror the system.
    if rc.direction is not None and rc.direction < 0:
        for c in schema:
            if c.sign_flip:
                cols[c.name] = -cols[c.name]
    df = pd.DataFrame({"condition": cond["id"], "replica": rep["id"], "time_ps": times})
    for c in schema:
        df[c.name] = cols[c.name]
    # flipped_columns lists every direction-bearing column of the schema; whether the
    # flip actually fired is told by direction, which is recorded beside it.
    diag = {"vdw_ladder_steps": rc.vdw_steps,
            "flipped_columns": [c.name for c in schema if c.sign_flip],
            "direction": rc.direction}
    return df, diag


def _worker(cfg_path: str, cond_id: str, rep_id: str) -> tuple[pd.DataFrame, dict]:
    """One (condition, replica) job, run in its own process.

    Takes the config by path and re-reads it, and reads the events artifact itself, so
    that nothing but plain strings has to cross the process boundary and each worker
    opens exactly one trajectory.
    """
    cfg = load_config(cfg_path)
    cond = next(c for c in cfg["data.conditions"] if c["id"] == cond_id)
    rep = next(r for r in cond["replicas"] if r["id"] == rep_id)
    ev_path = run_dir(cfg) / "events" / "events.parquet"
    ev = pd.read_parquet(ev_path)
    sel = ev[(ev["condition"] == cond_id) & (ev["replica"] == rep_id)]
    direction = cond.get("direction")
    # Events of the condition's own direction, plus those marked 0. Direction 0
    # means an
    # event read from a provided annotation, which records times only and cannot
    # state a
    # direction; dropping them would discard the annotated crossings of exactly the
    # systems where an external annotation exists.
    if direction is not None:
        sel = sel[(sel["direction"] == direction) | (sel["direction"] == 0)]
    return compute_replica(cfg, cond, rep, sel["t_exit_ps"].to_numpy())


def resolve_applicability(cfg: Config, table: pd.DataFrame,
                          schema: list[ColSpec]) -> dict[str, str]:
    """active / structurally inapplicable / measured and constant.

    The structural part is resolved FROM COORDINATES (the structure check runs
    in compute_replica and fails the run on any conflict with the config);
    here its results are read from the config, whose agreement with the
    structure was established there.

    The three verdicts are different statements and the missing-value budget treats
    them differently. "inapplicable_structural" means the subject of the measurement
    does not exist in this system, so the column is empty by construction and must not
    count against the budget: S3 of a three-site filter, or the subunit spread of a
    single-chain channel. "constant" means the subject was there and was measured and
    did not change over the run, which is information about the system rather than a
    defect. Only "active" columns are held to the budget, and only they can fail it.
    """
    n_sites = cfg["system.arch_profile.n_sites"] or 0
    filter_present = bool(cfg["system.arch_profile.filter_present"])
    n_sub = cfg["system.arch_profile.n_subunits"] or 0
    verdict: dict[str, str] = {}
    for c in schema:
        v = "active"
        if c.block == "named_sites":
            if not filter_present:
                v = "inapplicable_structural"
            elif c.name.startswith("ns_S"):
                # Site index out of the column name, ns_S<k>_occ: the union
                # schema always
                # carries NS_UNION_SITES of them, and a filter with fewer rings
                # leaves the
                # tail inapplicable rather than absent.
                s = int(c.name[4:c.name.index("_", 4)])
                if s >= n_sites:
                    v = "inapplicable_structural"
        if c.block == "symmetry":
            # A spread needs two subunits; the ring of circular neighbours needs
            # three,
            # since with two the two distances of the ring are the same distance
            # twice.
            if n_sub < 2 or (c.name == "sym_nn_dist_cv" and n_sub < 3):
                v = "inapplicable_structural"
        if c.block == "occupancy" and c.name.startswith("occ_site") and n_sites == 0:
            v = "inapplicable_structural"
        if v == "active":
            vals = table[c.name].dropna()
            if len(vals) and vals.nunique() <= 1:
                v = "constant"  # subject present, measured, unchanged over the run
            elif len(vals) == 0 and not c.conditional:
                pass  # left to the missing-value budget
        verdict[c.name] = v
    return verdict


def run_step(cfg: Config) -> None:
    """Build the run's feature table and write the features step's artifacts.

    One job per (condition, replica), each in its own process because each opens a
    trajectory of its own; PCM2_MAX_WORKERS caps the pool. The checks that follow the
    concatenation are the reason this step can fail: schema against table, applicability,
    the missing-value budget, and the requirement that every gap be explained by an
    indicator. Writes features.parquet, schema.json, applicability.json, resolution.json
    and diagnostics.json, and returns nothing — later steps read the artifacts.
    """
    from concurrent.futures import ProcessPoolExecutor

    schema = build_schema(cfg)
    with step_output(cfg, "features") as out:
        log = StepLog(out)
        jobs = [(c["id"], r["id"]) for c in cfg["data.conditions"] for r in c["replicas"]]
        dfs, diags = [], {}
        with ProcessPoolExecutor(max_workers=min(len(jobs), int(os.environ.get("PCM2_MAX_WORKERS", "14")))) as ex:
            futs = {ex.submit(_worker, str(cfg.source_path), ci, ri): (ci, ri)
                    for ci, ri in jobs}
            for fut, key in futs.items():
                df, diag = fut.result()
                dfs.append(df)
                diags["/".join(key)] = diag
                log.say(f"[{key[0]}/{key[1]}] {len(df)} frames, "
                        f"vdw ladder: {diag['vdw_ladder_steps']}")
        table = pd.concat(dfs, ignore_index=True)

        # Self-consistency check inside the build: schema == table.
        # Compared as ordered lists, not as sets, so a column that moved is
        # caught as well
        # as one that went missing; the three identity columns are not part of
        # the schema.
        schema_names = [c.name for c in schema]
        table_cols = [c for c in table.columns if c not in ("condition", "replica", "time_ps")]
        if schema_names != table_cols:
            raise RuntimeError("schema and table diverged: "
                               f"{set(schema_names) ^ set(table_cols)}")

        verdict = resolve_applicability(cfg, table, schema)
        n_active = sum(1 for v in verdict.values() if v == "active")
        n_inap = sum(1 for v in verdict.values() if v == "inapplicable_structural")
        n_const = sum(1 for v in verdict.values() if v == "constant")
        assert n_active + n_inap + n_const == len(schema), "verdict counts do not add up"
        log.say(f"applicability: {len(schema)} columns total = {n_active} active "
                f"+ {n_inap} structurally inapplicable + {n_const} constant")
        for name, v in verdict.items():
            if v != "active":
                log.say(f"  {name}: {v}")

        # Missing-value budget AFTER masking and only over active columns;
        # conditional
        # columns are exempt by name, but "missing on all frames" still applies.
        # A conditional column is missing whenever its subject is absent, and
        # how often
        # that happens is a property of the trajectory rather than of the code, so a
        # threshold on it would only measure the system. Missing on every frame
        # is still
        # a failure: it means the subject was never once present, and the column
        # is then
        # carrying nothing at all.
        budget = cfg["accept.missing_budget_frac"]
        cond_names = {c.name for c in schema if c.conditional}
        problems = []
        for c in schema:
            if verdict[c.name] != "active":
                continue
            frac = float(table[c.name].isna().mean())
            if c.name in cond_names:
                if frac >= 1.0:
                    problems.append(f"{c.name}: conditional column missing on all frames")
            elif frac > budget:
                problems.append(f"{c.name}: missing fraction {frac:.2f} > budget {budget}")
        if problems:
            raise RuntimeError("missing-value budget: " + "; ".join(problems))

        # Indicators must explain gaps: subject present but value absent is a
        # failure.
        # The column → indicator → THRESHOLD map is verified on every rebuild.
        # Non-model (baseline) columns combine two estimators' inputs and are not
        # seen by the model; their gaps are explained by the input indicators.
        for c in schema:
            if c.indicator is None or verdict[c.name] != "active" or not c.to_model:
                continue
            # A frame is bad when the value is absent although its indicator says the
            # subject was there. Frames where the indicator itself is missing are not
            # counted, since they make no claim about the subject.
            miss = table[c.name].isna()
            present = table[c.indicator] >= c.indicator_min
            bad = miss & present & table[c.indicator].notna()
            frac_bad = float(bad.mean())
            # A small tolerance rather than zero: an indicator and its value are
            # measured
            # from the same frame but not by the same expression, so a boundary
            # case can
            # fall on either side of a comparison. A systematic defect is far
            # above this.
            if frac_bad > 0.001:
                raise RuntimeError(
                    f"{c.name}: missing while indicator {c.indicator}>={c.indicator_min} "
                    f"on {frac_bad:.1%} of frames — a build failure, not a gap")

        # Resolving power of each column, reported as a table.
        # A column can pass every check above and still be useless: if it takes
        # two values
        # over a whole run it cannot describe the approach to a crossing. This
        # is reported,
        # not enforced, because low resolution is a fact about the descriptor and the
        # system together. frac_changed is taken over the concatenated table, so one
        # comparison per replica boundary reaches across into the previous
        # trajectory.
        resolution = {}
        for c in schema:
            v = table[c.name]
            resolution[c.name] = {
                "n_distinct": int(v.nunique(dropna=True)),
                "frac_changed": float((v.diff() != 0).mean()),
                "missing_frac": float(v.isna().mean()),
            }
        (out / "resolution.json").write_text(json.dumps(resolution, indent=1))
        # The two exempt columns are binary by definition, so two distinct
        # values is what
        # they are supposed to have and reporting them would train the reader to
        # ignore
        # this line.
        low = [n for n, r in resolution.items()
               if r["n_distinct"] <= 2 and verdict[n] == "active"
               and n not in ("ww_continuous", "dlv_crossed_before")]
        if low:
            log.say(f"columns with <=2 distinct values (do not resolve dynamics): {low}")

        table.to_parquet(out / "features.parquet")
        (out / "schema.json").write_text(json.dumps(
            {"columns": [c.as_dict() for c in schema],
             "n_columns_is_derived": "len(columns); a literal count is forbidden"},
            indent=1, ensure_ascii=False))
        (out / "applicability.json").write_text(json.dumps(verdict, indent=1, ensure_ascii=False))
        (out / "diagnostics.json").write_text(json.dumps(diags, indent=1, ensure_ascii=False, default=str))
        log.say(f"table: {len(table)} rows x {len(schema)} columns")
        log.close()


def load_features(cfg: Config) -> tuple[pd.DataFrame, list[ColSpec], dict]:
    """Read a finished features step back as (table, schema, applicability verdicts).

    The single entry point for every later step, so no consumer rebuilds the schema from
    the config: what is read is what was written, including the column order. The schema
    comes back as ColSpec objects, which is what lets a model step know which columns are
    conditional, which carry a sign and which are not meant for the model at all.
    """
    root = run_dir(cfg) / "features"
    table = pd.read_parquet(root / "features.parquet")
    doc = json.loads((root / "schema.json").read_text())
    schema = [ColSpec(**d) for d in doc["columns"]]
    verdict = json.loads((root / "applicability.json").read_text())
    return table, schema, verdict
