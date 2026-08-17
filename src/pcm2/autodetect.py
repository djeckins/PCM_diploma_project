"""Autodetection: every parameter is measured from the input files, or the step refuses.

Each measured quantity carries a basis: what was measured, over how many frames,
by what margin over the alternative. A refusal names the config key by which a
human settles the question. Measured values are appended to the config file the
run was launched with.
"""

from __future__ import annotations

import json
import re
from collections import Counter

import numpy as np

from . import io as pio
from .config import Config, save
from .events import ladder_pass, relative_axial, zone_planes
from .pore import PoreProfiler
from .runtime import StepLog, step_output
from .vdw import resolve_radii

# Charge separating an ion from a neutral particle: even under ECC scaling
# ~0.7 the charge of a monovalent ion stays well above the partial charges of
# virtual sites; the threshold is half the smallest ECC scale encountered.
ION_MIN_ABS_CHARGE = 0.35

ONE_LETTER = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "HSD": "H", "HSE": "H", "HSP": "H",
    "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # D-amino acids and modifications (gramicidin): same side chain, same letter.
    "DVAL": "V", "DLEU": "L", "FVA": "V",
}
# Canonical K-channel signature; engineered filter mutants may declare an
# extended pattern in the config (geometric ring confirmation still applies).
DEFAULT_FILTER_MOTIF = "T[VIL]G[YF]G"


class Refusal:
    """A measurement that did not decide, named by the config key that would settle it.

    key is the dotted config path a human can declare; message is the measurement that
    failed, with its numbers, so the declaration is an informed one. A refusal is a value
    and not an exception: one stage can report several at once, and report.json names all
    of them alongside the origins of everything that was measured.
    """

    def __init__(self, key: str, message: str):
        self.key = key
        self.message = message

    def as_dict(self):
        return {"key": self.key, "message": self.message}


def detect_selections(u) -> tuple[dict[str, str], dict[str, str], list[Refusal]]:
    """Group residues by COMPOSITION, not by residue names.

    Residue names are a property of the force-field setup (SOL, TIP3, TIP3P all name the
    same water), so a name-based rule would have to be extended for every new input. What
    a water, an ion and a lipid are made of does not change. → (selections, bases,
    refusals): the MDAnalysis selection strings, the text recorded as their origin, and
    the roles that could not be filled.
    """
    sels: dict[str, str] = {}
    bases: dict[str, str] = {}
    refusals: list[Refusal] = []
    res_kinds: dict[str, dict] = {}
    for res in u.residues:
        rn = res.resname
        # Composition is read from the first residue carrying a name; later
        # residues of
        # that name only raise its count. Within one topology a name means one
        # chemistry,
        # so the stoichiometry and the total charge do not need re-reading.
        if rn in res_kinds:
            res_kinds[rn]["count"] += 1
            continue
        atoms = res.atoms
        elems = Counter(str(e).upper() for e in atoms.elements)
        res_kinds[rn] = {"count": 1, "n_atoms": len(atoms),
                         "charge": float(atoms.charges.sum()),
                         "elements": elems}
    ion_pos, ion_neg, water, lipid, rest = [], [], [], [], []
    for rn, k in res_kinds.items():
        # One atom carrying most of an elementary charge is a monovalent ion;
        # the sign of
        # the charge separates the cations from the anions.
        if k["n_atoms"] == 1 and abs(k["charge"]) >= ION_MIN_ABS_CHARGE:
            (ion_pos if k["charge"] > 0 else ion_neg).append(rn)
        # Three atoms with the stoichiometry O + 2H: a rigid three-site water
        # model. A four-
        # or five-site model carries virtual sites and would not match, which is
        # what the
        # count of matched waters in the basis text below is there to expose.
        elif (k["n_atoms"] == 3 and k["elements"].get("O") == 1
              and k["elements"].get("H") == 2):
            water.append(rn)
        # A lipid is large AND numerous: a phospholipid runs to well over a
        # hundred atoms
        # and a bilayer holds tens to hundreds of copies of each species. An
        # amino acid
        # stays under 50 atoms, so no part of the protein can match, and a
        # cofactor or a
        # detergent present in a handful of copies fails the count.
        elif k["n_atoms"] > 50 and k["count"] > 20:
            lipid.append(rn)
        # Everything else is the channel, by subtraction: an unrecognized
        # residue stays
        # with the protein, so it still occupies volume in the pore profile
        # instead of
        # being invisible to every selection.
        else:
            rest.append(rn)
    if not ion_pos:
        refusals.append(Refusal("data.selections.ion",
                                "no single-atom residues with a substantial positive charge"))
    if not water:
        refusals.append(Refusal("data.selections.water_oxygen",
                                "no residue with O+2H stoichiometry found"))
    if not lipid:
        refusals.append(Refusal("data.selections.lipid",
                                "no large residue repeated many times found"))
    if not rest:
        refusals.append(Refusal("data.selections.channel",
                                "nothing left after subtracting water/ions/lipids"))
    if refusals:
        return sels, bases, refusals

    roles = {"ion_pos": sorted(ion_pos), "ion_neg": sorted(ion_neg),
             "water": sorted(water), "lipid": sorted(lipid), "rest": sorted(rest)}
    sels = selections_from_roles(roles)
    counts = {rn: res_kinds[rn]["count"] for rn in res_kinds}
    bases["selections"] = (
        f"by composition: cations {ion_pos} "
        f"(|q|={[round(res_kinds[r]['charge'],2) for r in ion_pos]}), "
        f"anions {ion_neg}, water {water} (O+2H, n={sum(counts[r] for r in water)}), "
        f"lipid {lipid} ({[res_kinds[r]['n_atoms'] for r in lipid]} atoms × "
        f"{[counts[r] for r in lipid]}), channel by subtraction: {sorted(rest)[:8]}…")
    bases["roles"] = roles
    return sels, bases, refusals


def selections_from_roles(roles: dict[str, list[str]]) -> dict[str, str]:
    """Turn the composition roles into the MDAnalysis selection strings of the config.

    Names are sorted inside every string so that the same roles always produce the same
    text: the selections are written into the config and compared between runs.

    Two of the strings are worth reading twice. ion_negative falls back to a selection
    that parses and matches nothing when the system has no anions, so the key is always a
    selection rather than sometimes null. channel is written as a subtraction, which is
    what puts anything unclassified on the protein's side.
    """
    def rl(names):
        return " ".join(sorted(names))
    ion_pos, ion_neg = roles["ion_pos"], roles["ion_neg"]
    water, lipid = roles["water"], roles["lipid"]
    sels = {
        "ion": f"resname {rl(ion_pos)}",
        "ion_negative": (f"resname {rl(ion_neg)}" if ion_neg
                         else f"resname {rl(ion_pos)} and name NONE"),
        "water_oxygen": f"resname {rl(water)} and element O",
        "water_hydrogen": f"resname {rl(water)} and element H",
        "lipid": f"resname {rl(lipid)}",
        # The headgroup marker is the phosphorus of the lipids: one atom per
        # lipid, at a
        # fixed place in the molecule, and it is what the leaflet planes are
        # measured from.
        "phosphate": f"resname {rl(lipid)} and element P",
        "channel": f"not resname {rl(ion_pos + ion_neg + water + lipid)}",
    }
    return sels


def merge_roles(per_topology: list[dict[str, list[str]]]) -> tuple[dict[str, list[str]], str | None]:
    """Union the composition classification across topologies.

    Residue names vary between force-field setups (SOL vs TIP3 water), so the
    invariant is the composition class. A residue name assigned to different
    roles in different topologies is a genuine composition conflict and is
    returned as an error message.

    → (union of roles, None) on success, or ({}, message) on a conflict. The union means
    one set of selections works on every replica of the run, so the descriptor tables of
    different conditions describe the same groups.
    """
    union: dict[str, set[str]] = {k: set() for k in ("ion_pos", "ion_neg", "water",
                                                     "lipid", "rest")}
    assigned: dict[str, str] = {}
    for roles in per_topology:
        for role, names in roles.items():
            for rn in names:
                if rn in assigned and assigned[rn] != role:
                    return {}, (f"residue {rn} is classified as {assigned[rn]} in one "
                                f"topology and as {role} in another")
                assigned[rn] = role
                union[role].add(rn)
    return {k: sorted(v) for k, v in union.items()}, None


def detect_axis(u, phosphate, n_frames_sample: int) -> tuple[str | None, str]:
    """Membrane normal from leaflet separation; an analytic floor is mandatory.

    The estimator is the emptiness of the central third between leaflets: along
    the normal no phosphates sit in the hydrophobic core at all, while along the
    in-plane axes about a third of them do.
    Score = 1/(phosphate fraction in the central third + 0.01): normal ≫ plane.

    → (axis letter, basis text) or (None, why it could not be decided). The sampled frames
    are spread over the whole trajectory rather than taken from its start, so no single
    unrepresentative stretch of the run decides the axis.
    """
    idx = np.linspace(0, u.trajectory.n_frames - 1, n_frames_sample, dtype=int)
    scores = {a: [] for a in "xyz"}
    for fi in idx:
        u.trajectory[int(fi)]
        pos = phosphate.positions
        for a, col in zip("xyz", range(3)):
            v = pos[:, col]
            # Split the phosphates at their own median and take the mean of each
            # half:
            # along the normal these are the two headgroup planes, along an
            # in-plane axis
            # they are two arbitrary halves of one flat sheet.
            med = np.median(v)
            m_lo, m_hi = v[v < med].mean(), v[v >= med].mean()
            band = (m_hi - m_lo) / 3.0
            central = np.mean((v > m_lo + band) & (v < m_hi - band))
            # +0.01 keeps the score finite; a completely empty central third
            # therefore
            # scores exactly 100, which is the ceiling of this estimator.
            scores[a].append(1.0 / (float(central) + 0.01))
    mean_scores = {a: float(np.mean(s)) for a, s in scores.items()}
    ranked = sorted(mean_scores.items(), key=lambda kv: -kv[1])
    best, second = ranked[0], ranked[1]
    # No-information floor: the best axis must clearly separate from the runner-up.
    # The absolute floor of 10 means at most 9% of the phosphates in the central
    # third;
    # the factor of two against the runner-up rejects a nearly isotropic
    # distribution,
    # where all three axes score alike and none of them is a membrane normal.
    if best[1] < 2.0 * second[1] or best[1] < 10.0:
        return None, (f"axis undetermined: central emptiness {mean_scores} "
                      f"(best {best[0]}={best[1]:.1f} does not separate from "
                      f"{second[0]}={second[1]:.1f})")
    basis = (f"emptiness of the central third between leaflets from phosphates on "
             f"{n_frames_sample} frames: { {k: round(v, 1) for k, v in mean_scores.items()} }; "
             f"best-over-second margin {best[1] / max(second[1], 1e-9):.0f}×")
    return best[0], basis


def fragment_sequences(u, channel_sel: str):
    """One sequence per connected fragment of the channel; chain ids are unused.

    A fragment is a connected component of the bond graph, so it is a chain whether or not
    the topology carries chain identifiers — a GROMACS .top read through the .itp parser
    does not. → [(fragment, residues sorted by resid, one-letter sequence)]; a residue
    outside ONE_LETTER becomes "X", which no motif pattern matches.
    """
    ch = u.select_atoms(channel_sel)
    frags = []
    for frag in ch.fragments:
        # The channel selection is made by subtraction, so anything unclassified
        # lands in
        # it. Fragments this small are not chains and carry no sequence to search.
        if len(frag) < 30:
            continue
        residues = sorted(frag.residues, key=lambda r: r.resid)
        seq = "".join(ONE_LETTER.get(r.resname, "X") for r in residues)
        frags.append((frag, residues, seq))
    return frags


def detect_filter(u, channel_sel: str, ax: int,
                  motif_pattern: str = DEFAULT_FILTER_MOTIF) -> dict:
    """Sequence motif plus geometric confirmation that its copies form a ring.

    The sequence alone is not evidence: the pattern is five residues long and can occur
    anywhere in a large protein. In a real selectivity filter the copies of the motif sit
    around the pore axis at the same height, one per subunit, so the geometry is what
    decides. Everything is read on frame 0, and the basis text says so.

    ax is the index of the membrane normal among the Cartesian axes. → a dict with
    filter_present, the matched motif as one-letter codes, the number of copies, the motif
    resids pooled and per copy, and the basis text recorded in the config.
    """
    frags = fragment_sequences(u, channel_sel)
    u.trajectory[0]
    matches = []
    for frag, residues, seq in frags:
        for m in re.compile(motif_pattern).finditer(seq):
            motif_res = residues[m.start():m.end()]
            ca = sum((r.atoms.select_atoms("name CA") for r in motif_res), start=u.atoms[[]])
            # Every residue of the match must have a C-alpha, or the stretch is not
            # backbone and its centroid would not mean anything. The centroid of
            # those
            # C-alphas is where this copy of the motif sits.
            if len(ca) == len(motif_res):
                matches.append((motif_res, ca.positions.mean(axis=0)))
    result = {"filter_present": False, "motif": None, "n_copies": 0,
              "motif_resids": [], "motif_resids_by_copy": [], "basis": ""}
    if not matches:
        result["basis"] = (f"motif {motif_pattern} not found in "
                           f"{len(frags)} connected fragments")
        return result
    ch = u.select_atoms(channel_sel)
    com = ch.center_of_mass()
    plane = [i for i in range(3) if i != ax]
    centroids = np.array([c for _, c in matches])
    # In-plane distance of each copy from the channel's centre of mass, and the
    # scatter of
    # the copies along the normal, both in A.
    r_off = np.hypot(*(centroids[:, plane].T - com[plane][:, None]))
    z_spread = centroids[:, ax].std()
    # A ring: at least two copies, since a single match arranges nothing; each
    # copy within
    # 10 A of the axis, so it lines the lumen and not the protein's outside; and
    # an axial
    # scatter under 6 A, so the copies line the same short stretch of the pore
    # rather than
    # sitting at different depths.
    ring_ok = len(matches) >= 2 and np.all(r_off < 10.0) and z_spread < 6.0
    if not ring_ok:
        result["basis"] = (f"motif found {len(matches)} times, but geometry does not "
                           f"confirm a ring: r_off={np.round(r_off,1)}, z_spread={z_spread:.1f} A")
        return result
    result.update({
        "filter_present": True,
        "motif": matches[0][0][0].resname and "".join(
            ONE_LETTER.get(r.resname, "X") for r in matches[0][0]),
        "n_copies": len(matches),
        "motif_resids": sorted({r.resid for mr, _ in matches for r in mr}),
        "motif_resids_by_copy": [[r.resid for r in mr] for mr, _ in matches],
        "basis": (f"motif in {len(matches)} copies; ring: centroid offsets from the axis "
                  f"{np.round(r_off,1).tolist()} A (<10), axial spread {z_spread:.1f} A (<6); "
                  f"frame 0"),
    })
    return result


def filter_oxygen_layers(u, motif_resids_by_copy: list[list[int]], channel_sel: str,
                         ax: int) -> list:
    """Oxygen rings of the motif: threonine OG1 + carbonyl O of every position.

    A ring is the atoms of ONE motif position across ALL copies; homo-oligomer
    resids run through the chains without repeating, so grouping is by position
    within the motif. Sole estimator of site layers, called by both autodetection
    and the feature block. Returns rings ordered by increasing mean axial
    coordinate.

    The oxygens that coordinate a permeating cation in a K-type filter are the backbone
    carbonyls of the motif plus the hydroxyl of its threonine (Zhou Y., Morais-Cabral
    J.H., Kaufman A., MacKinnon R. (2001) Nature 414:43–48, doi:10.1038/35102009), which
    is why OG1 is taken from motif position 0 and a carbonyl O from every position. Sites
    are then the gaps between neighbouring rings, and the count follows from the motif
    length rather than being asserted.

    → a list of AtomGroups, one per ring. The ordering is read from the frame currently
    loaded in u, so the caller decides which frame defines it.
    """
    ch = u.select_atoms(channel_sel)
    n_pos = len(motif_resids_by_copy[0])
    layers = []
    for name, positions in [("OG1", [0]), ("O", list(range(n_pos)))]:
        for j in positions:
            groups = []
            for copy in motif_resids_by_copy:
                # One atom per subunit at motif position j. A copy contributing
                # nothing is
                # skipped rather than refused: the ring at that position is
                # thinner, but
                # its mean axial coordinate is still defined.
                sel = ch.select_atoms(f"resid {copy[j]} and name {name}")
                if len(sel):
                    groups.append(sel)
            if groups:
                merged = groups[0]
                for g in groups[1:]:
                    merged = merged + g
                layers.append(merged)
    # Ascending mean axial coordinate, so consecutive entries are neighbouring
    # rings and
    # the gaps between them are the sites. The OG1 ring is built first but need
    # not stay
    # first: its place in the filter is decided by the structure, not by this loop.
    layers.sort(key=lambda ag: float(ag.positions[:, ax].mean()))
    return layers


def detect_n_subunits(u, channel_sel: str, ax: int, cylinder_A: float) -> tuple[int, str]:
    """Number of chains that reach the pore → (count, basis text).

    Counted as connected fragments of the bond graph, so it does not rely on chain
    identifiers, and only fragments that come near the axis are counted: a peripheral
    chain of the same assembly does not line the permeation pathway. cylinder_A is the
    confinement cylinder from the config; the test uses 1.5 times that radius, so a
    subunit counts when it reaches the neighbourhood of the pore and not only when it
    lines the lumen itself. Read on frame 0.
    """
    u.trajectory[0]
    ch = u.select_atoms(channel_sel)
    com = ch.center_of_mass()
    plane = [i for i in range(3) if i != ax]
    n = 0
    for frag in ch.fragments:
        if len(frag) < 30:
            continue
        pos = frag.positions
        r = np.hypot(pos[:, plane[0]] - com[plane[0]], pos[:, plane[1]] - com[plane[1]])
        # A single atom near the axis is enough: the question is whether the
        # chain reaches
        # the pore at all, not how much of it lies there.
        if np.any(r < cylinder_A * 1.5):
            n += 1
    return n, (f"connected fragments of ≥30 atoms touching the pore "
               f"(r<{cylinder_A * 1.5:.0f} A from the axis): {n}; frame 0")


def strided_crossing_scan(cfg: Config, rep: dict, candidates: dict[str, str]) -> dict:
    """Per-particle crossings for candidates (permeant) + flux sign + density.

    The sampling interval is set in TIME units (data.autodetect_scan.sample_ps):
    with multi-nanosecond gaps a bulk ion's rms displacement approaches half the
    cell height, and reservoir-to-reservoir diffusion through the periodic seam
    then fabricates full ladder passes that finer sampling would veto. The scan
    covers a contiguous window from the start (window_ps), which is enough for
    the rates that species ranking and flux signs are read from.

    candidates maps a name to a selection string; each is passed through the same ladder
    as the events step, so the ranking compares species measured identically. → a dict
    with the stride and the number of frames scanned, one sub-dict per candidate
    (particles, crossings, crossings per particle, net flux sign) and permeant_z_rel, the
    axial positions of the cations inside the pore that the site detection works on.
    """
    u = pio.open_replica(rep)
    ch = u.select_atoms(cfg["data.selections.channel"])
    pio.attach_prep(u, ch)
    phos = u.select_atoms(cfg["data.selections.phosphate"])
    ax = pio.AXIS_INDEX[cfg["data.pore.axis"]]
    plane = [i for i in range(3) if i != ax]
    cyl = cfg["data.pore.crossing_cylinder_radius_A"]
    n_frames = u.trajectory.n_frames
    _t0, dt = pio.frame_times_ps(u)
    # Both scan parameters are declared in ps and converted here with this
    # file's own frame
    # step, so the scan covers the same physical interval in trajectories written at
    # different rates. A requested interval finer than the frame step falls back
    # to every
    # frame rather than to a stride of zero.
    stride = max(1, int(round(cfg["data.autodetect_scan.sample_ps"] / dt)))
    span_frames = min(n_frames,
                      int(round(cfg["data.autodetect_scan.window_ps"] / dt)) + 1)
    frames = list(range(0, span_frames, stride))
    groups = {name: u.select_atoms(sel) for name, sel in candidates.items()}
    T = len(frames)
    data = {name: (np.empty((T, len(ag))), np.empty((T, len(ag))))
            for name, ag in groups.items()}
    times = np.empty(T)
    p_low = np.empty(T)
    p_high = np.empty(T)
    box_z = np.empty(T)
    ch_z = []
    for k, fi in enumerate(frames):
        ts = u.trajectory[fi]
        times[k] = ts.time
        box_z[k] = ts.dimensions[ax]
        # The reference here is the channel centre of mass, not the anchor
        # selection: the
        # ladder planes and the particle coordinates are all measured against
        # it, so the
        # zones are the same as long as the reference is one and the same.
        com = ch.center_of_mass()
        # Anchor-relative minimum-image axial coordinates with a median leaflet
        # split: the same shift-invariant formulation the events collector uses.
        pz = relative_axial(phos.positions[:, ax], com[ax], box_z[k])
        med = np.median(pz)
        p_low[k] = pz[pz < med].mean()
        p_high[k] = pz[pz >= med].mean()
        # Axial extent of the channel as the 2nd and 98th percentiles of its
        # atoms, in A:
        # a diagnostic robust to a single stretched-out loop, unlike min and max.
        chpos = relative_axial(ch.positions[:, ax], com[ax], box_z[k])
        ch_z.append(np.percentile(chpos, [2, 98]))
        for name, ag in groups.items():
            pos = ag.positions
            data[name][0][k] = relative_axial(pos[:, ax], com[ax], box_z[k])
            r = np.hypot(pos[:, plane[0]] - com[plane[0]], pos[:, plane[1]] - com[plane[1]])
            data[name][1][k] = r
    ch_z = np.array(ch_z).mean(axis=0)
    p_low, chan_lo, chan_hi, p_high = zone_planes(p_low, p_high)
    mid = 0.5 * (p_low + p_high)
    out = {"stride": stride, "n_frames_scanned": T, "channel_z_extent": ch_z.tolist()}
    for name, (z, r) in data.items():
        ev, n_cross, _ = ladder_pass(times, z, r, p_low, p_high, chan_lo, chan_hi, cyl,
                                     box_z=box_z)
        n_particles = z.shape[1]
        out[name] = {
            "n_particles": int(n_particles),
            "crossings": int(n_cross.sum()),
            # Per particle, because water outnumbers the ions by orders of
            # magnitude and
            # would win any comparison of raw crossing counts.
            "per_particle": float(n_cross.sum() / max(n_particles, 1)),
            # Sum of the signed directions: +1 for each upward pass, −1 for each
            # downward
            # one, so the sign is the net flux and it is zero when the passes cancel.
            "net_sign": int(np.sign(sum(e["direction"] for e in ev))) if ev else 0,
        }
    # Axial permeant density inside the cylinder, for the site peaks. Positions
    # are taken
    # relative to the membrane midplane, so the peaks of several replicas can be
    # pooled
    # even though each has its own channel centre of mass.
    z, r = data["ion"]
    mask = (r <= cyl) & (z >= chan_lo[:, None]) & (z <= chan_hi[:, None])
    perm_rel = (z - mid[:, None])[mask]
    out["permeant_z_rel"] = perm_rel.tolist()
    return out


def detect_coordination_cutoff(cfg: Config, rep: dict, n_frames_sample: int) -> tuple[float | None, str]:
    """First minimum of the ion-water g(r) within a physically sensible window.

    The cutoff separates the first hydration shell from the second, and it is what the
    coordination number of the permeating ion is counted with. Measuring it per system
    rather than declaring one value keeps the count comparable between a K+ and a Na+
    trajectory, whose shells sit at different radii. → (cutoff in A, basis text), or
    (None, why) when the histogram is too sparse to have a minimum.
    """
    from MDAnalysis.lib.distances import capped_distance
    u = pio.open_replica(rep)
    ions = u.select_atoms(cfg["data.selections.ion"])
    wox = u.select_atoms(cfg["data.selections.water_oxygen"])
    idx = np.linspace(0, u.trajectory.n_frames - 1, n_frames_sample, dtype=int)
    # Ion-oxygen distances from 2 to 6 A in 0.05 A bins, so the cutoff is
    # located to within
    # one bin: fine enough for a coordination count, coarse enough to hold
    # counts per bin.
    edges = np.arange(2.0, 6.0, 0.05)
    hist = np.zeros(len(edges) - 1)
    for fi in idx:
        ts = u.trajectory[int(fi)]
        # Distances through the box: the pairs are collected over the whole
        # system, most
        # of them in bulk water far from the channel, where the cell faces matter.
        pairs, dists = capped_distance(ions.positions, wox.positions, max_cutoff=6.0,
                                       box=ts.dimensions, return_distances=True)
        h, _ = np.histogram(dists, bins=edges)
        hist += h
    centers = 0.5 * (edges[1:] + edges[:-1])
    # Divide out the 4*pi*r^2 growth of the shell volume. The remaining
    # constants (bin
    # width, bulk density, number of frames) scale the whole curve and do not
    # move its
    # extrema, so this is g(r) up to a factor, which is all that is needed here.
    shell = hist / (4 * np.pi * centers ** 2)
    # The search window brackets the first hydration shell of a monovalent
    # cation: below
    # 2.2 A lies the repulsive core where no oxygen sits, and beyond 4.8 A the second
    # shell begins and would offer its own minimum.
    window = (centers > 2.2) & (centers < 4.8)
    # Fewer than a hundred pairs in the window: the minimum would be a gap
    # between counts
    # rather than a feature of the distribution.
    if hist[window].sum() < 100:
        return None, "too few ion-water pairs for g(r)"
    # window is boolean, so the product zeroes everything outside it and the
    # peak is the
    # first-shell maximum. The minimum is then searched strictly beyond that
    # peak; the
    # excluded bins are set to +inf so that argmin cannot return one of them.
    peak_i = np.argmax(shell * window)
    after = shell.copy()
    after[centers <= centers[peak_i]] = np.inf
    after[centers > 4.8] = np.inf
    min_i = int(np.argmin(after))
    cutoff = float(centers[min_i])
    basis = (f"ion-water g(r) on {n_frames_sample} frames: peak {centers[peak_i]:.2f} A, "
             f"first minimum {cutoff:.2f} A (window 2.2–4.8 A)")
    return cutoff, basis


def detect_sites_from_density(z_rel: np.ndarray, bin_A: float = 0.5
                              ) -> tuple[int | None, list[float], str]:
    """Axial density peaks with a self-check across trajectory halves.

    The fallback for a channel with no sequence-defined filter: where the permeant lingers
    is read off its own axial density. z_rel are the pooled axial positions in A of the
    permeant inside the pore, from strided_crossing_scan. → (number of peaks, their
    positions in A, basis text) or (None, [], why). The basis says in as many words that
    this discretizes a density and does not count physical binding sites.
    """
    from scipy.signal import find_peaks
    if len(z_rel) < 200:
        return None, [], f"too few permeant observations in the pore ({len(z_rel)})"
    z_rel = np.asarray(z_rel)
    # The two halves are consecutive in time. A peak count that differs between
    # them is
    # either drift or noise, and in both cases it is not a property of the system.
    halves = np.array_split(z_rel, 2)
    # 1st to 99th percentile: the rare excursion to the mouth would otherwise
    # stretch the
    # histogram range and leave the pore itself in a handful of bins.
    lo, hi = np.percentile(z_rel, [1, 99])
    edges = np.arange(lo, hi + bin_A, bin_A)
    centers = 0.5 * (edges[1:] + edges[:-1])

    def peaks_of(v):
        """Indices of the density peaks of one sample, on the shared bin edges."""
        h, _ = np.histogram(v, bins=edges)
        if h.max() == 0:
            return []
        # Three-bin moving average over 0.5 A bins: it removes single-bin
        # fluctuations
        # while leaving structure on the 1.5 A scale intact. Prominence must
        # reach a fifth
        # of the highest smoothed count, so the shoulder of a strong peak is not
        # counted as
        # a position of its own.
        sm = np.convolve(h, np.ones(3) / 3, mode="same")
        pk, _ = find_peaks(sm, prominence=0.2 * sm.max())
        return list(pk)

    p1, p2 = peaks_of(halves[0]), peaks_of(halves[1])
    if len(p1) == 0 or len(p1) != len(p2):
        return None, [], (f"density peaks do not reproduce across halves: "
                          f"{len(p1)} vs {len(p2)}")
    # Positions come from the pooled sample, which locates them best. Should
    # pooling yield
    # a different number of peaks than the halves agreed on, the reproducible
    # count wins
    # and the first half supplies the positions.
    pos = [round(float(centers[k]), 2) for k in peaks_of(z_rel)]
    if len(pos) != len(p1):
        pos = [round(float(centers[k]), 2) for k in p1]
    return len(p1), pos, (f"{len(p1)} axial permeant density peaks at {pos} A from the anchor; "
                          f"consistent across halves ({len(z_rel)} observations). This "
                          f"discretizes the density; it does not count physical sites")


def run_step(cfg: Config) -> None:
    """Measure every open parameter of the run and write the measurements into the config.

    The numbered stages below run in a fixed order because each uses what the previous
    ones established: selections come first because everything is expressed in them, then
    the axis, then the filter and the anchor that all axial coordinates are counted from,
    and only then the quantities that need a geometry. A stage whose result later stages
    depend on ends the step with a refusal instead of leaving them to measure something
    against an unknown reference.

    Every value is written with set_key together with the basis text that produced it, so
    the config records not just what was used but why. Nothing here has a default to fall
    back on: an undetermined parameter is named as a question for a human.
    """
    report: dict = {"measured": {}, "refusals": []}
    with step_output(cfg, "autodetect") as out:
        log = StepLog(out)
        refusals: list[Refusal] = []

        # Measurements that do not need every replica are made on the first one;
        # those that
        # must hold for the whole run (topologies, frame steps) loop over all of
        # them.
        first_cond = cfg["data.conditions"][0]
        first_rep = first_cond["replicas"][0]
        u = pio.open_replica(first_rep)

        # 1. Atom groups by composition: classified on every unique topology and
        # merged as a union of names per role, since residue NAMES vary between
        # force-field setups (SOL vs TIP3 water) while the composition does not.
        # A name that changes role between topologies is a genuine conflict.
        sels, bases, refs = detect_selections(u)
        refusals += refs
        if refs:
            _finish(cfg, out, log, report, refusals)
            return
        per_top_roles = [bases["roles"]]
        seen_tops = {first_rep["topology"]}
        for cond_id, rep in cfg.replicas():
            # One classification per distinct topology file. Replicas that share
            # a topology
            # share their residue names, so re-reading them would add nothing.
            if rep["topology"] in seen_tops:
                continue
            seen_tops.add(rep["topology"])
            u_check = pio.open_replica(rep)
            _s2, b2, r2 = detect_selections(u_check)
            refusals += r2
            if r2:
                _finish(cfg, out, log, report, refusals)
                return
            per_top_roles.append(b2["roles"])
        union, err = merge_roles(per_top_roles)
        if err:
            refusals.append(Refusal("data.selections", err))
            _finish(cfg, out, log, report, refusals)
            return
        sels = selections_from_roles(union)
        union_note = (f"; names unioned over {len(per_top_roles)} topologies"
                      if len(per_top_roles) > 1 else "")
        for name, sel in sels.items():
            cfg.set_key(f"data.selections.{name}", sel, "detected",
                        bases["selections"] + union_note)
        # Frame preparation: done once, before all subsequent measurements.
        pio.attach_prep(u, u.select_atoms(sels["channel"]))
        log.say(f"selections: {json.dumps(sels, ensure_ascii=False)}")

        # 2. Frame step: from the first two frame times of EVERY file; box cell.
        dts = {}
        for cond_id, rep in cfg.replicas():
            uu = pio.open_replica(rep)
            t0, dt = pio.frame_times_ps(uu)
            dts[f"{cond_id}/{rep['id']}"] = dt
            del uu
        # Windows and horizons are declared in ps and converted per file, so
        # different
        # frame steps are workable — but they change what "one frame" means in a
        # pooled
        # table, so they are said out loud.
        if len(set(dts.values())) > 1:
            log.say(f"WARNING: frame step differs between files: {dts}")
        report["measured"]["dt_ps_per_file"] = dts
        report["measured"]["triclinic"] = pio.triclinic_flag(u)
        log.say(f"frame step, ps: {dts}; triclinic cell: {report['measured']['triclinic']}")

        # 3. Pore axis: from leaflet separation, with a no-information floor.
        phos = u.select_atoms(sels["phosphate"])
        axis, axis_basis = detect_axis(u, phos, n_frames_sample=7)
        if axis is None:
            refusals.append(Refusal("data.pore.axis", axis_basis))
            _finish(cfg, out, log, report, refusals)
            return
        cfg.set_key("data.pore.axis", axis, "detected", axis_basis)
        ax = pio.AXIS_INDEX[axis]
        log.say(f"pore axis: {axis} ({axis_basis})")

        # 4. Filter: motif + geometry; subunits; anchor.
        filt = detect_filter(u, sels["channel"], ax,
                     cfg["system.arch_profile.filter_motif_pattern"])
        cfg.set_key("system.arch_profile.filter_present", filt["filter_present"],
                    "detected", filt["basis"])
        log.say(f"filter: {filt['filter_present']} ({filt['basis']})")
        # The confinement cylinder is part of the problem definition and is normally
        # declared; 10 A stands in only for counting subunits, which needs a
        # length scale
        # for "near the axis" and not the definition itself.
        cyl_declared = cfg["data.pore.crossing_cylinder_radius_A"]
        n_sub, sub_basis = detect_n_subunits(u, sels["channel"], ax, cyl_declared or 10.0)
        cfg.set_key("system.arch_profile.n_subunits", n_sub, "detected", sub_basis)
        # The anchor is the zero of every axial offset in the run. With a filter
        # it is the
        # filter itself, so a bin means the same part of the pore in every frame
        # and in
        # every system; without one there is no such landmark and the whole
        # backbone is
        # used, which at least moves with the protein rather than with the box.
        if filt["filter_present"]:
            cfg.set_key("system.arch_profile.filter_motif", filt["motif"], "detected", filt["basis"])
            anchor_sel = (f"({sels['channel']}) and name CA and resid "
                          + " ".join(str(r) for r in filt["motif_resids"]))
            anchor_basis = f"CA atoms of the filter motif residues (resid {filt['motif_resids']})"
            cfg.set_key("system.architecture", "k_filter_channel", "detected",
                        "K-type selectivity filter found")
        else:
            anchor_sel = f"({sels['channel']}) and name CA"
            anchor_basis = "CA of the whole channel: no filter, anchor = backbone center of mass"
            cfg.set_key("system.architecture", "single_file_channel", "detected",
                        "no filter; the channel is a single-file pore (water wire)")
        cfg.set_key("data.selections.anchor", anchor_sel, "detected", anchor_basis)
        log.say(f"anchor: {anchor_sel}")

        # 5. Permeant from per-particle crossings; direction per condition.
        # All three mobile species are scanned, cations, anions and water, so
        # which of them
        # the channel actually conducts is a measurement rather than an
        # assumption about
        # what a channel of this name is supposed to do.
        candidates = {"ion": sels["ion"], "ion_negative": sels["ion_negative"],
                      "water": sels["water_oxygen"]}
        scan_by_cond: dict[str, list] = {}
        for cond in cfg["data.conditions"]:
            scan_by_cond[cond["id"]] = []
            for rep in cond["replicas"]:
                sc = strided_crossing_scan(cfg, rep, candidates)
                scan_by_cond[cond["id"]].append(sc)
                log.say(f"scan {cond['id']}/{rep['id']} (stride {sc['stride']}): "
                        + ", ".join(f"{n}: {sc[n]['crossings']} crossings over "
                                    f"{sc[n]['n_particles']} particles"
                                    for n in candidates))
        per_particle = {n: float(np.mean([sc[n]["per_particle"]
                                          for scans in scan_by_cond.values() for sc in scans]))
                       for n in candidates}
        pp_txt = {k: round(v, 3) for k, v in per_particle.items()}
        declared_perm = cfg["system.permeant"]
        # A human declaration wins, but the measurement is still made and logged
        # next to it:
        # the declaration is then visible as a choice, not as the only available
        # number.
        if declared_perm is not None and cfg.origin_of("system.permeant")[0] == "declared":
            best_name = next((n for n, sel in candidates.items()
                              if f"resname {declared_perm}" in sel), "ion")
            log.say(f"permeant declared by a human: {declared_perm}; measured "
                    f"per-particle crossings: {pp_txt}")
        else:
            ranked = sorted(per_particle.items(), key=lambda kv: -kv[1])
            if ranked[0][1] <= 0:
                refusals.append(Refusal("system.permeant",
                                        f"no candidate crosses the membrane: {pp_txt}"))
                _finish(cfg, out, log, report, refusals)
                return
            best_name = ranked[0][0]
            if best_name == "water":
                # The per-particle crossing rule keeps water from winning on
                # sheer numbers. Water winning per particle too means the
                # channel conducts both water and ions; the measurement is
                # ambiguous, so refuse and print the full table.
                refusals.append(Refusal(
                    "system.permeant",
                    f"water beats the ions even in per-particle crossings: {pp_txt}; "
                    f"declare the permeant in the config (system.permeant)"))
                _finish(cfg, out, log, report, refusals)
                return
            perm_resnames = candidates[best_name].split("resname ")[1].split(" and")[0]
            cfg.set_key("system.permeant", perm_resnames, "detected",
                        f"per-particle crossings: {pp_txt} — leader {best_name} "
                        f"with a {ranked[0][1] / max(ranked[1][1], 1e-9):.1f}× margin over second")
            log.say(f"permeant: {perm_resnames} (per particle: {pp_txt})")
        permeant_sel = candidates[best_name]

        # The conduction direction of a condition is what sign-carrying columns are
        # canonicalized against, so it must be one direction for the whole condition.
        # Replicas that produced no crossing say nothing and are left out of the
        # vote;
        # replicas that disagree make the condition undecided rather than averaged.
        for cond in cfg["data.conditions"]:
            signs = [sc[best_name]["net_sign"] for sc in scan_by_cond[cond["id"]]]
            nonzero = sorted(set(s for s in signs if s != 0))
            if len(nonzero) == 1:
                cond["direction"] = nonzero[0]
                cfg.origins[f"data.conditions[{cond['id']}].direction"] = {
                    "origin": "detected",
                    "basis": f"sign of the net permeant flux across replicas: {signs}"}
                log.say(f"direction {cond['id']}: {nonzero[0]} (replica signs {signs})")
            else:
                refusals.append(Refusal(
                    f"data.conditions[{cond['id']}].direction",
                    f"replica flux signs disagree or are all zero: {signs}"))

        # 6. Charge scaling per replica: from the topology.
        for cond in cfg["data.conditions"]:
            for rep in cond["replicas"]:
                uu = pio.open_replica(rep)
                # Mean |q| of the permeant in e, read from the topology rather
                # than from the
                # run notes: under electronic continuum correction the ion
                # charge is scaled
                # (about 0.7 e for a monovalent ion), and that scaling belongs to the
                # protocol record of the replica it was applied in.
                q = float(np.abs(uu.select_atoms(permeant_sel).charges).mean())
                rep["protocol"]["charge_scaling"] = round(q, 4)
                cfg.origins[f"data.conditions[{cond['id']}].replicas[{rep['id']}]"
                            ".protocol.charge_scaling"] = {
                    "origin": "detected",
                    "basis": f"|q| of the monovalent permeant in the topology = {q:.4f} e"}
                del uu

        # 7. Coordination cutoff: median over conditions, spread is printed.
        # One value for the whole run, since a coordination number counted with
        # different
        # cutoffs in different conditions is not one column. The median resists
        # a condition
        # whose g(r) came out noisy, and the spread is printed so a large
        # disagreement
        # between conditions is visible rather than averaged away.
        cuts = []
        for cond in cfg["data.conditions"]:
            c, cb = detect_coordination_cutoff(cfg, cond["replicas"][0], n_frames_sample=5)
            if c is not None:
                cuts.append(c)
                log.say(f"coordination cutoff {cond['id']}: {c:.2f} A ({cb})")
        if cuts:
            cfg.set_key("features.coordination_cutoff_A", float(np.median(cuts)), "detected",
                        f"median over conditions from g(r): {cuts}; spread "
                        f"{max(cuts) - min(cuts):.2f} A")
        else:
            refusals.append(Refusal("features.coordination_cutoff_A", "g(r) could not be built"))

        # 8. Pore extent and bins: fixed offsets from the anchor.
        anchor_ag = u.select_atoms(anchor_sel)
        ch = u.select_atoms(sels["channel"])
        exts = []
        for fi in np.linspace(0, u.trajectory.n_frames - 1, 5, dtype=int):
            u.trajectory[int(fi)]
            az = anchor_ag.center_of_mass()[ax]
            # Extent measured from the anchor, so the numbers written to the
            # config are the
            # offsets the profile and the bins are later expressed in. The 5th
            # and 95th
            # percentiles rather than the extremes: a single protruding loop
            # would push the
            # range out by several angstroms without widening the pore.
            chz = ch.positions[:, ax] - az
            exts.append(np.percentile(chz, [5, 95]))
        ext = np.mean(exts, axis=0)
        # Rounded outwards to whole angstroms, so the bin edges are readable
        # numbers and
        # small changes in the sampled frames do not shift them.
        z_lo = float(np.floor(ext[0]))
        z_hi = float(np.ceil(ext[1]))
        # The bin count is a schema constant, the same for every system: it is
        # the bin WIDTH
        # that adapts to the channel, so tables of different proteins have the
        # same columns
        # and each bin means the same fraction of that channel's pore.
        n_bins = cfg["data.pore.n_bins"]
        bin_w = round((z_hi - z_lo) / n_bins, 2)
        ext_basis = (f"5–95 percentiles of channel axial coordinates from the anchor "
                     f"on 5 frames: [{ext[0]:.1f}, {ext[1]:.1f}] A")
        # The profile is measured 4 A beyond the channel at each end while the
        # bins stop at
        # the channel: the widening at the mouth is then inside the measured
        # range, where
        # the fence can reject it, instead of falling on the last bin of the table.
        cfg.set_key("data.pore.z_low_offset_A", z_lo - 4.0, "detected",
                    ext_basis + "; −4 A at the mouth")
        cfg.set_key("data.pore.z_high_offset_A", z_hi + 4.0, "detected",
                    ext_basis + "; +4 A at the mouth")
        cfg.set_key("data.pore.bin_low_offset_A", z_lo, "detected", ext_basis)
        cfg.set_key("data.pore.bin_width_A", bin_w, "detected",
                    f"({z_hi}-{z_lo})/{n_bins} bins")
        log.say(f"pore extent from the anchor: [{z_lo}, {z_hi}] A; bin {bin_w} A × {n_bins}")

        # 9. Pore cylinder for the columns: from the R(z) profile on sampled frames.
        # Radii are resolved for the same concatenation that becomes occl below,
        # so the
        # radii array and the position array share one atom order — the profile
        # reads them
        # by index, and a mismatch would assign each atom a neighbour's radius.
        radii, steps = resolve_radii(
            (ch + u.select_atoms(sels["lipid"])).resnames,
            (ch + u.select_atoms(sels["lipid"])).names,
            (ch + u.select_atoms(sels["lipid"])).elements)
        report["measured"]["vdw_ladder_steps"] = steps
        prof = PoreProfiler(
            z_offsets=np.arange(z_lo, z_hi + 0.01, cfg["data.pore.z_step_A"]),
            search_radius_A=cfg["data.pore.probe_search_radius_A"],
            dr_A=cfg["data.pore.probe_grid_dr_A"],
            n_theta=cfg["data.pore.probe_grid_n_theta"],
            fence_max_gap_deg=cfg["data.pore.fence_max_gap_deg"],
            lipschitz_tol_A=cfg["data.pore.lipschitz_tol_A"],
            slab_pad_A=cfg["data.pore.slab_pad_A"])
        # Channel and lipids together: both occupy volume in the pore, and acyl
        # chains
        # entering the lumen are a documented closure mechanism (see the module
        # docstring
        # of pore.py). Water and ions are left out — the quantity is the space
        # available
        # to the mobile phase.
        occl = ch + u.select_atoms(sels["lipid"])
        r_maxes = []
        for fi in np.linspace(0, u.trajectory.n_frames - 1, 3, dtype=int):
            u.trajectory[int(fi)]
            anchor = anchor_ag.center_of_mass()
            res = prof.profile(occl.positions, radii, anchor)
            # Widest radius the estimator actually measured on this frame;
            # slices it refused
            # are excluded, so a mouth slice cannot set the width of the cylinder.
            if np.any(res["search"]):
                r_maxes.append(float(np.nanmax(res["R"][res["search"]])))
        if not r_maxes:
            refusals.append(Refusal("data.pore.cylinder_radius_A",
                                    "the profile estimator refused on all sampled frames"))
        else:
            # This cylinder is the volume the water and ion columns are counted
            # in, so it has
            # to contain the whole lumen: the widest radius the profile
            # measured, plus a
            # 1.5 A margin, rather than a radius that clips the vestibule. The
            # median over
            # frames keeps one unusually open frame from setting it.
            cyl_r = round(float(np.median(r_maxes)) + 1.5, 1)
            cfg.set_key("data.pore.cylinder_radius_A", cyl_r, "detected",
                        f"median of the maximum search R(z) over 3 frames "
                        f"{np.round(r_maxes, 2).tolist()} + 1.5 A margin")
            log.say(f"pore cylinder: {cyl_r} A")

        # 10. Sites: motif rings, otherwise density peaks with a self-check.
        # With a filter the sites are structural — the gaps between neighbouring
        # oxygen
        # rings — and their number follows from the motif. Without one there is
        # nothing to
        # count rings on, and the peaks of the permeant's own density stand in;
        # the two
        # cases are told apart in the config by sites_basis, since they are not
        # the same
        # kind of quantity and a reader of the table must know which one it is.
        if filt["filter_present"]:
            layers = filter_oxygen_layers(u, filt["motif_resids_by_copy"], sels["channel"], ax)
            n_sites = len(layers) - 1
            cfg.set_key("system.arch_profile.n_sites", n_sites, "detected",
                        f"{len(layers)} motif oxygen rings → {n_sites} inter-ring sites")
            cfg.set_key("system.arch_profile.sites_basis", "motif_rings", "detected",
                        "selectivity filter rings")
            # Site occupancy is a Gaussian kernel of width sigma centred on a
            # site (see
            # features/named_sites.py). Taking sigma as half the ring spacing
            # puts the
            # neighbouring site two sigma away, so the kernels of adjacent sites
            # hardly
            # overlap and one ion is not counted in two sites at once.
            gaps = np.diff([float(g.positions[:, ax].mean()) for g in layers])
            cfg.set_key("features.site_sigma_A", round(float(np.median(np.abs(gaps))) / 2, 2),
                        "detected", f"half the median ring spacing {np.round(gaps, 2).tolist()}")
        else:
            # Pooled over every replica of every condition: one replica rarely
            # holds enough
            # in-pore observations for a peak count that reproduces across its
            # own halves.
            z_rel = np.concatenate([np.asarray(sc["permeant_z_rel"])
                                    for scans in scan_by_cond.values() for sc in scans])
            n_sites, site_pos, ns_basis = detect_sites_from_density(z_rel)
            if n_sites is None:
                refusals.append(Refusal("system.arch_profile.n_sites", ns_basis))
            else:
                cfg.set_key("system.arch_profile.n_sites", n_sites, "detected", ns_basis)
                cfg.set_key("system.arch_profile.sites_basis", "density_peaks", "detected",
                            ns_basis)
                cfg.set_key("features.site_centers_offset_A", site_pos, "detected", ns_basis)
                # Without rings there is no ring spacing to halve, so the width
                # comes from
                # the occupied range divided among the peaks: the same
                # convention of half a
                # spacing, with the spacing estimated from the extent of the density.
                lo, hi = np.percentile(z_rel, [1, 99])
                cfg.set_key("features.site_sigma_A",
                            round(float(hi - lo) / max(n_sites, 1) / 2, 2), "detected",
                            f"half the peak spacing over the range [{lo:.1f},{hi:.1f}] A")
            log.say(f"sites from density: {n_sites} ({ns_basis})")

        _finish(cfg, out, log, report, refusals)


def _finish(cfg: Config, out, log: StepLog, report: dict, refusals: list[Refusal]) -> None:
    """Write report.json, print the refusals and save the config, whether or not it refused.

    Called from every exit of the step, including the early returns. A partial run still
    leaves its report and still saves what it did measure: the values already written are
    correct measurements, and the refusals name exactly what is missing.
    """
    report["refusals"] = [r.as_dict() for r in refusals]
    # The origins of every key after the step: value, who set it (declared or
    # detected) and
    # the basis. This is the record that makes a run's parameters auditable
    # afterwards.
    report["origins_after"] = cfg.origins
    (out / "report.json").write_text(json.dumps(report, indent=1, ensure_ascii=False, default=str))
    if refusals:
        log.say("Autodetection REFUSALS (each names a config key):")
        for r in refusals:
            log.say(f"  {r.key}: {r.message}")
    else:
        log.say("autodetection finished without refusals")
    save(cfg)
    log.say(f"config updated in place: {cfg.source_path}")
    log.close()
