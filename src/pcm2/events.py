"""Crossings: read a provided annotation OR collect events from the trajectory.

Own collection is the primary path: a finite-state machine over zones along
the pore axis. Zones are separated by planes derived from lipid headgroups:
  0: below the lower headgroups; 1: inside the lower leaflet; 2: inside the
  channel (axial channel slab AND radially inside the confinement cylinder);
  3: inside the upper leaflet; 4: above the upper headgroups.
A crossing counts only on a complete one-way pass of the ladder; the
completion time is the frame time of the last transition (entry into the
extreme zone). A particle that traverses the axial channel slab OUTSIDE the
cylinder (through the bilayer past the pore) gets an invalid origin and
produces no event.

An ion that is not in an extreme zone on the first frame is left-censored:
its entry time was never observed and must not be replaced by the exit time.

With source=auto and both sources available the step REFUSES: it prints the
comparison and stops; a human must declare the choice in the config.
"""

from __future__ import annotations

import os

import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import io as pio
from .config import Config
from .runtime import StepLog, run_dir, step_output

INVALID = -9  # particle origin after traversing the channel slab outside the cylinder


def relative_axial(z: np.ndarray, anchor_z: float, box_z: float) -> np.ndarray:
    """Axial coordinate relative to the anchor, minimum-image in the periodic cell.

    The ladder must be invariant to a rigid shift of the whole frame: frame
    preparation re-centres on the channel's centre of mass, and when the
    unwrapping of a multi-chain channel picks a different periodic image the
    entire frame jumps by one cell. Absolute wrapped coordinates then teleport
    and reset particle origins mid-pass; anchor-relative minimum-image
    coordinates cancel any global shift identically.

    z, anchor_z and box_z are in A along one Cartesian axis (the membrane normal);
    box_z is that axis's box length in the frame the coordinates were read from. The
    result lies within half a box of zero and is positive on the anchor's high side.
    """
    d = np.asarray(z) - anchor_z
    return d - box_z * np.round(d / box_z)


def zone_planes(p_low: np.ndarray, p_high: np.ndarray):
    """Ladder planes derived from lipid headgroups only.

    The channel zone is the central half of the bilayer: midpoints between each
    leaflet plane and the membrane midplane. Single estimator of zone bounds:
    both autodetect (permeant selection) and the events step call it.

    Deriving them from the lipids alone keeps the ladder a property of the membrane: the
    zones are defined the same way for a channel of any length, and they do not move
    with the protein conformation whose readiness the labels are about.

    p_low, p_high: [T] mean axial coordinate of the lower and upper phosphates in A,
    anchor-relative. → (p_low, chan_lo, chan_hi, p_high), the four ladder planes in
    ascending order; the two inputs are passed through so one call yields all of them.
    """
    mid = 0.5 * (p_low + p_high)
    return p_low, 0.5 * (p_low + mid), 0.5 * (mid + p_high), p_high


def ladder_pass(times_ps: np.ndarray, z: np.ndarray, rxy: np.ndarray,
                p_low: np.ndarray, p_high: np.ndarray,
                chan_lo: np.ndarray, chan_hi: np.ndarray,
                cylinder_A: float,
                box_z: np.ndarray | None = None) -> tuple[list[dict], np.ndarray, np.ndarray]:
    """Zone state machine. z, rxy: [T, N] — axial coordinate and distance to the axis.

    box_z — axial box size per frame. The upper and lower solutions are
    connected through the periodic boundary: a particle jump larger than half
    the box is an image transfer through the boundary, not a pass through the
    membrane. Such a jump resets the particle origin WITHOUT counting an
    event; otherwise every ion wandering in bulk would yield spurious
    "complete passes" (zones 0 and 4 are adjacent across the box boundary).

    All lengths are in A and all times in ps. z and rxy are anchor-relative (see
    relative_axial); p_low, p_high, chan_lo, chan_hi are [T] and may move from frame to
    frame with the membrane. cylinder_A is data.pore.crossing_cylinder_radius_A, part of
    the problem definition rather than a geometry detail.

    Each event carries t_exit_ps — the frame time at which the particle entered the far
    extreme zone, i.e. the completion time the label is built on — t_entry_ps, the time
    of its last entry into the channel zone, entrance_observed, and direction: +1 for a
    pass from zone 0 to zone 4 (toward increasing anchor-relative axial coordinate), −1
    for the reverse.

    → (events, crossings_per_particle [N], left_censored_at_start [N]).
    """
    T, N = z.shape

    def zones_at(t: int) -> np.ndarray:
        """Zone of every particle in frame t; see the ladder in the module docstring."""
        zt = z[t]
        # Zone 1 (inside a leaflet) is the default; the masks below relabel the
        # extreme
        # zones, the channel slab and the upper leaflet, so what stays 1 is the band
        # between the lower headgroup plane and the lower edge of the channel slab.
        out = np.full(N, 1, dtype=np.int8)
        out[zt < p_low[t]] = 0
        out[zt > p_high[t]] = 4
        in_slab = (zt >= chan_lo[t]) & (zt <= chan_hi[t])
        out[in_slab & (rxy[t] <= cylinder_A)] = 2
        out[in_slab & (rxy[t] > cylinder_A)] = INVALID  # past the pore
        upper = (zt > chan_hi[t]) & (zt <= p_high[t])
        out[upper] = 3
        return out

    z0 = zones_at(0)
    # origin: the extreme zone the particle was last seen in, and the only thing
    # that can
    # authorise an event. A particle that starts anywhere else has no observed
    # origin, so
    # its first arrival in an extreme zone sets the origin instead of closing a pass.
    origin = np.where(z0 == 0, 0, np.where(z0 == 4, 4, INVALID)).astype(np.int8)
    left_censored_start = (z0 != 0) & (z0 != 4)
    # t_enter is the time of the last entry into the channel zone; enter_seen
    # says whether
    # that entry was observed as an entry rather than inferred. The entrance-anchored
    # labels are a control ablation (see labels.py) and use only the observed ones.
    t_enter = np.full(N, np.nan)
    enter_seen = np.zeros(N, dtype=bool)
    # A complete ladder pass, taken literally: between the extreme zones the
    # particle MUST be observed in the channel zone (inside the cylinder).
    # A direct 0↔4 jump between sampled frames is diffusion through the
    # periodic boundary, across which the extreme zones are adjacent; it
    # yields no event.
    seen_channel = np.zeros(N, dtype=bool)
    prev_zone = z0
    events: list[dict] = []
    n_cross = np.zeros(N, dtype=np.int64)

    for t in range(1, T):
        cur = zones_at(t)
        if box_z is not None:
            # A step longer than half the cell cannot be a physical displacement
            # between
            # two sampled frames: it is the same particle re-entering through
            # the opposite
            # face. The particle keeps its new position but loses its history.
            wrapped = np.abs(z[t] - z[t - 1]) > 0.5 * box_z[t]
            for i in np.flatnonzero(wrapped):
                c = int(cur[i])
                origin[i] = c if c in (0, 4) else INVALID
                enter_seen[i] = False
                seen_channel[i] = False
                t_enter[i] = np.nan
        else:
            wrapped = np.zeros(N, dtype=bool)
        # Zone changes across a wrap are excluded: the transition below would
        # read as a
        # move between zones when in fact the particle only changed periodic image.
        changed = np.flatnonzero((cur != prev_zone) & ~wrapped)
        for i in changed:
            c, p = int(cur[i]), int(prev_zone[i])
            if c == INVALID:
                # The particle crossed the channel slab outside the confinement
                # cylinder:
                # it went through the bilayer past the pore. Its origin is
                # discarded, so
                # whatever it does next cannot be closed as a pore crossing.
                origin[i] = INVALID
                seen_channel[i] = False
            elif c == 2:
                # entry into the channel zone: entry time observed only if the
                # particle came from a leaflet
                t_enter[i] = times_ps[t]
                enter_seen[i] = p in (1, 3)
                # The pass under construction counts as passing through the pore
                # only if
                # the particle reached the channel zone from a known extreme zone.
                seen_channel[i] = origin[i] in (0, 4)
            elif c == 0:
                # Arrival in an extreme zone closes a pass when it came from the
                # OTHER
                # extreme zone through the pore; either way it becomes the new
                # origin, and
                # the channel observation is cleared so the next pass starts
                # from scratch.
                if origin[i] == 4 and seen_channel[i]:
                    events.append({"ion_index": int(i), "t_exit_ps": float(times_ps[t]),
                                   "t_entry_ps": float(t_enter[i]),
                                   "entrance_observed": bool(enter_seen[i]),
                                   "direction": -1})
                    n_cross[i] += 1
                origin[i] = 0
                seen_channel[i] = False
            elif c == 4:
                if origin[i] == 0 and seen_channel[i]:
                    events.append({"ion_index": int(i), "t_exit_ps": float(times_ps[t]),
                                   "t_entry_ps": float(t_enter[i]),
                                   "entrance_observed": bool(enter_seen[i]),
                                   "direction": +1})
                    n_cross[i] += 1
                origin[i] = 4
                seen_channel[i] = False
        prev_zone = cur
    return events, n_cross, left_censored_start


def read_provided(path: Path, direction: int | None) -> pd.DataFrame:
    """Annotation in transit-times.dat format: ion, entry [ps], exit [ps].

    One whitespace-separated record per line; further columns are ignored, and a line whose
    first three fields do not parse as numbers is skipped, so a header or a comment line
    does not break the read. direction is taken from the condition, since the annotation
    records times only, and is written as 0 when the condition has none.
    The frame is marked source="provided" so a later comparison can tell the two apart.
    """
    rows = []
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) >= 3:
            try:
                ion, te, tx = int(parts[0]), float(parts[1]), float(parts[2])
            except ValueError:
                continue
            rows.append({"ion_index": ion, "t_entry_ps": te, "t_exit_ps": tx,
                         # an entry coinciding with the record start is truncated
                         "entrance_observed": te > 0.0,
                         "direction": direction if direction is not None else 0})
    df = pd.DataFrame(rows)
    df["source"] = "provided"
    return df


def collect_own_for_replica(cfg: Config, rep: dict, stride: int = 1,
                            max_frames: int | None = None) -> dict:
    """Full pass over the trajectory: only the arrays the state machine needs.

    → {"events", "n_cross_per_ion" [N], "left_censored_start" [N], "times_ps" [T],
    "n_ions"}. The step calls this with the defaults, so crossings are collected on every
    stored frame regardless of data.stride: a strided descriptor table then still refers
    to labels built from the full event record. stride and max_frames exist to shorten
    the pass over a long trajectory.
    """
    u = pio.open_replica(rep)
    groups = pio.resolve_groups(u, cfg)
    pio.attach_prep(u, groups["channel"])
    anchor_ag = u.select_atoms(cfg["data.selections.anchor"])
    if len(anchor_ag) == 0:
        raise pio.IoError("anchor selection is empty")
    ions = u.select_atoms(f"resname {cfg['system.permeant']}")
    # The permeant is a measured quantity, not a default: without it there is
    # nothing to
    # count crossings of, and a guess would silently decide which species the
    # run is about.
    if len(ions) == 0:
        raise pio.IoError("permeant not selected: run the autodetect step first")
    phos = groups["phosphate"]
    ax = pio.AXIS_INDEX[cfg["data.pore.axis"]]
    cyl = cfg["data.pore.crossing_cylinder_radius_A"]
    frames = range(0, u.trajectory.n_frames, stride)
    if max_frames is not None:
        frames = list(frames)[:max_frames]
    T = len(list(frames)) if not isinstance(frames, list) else len(frames)
    frames = list(frames)
    N = len(ions)
    times = np.empty(T)
    z = np.empty((T, N))
    rxy = np.empty((T, N))
    p_low = np.empty(T)
    p_high = np.empty(T)
    box_z = np.empty(T)
    plane_axes = [i for i in range(3) if i != ax]
    for k, fi in enumerate(frames):
        ts = u.trajectory[fi]
        times[k] = ts.time
        box_z[k] = ts.dimensions[ax]
        anchor = anchor_ag.center_of_mass()
        ipos = ions.positions
        # Anchor-relative minimum-image axial coordinates: invariant to rigid
        # frame shifts (see relative_axial).
        z[k] = relative_axial(ipos[:, ax], anchor[ax], box_z[k])
        # The pore axis is the line through the anchor centre of mass parallel to the
        # membrane normal, so rxy is the in-plane distance to that line, in A.
        # No minimum
        # image is needed in the plane: frame preparation centres the channel
        # and wraps
        # everything else into the cell, so each in-plane offset is already at
        # most half
        # a box and is the minimum image by construction.
        rxy[k] = np.hypot(ipos[:, plane_axes[0]] - anchor[plane_axes[0]],
                          ipos[:, plane_axes[1]] - anchor[plane_axes[1]])
        pz = relative_axial(phos.positions[:, ax], anchor[ax], box_z[k])
        # Leaflet split at the phosphate median: robust to where the anchor sits
        # along the pore (a filter anchor lies near one membrane face).
        med = np.median(pz)
        p_low[k] = pz[pz < med].mean() if np.any(pz < med) else np.nan
        p_high[k] = pz[pz >= med].mean() if np.any(pz >= med) else np.nan
    # A frame with all phosphates on one side of their own median has no second
    # leaflet
    # to place a plane on. The ladder would then be built on a plane derived
    # from nothing,
    # so the replica is refused instead of measured.
    if np.any(~np.isfinite(p_low)) or np.any(~np.isfinite(p_high)):
        raise pio.IoError("leaflet without phosphates on some frame: ladder zones undefined")
    p_low, chan_lo, chan_hi, p_high = zone_planes(p_low, p_high)
    events, n_cross, left0 = ladder_pass(times, z, rxy, p_low, p_high, chan_lo, chan_hi,
                                         cyl, box_z=box_z)
    return {"events": events, "n_cross_per_ion": n_cross, "left_censored_start": left0,
            "times_ps": times, "n_ions": N}


def compare_sources(own: pd.DataFrame, provided: pd.DataFrame,
                    taus_ps: list[float], dt_ps: float) -> dict:
    """Source diagnostics: event counts, exit shifts, matched fraction, Jaccard by τ.

    The two sources are never merged and neither is corrected against the other. This
    only measures how far apart they are, so that a human choosing events.source knows
    what the choice costs. dt_ps is the frame step of the replica: it is the finest
    disagreement the trajectory can resolve, hence the tolerance for "the same event".
    """
    out: dict = {"n_own": int(len(own)), "n_provided": int(len(provided))}
    if len(own) == 0 or len(provided) == 0:
        return out
    # Distance from every provided exit to the nearest own exit, in ps. Asymmetric on
    # purpose: it asks whether each annotated crossing has a counterpart here,
    # and the
    # count difference above already says whether one source has extra events.
    shifts = []
    for tx in provided["t_exit_ps"]:
        shifts.append(float(np.min(np.abs(own["t_exit_ps"].to_numpy() - tx))))
    shifts = np.asarray(shifts)
    out["exit_shift_ps"] = {"median": float(np.median(shifts)),
                            "p90": float(np.percentile(shifts, 90))}
    out["matched_within_frame_frac"] = float(np.mean(shifts <= dt_ps))
    out["jaccard_by_tau"] = {}
    t_max = max(own["t_exit_ps"].max(), provided["t_exit_ps"].max())
    grid = np.arange(0.0, t_max + dt_ps, dt_ps)
    # Agreement at the level that matters: each source is turned into the label
    # it would
    # produce at horizon τ, by exactly the rule of labels.py — a grid point is
    # positive
    # when a crossing completes in (t, t+τ]. The Jaccard index of the two label
    # vectors
    # says how far apart the two sources would push the trained monitor, which a
    # difference in raw event counts does not.
    for tau in taus_ps:
        def label(df):
            y = np.zeros(len(grid), dtype=bool)
            for tx in df["t_exit_ps"]:
                y |= (grid < tx) & (grid >= tx - tau)
            return y
        a, b = label(own), label(provided)
        # With no positive grid point in either source the index is undefined;
        # NaN, not
        # zero, since zero would read as complete disagreement.
        denom = np.sum(a | b)
        out["jaccard_by_tau"][str(tau)] = float(np.sum(a & b) / denom) if denom else np.nan
    return out


def run_step(cfg: Config) -> None:
    """Collect crossings for every replica and write events.parquet plus summary.json.

    Replicas are independent, so they are collected in separate processes; each worker
    reopens the config and its own trajectory. The pool size is capped by the number of
    replicas and by PCM2_MAX_WORKERS for machines with less memory than cores.
    """
    from concurrent.futures import ProcessPoolExecutor

    source = cfg["events.source"]
    with step_output(cfg, "events") as out:
        log = StepLog(out)
        all_events: list[pd.DataFrame] = []
        summary: dict = {"replicas": {}, "source_declared": source}
        jobs = list(cfg.replicas())
        with ProcessPoolExecutor(max_workers=min(len(jobs), int(os.environ.get("PCM2_MAX_WORKERS", "14")))) as ex:
            futs = {ex.submit(_collect_worker, cfg.source_path, cond_id, rep["id"]): (cond_id, rep)
                    for cond_id, rep in jobs}
            results = {}
            for fut, (cond_id, rep) in futs.items():
                results[(cond_id, rep["id"])] = fut.result()
        for cond in cfg["data.conditions"]:
            for rep in cond["replicas"]:
                key = (cond["id"], rep["id"])
                res = results[key]
                own_df = pd.DataFrame(res["events"])
                # A non-conducting replica is a legitimate result, and its frame must
                # still carry the full column set: the concatenation below is
                # what feeds
                # events.parquet, and a column missing here would be missing
                # everywhere.
                if len(own_df) == 0:
                    own_df = pd.DataFrame(columns=["ion_index", "t_exit_ps", "t_entry_ps",
                                                   "entrance_observed", "direction"])
                own_df["source"] = "own"
                provided_df = None
                if rep["events"]:
                    provided_df = read_provided(Path(rep["events"]), cond.get("direction"))
                # Frame step of this replica, from the first two collected
                # times; it is
                # the tolerance of the source comparison, not a config value.
                dt = float(np.diff(res["times_ps"][:2])[0])
                rep_sum = {
                    "n_own": int(len(own_df)),
                    "n_provided": None if provided_df is None else int(len(provided_df)),
                    "n_ions": res["n_ions"],
                    "left_censored_at_start": int(np.sum(res["left_censored_start"])),
                    "net_flux_sign": int(np.sign(own_df["direction"].sum())) if len(own_df) else 0,
                }
                if provided_df is not None:
                    rep_sum["comparison"] = compare_sources(
                        own_df, provided_df, cfg["labels.tau_ps"], dt)
                summary["replicas"][f"{key[0]}/{key[1]}"] = rep_sum
                log.say(f"[{key[0]}/{key[1]}] own={rep_sum['n_own']} "
                        f"provided={rep_sum['n_provided']} flux_sign={rep_sum['net_flux_sign']} "
                        f"left-censored at start: {rep_sum['left_censored_at_start']}")
                if source == "own":
                    chosen = own_df
                elif source == "provided":
                    if provided_df is None:
                        raise RuntimeError(f"{key}: source=provided, but no annotation file given")
                    chosen = provided_df
                else:
                    chosen = None
                if chosen is not None:
                    chosen = chosen.copy()
                    chosen["condition"] = key[0]
                    chosen["replica"] = key[1]
                    all_events.append(chosen)

        if source == "auto":
            # Both sources are available at least somewhere: refuse, printing
            # the comparison.
            (out / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
            (out / "REFUSED.txt").write_text(
                "events.source=auto: the system refuses to choose the label source.\n"
                "The source comparison is in summary.json. Declare events.source in the config "
                "(key events.source: own|provided) and rerun the step.\n")
            log.say("REFUSED: events.source=auto — a human must declare the source choice.")
            log.close()
            return

        ev = pd.concat(all_events, ignore_index=True) if all_events else pd.DataFrame()
        # Conduction direction: all replicas of a condition must agree on the
        # flux sign.
        for cond in cfg["data.conditions"]:
            signs = [summary["replicas"][f"{cond['id']}/{r['id']}"]["net_flux_sign"]
                     for r in cond["replicas"]]
            summary.setdefault("condition_flux_signs", {})[cond["id"]] = signs
            declared = cond.get("direction")
            # A replica with no crossings has sign 0 and says nothing about the
            # direction,
            # so it neither creates a disagreement nor contradicts the config.
            nonzero = [s for s in signs if s != 0]
            if nonzero and len(set(nonzero)) > 1:
                log.say(f"WARNING: condition {cond['id']} replicas disagree on flux sign: {signs}")
            if declared is not None and nonzero and any(s != declared for s in nonzero):
                raise RuntimeError(
                    f"condition {cond['id']}: measured flux sign {signs} contradicts "
                    f"direction {declared} in the config — a config silently contradicting "
                    "the data must fail the run")
        ev.to_parquet(out / "events.parquet")
        (out / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
        # Characteristic process times, reported next to the declared τ
        # horizons: a horizon
        # far shorter than the transit time or far longer than the interval between
        # crossings describes a different problem, and the two numbers make that
        # visible.
        if len(ev):
            transit = (ev["t_exit_ps"] - ev["t_entry_ps"]).dropna()
            inter = ev.sort_values("t_exit_ps").groupby(["condition", "replica"])["t_exit_ps"].diff().dropna()
            log.say(f"median transit time: {float(transit.median()):.0f} ps; "
                    f"median interval between crossings: {float(inter.median()):.0f} ps; "
                    f"declared horizons τ: {cfg['labels.tau_ps']} ps")
        log.close()


def _collect_worker(cfg_path: Path, cond_id: str, rep_id: str) -> dict:
    """One replica in a separate process: the config path travels, the Universe does not."""
    from . import config as cmod
    cfg = cmod.load(cfg_path)
    rep = next(r for c, r in cfg.replicas() if c == cond_id and r["id"] == rep_id)
    return collect_own_for_replica(cfg, rep)


def load_events(cfg: Config) -> pd.DataFrame:
    """The events artifact of this run; the only way later steps see crossings.

    A refusal leaves no events.parquet, so a step that depends on labels stops here
    rather than proceeding with an empty event record.
    """
    p = run_dir(cfg) / "events" / "events.parquet"
    if not p.exists():
        raise FileNotFoundError(
            f"{p} missing: the events step has not run or refused (see events/REFUSED.txt)")
    return pd.read_parquet(p)
