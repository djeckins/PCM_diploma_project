"""Label "a crossing completes within τ" plus censoring.

y(t) = 1 ⟺ the COMPLETION time of at least one crossing lies in (t, t + τ].
Built per trajectory; the window never crosses a trajectory boundary.

Right censoring: a frame with an incomplete window is unknown, not negative;
valid = window complete OR y = 1. A confirmed positive is never censored.
The number of dropped frames is printed: frames silently labelled zero would
deflate the base rate.

The "entrance" anchor is a circularity-control ablation: an event without an
observed entry time is EXCLUDED from it.
"""

from __future__ import annotations

import json

import numpy as np
import pandas as pd

from . import io as pio
from .config import Config
from .events import load_events
from .runtime import StepLog, run_dir, step_output


def label_one_trajectory(times_ps: np.ndarray, exit_times: np.ndarray,
                         entry_times: np.ndarray, tau_ps: float,
                         t_last_ps: float) -> pd.DataFrame:
    """Label and validity columns of ONE trajectory at ONE horizon, one row per frame.

    times_ps — frame times of this trajectory at the analysis stride, ascending;
    exit_times and entry_times — completion and entrance times of its crossings;
    t_last_ps — time of its last recorded frame. Everything in ps, and everything
    from the same trajectory: a horizon that reached into a neighbouring replica
    would join two unrelated pieces of dynamics.

    Returned in the frame order given: y and valid for the exit anchor, then the
    same pair for the entrance anchor.
    """
    y = np.zeros(len(times_ps), dtype=np.int8)
    for tx in exit_times:
        # Frames t for which tx lies in (t, t+τ]: the crossing is still ahead
        # (t < tx) and completes no later than the end of the horizon
        # (tx <= t + τ). The frame on which the crossing completes is left
        # negative — there the passage is over and there is nothing to forewarn.
        # OR-ing the 0/1 masks takes the union over crossings, so the overlapping
        # horizons of two ions leaving in quick succession count once.
        y |= ((times_ps < tx) & (times_ps >= tx - tau_ps)).astype(np.int8)
    # Right censoring. A frame whose horizon runs past the end of the record
    # cannot be called negative: the crossing that would have made it positive
    # may lie in the part of the trajectory that was never written.
    full_window = times_ps + tau_ps <= t_last_ps
    # A confirmed positive needs no complete window — an observed completion
    # inside the horizon settles the frame however little record follows it.
    valid = full_window | (y == 1)
    # Entrance anchor, same construction on the time an ion enters the pore
    # instead of the time it leaves. With the exit anchor a positive frame often
    # already holds the crossing ion inside the lumen, so occupancy descriptors
    # can restate the label; anchoring on the entrance puts the horizon before
    # the ion is in. Only entrances that were actually observed reach here — the
    # caller drops left-censored events, whose entry time is unknown rather than
    # early.
    y_ent = np.zeros(len(times_ps), dtype=np.int8)
    for te in entry_times:
        y_ent |= ((times_ps < te) & (times_ps >= te - tau_ps)).astype(np.int8)
    valid_ent = full_window | (y_ent == 1)
    # τ is stamped into every column name as a whole number of ps; the training,
    # evaluation and transfer steps look these columns up by that name.
    return pd.DataFrame({
        f"y_tau{int(tau_ps)}": y,
        f"valid_tau{int(tau_ps)}": valid,
        f"y_entr_tau{int(tau_ps)}": y_ent,
        f"valid_entr_tau{int(tau_ps)}": valid_ent,
    })


def frame_grid(cfg: Config, rep: dict) -> np.ndarray:
    """Frame times of the run at the declared stride, in ps, ascending.

    Times are read from the trajectory rather than reconstructed from a saving
    interval, and the file is strided from frame 0 exactly as the features step
    strides it, so the two steps produce the same grid. The training step merges
    labels onto features by (condition, replica, time_ps) with a one-to-one
    check, which is where any divergence of the grids would surface.
    """
    u = pio.open_replica(rep)
    stride = cfg["data.stride"]
    idx = range(0, u.trajectory.n_frames, stride)
    times = np.array([float(u.trajectory[i].time) for i in idx])
    return times


def run_step(cfg: Config) -> None:
    """Labels step: one row per frame of every replica, at every horizon in the sweep.

    Writes labels.parquet (the label and validity columns plus the in_transit
    flag) and summary.json, which holds the event counts, the number of censored
    frames dropped and the base rate on the valid frames for each replica and
    horizon. The horizon sweep is mandatory: a conclusion that holds only at one
    τ is a property of that τ.
    """
    ev = load_events(cfg)
    taus = cfg["labels.tau_ps"]
    with step_output(cfg, "labels") as out:
        log = StepLog(out)
        frames_list = []
        summary: dict = {"per_replica": {}, "tau_ps": taus}
        for cond in cfg["data.conditions"]:
            direction = cond.get("direction")
            for rep in cond["replicas"]:
                key = f"{cond['id']}/{rep['id']}"
                times = frame_grid(cfg, rep)
                # The horizon is judged against the last frame actually on the
                # grid, not against a nominal simulation length.
                t_last = float(times[-1])
                sel = ev[(ev["condition"] == cond["id"]) & (ev["replica"] == rep["id"])]
                # The label is defined by crossings in the condition's
                # conduction direction.
                # A passage against the applied field is a different event, and
                # counting it would make this an "any passage" label rather than a
                # conduction one. In the shipped runs the filter removes nothing:
                # every recorded crossing goes with the field.
                # direction = 0 marks an event whose sign was never recorded — a
                # provided annotation read while the condition's direction was still
                # undeclared. Those are kept rather than discarded by a filter that
                # cannot classify them.
                if direction is not None:
                    sel_dir = sel[(sel["direction"] == direction) | (sel["direction"] == 0)]
                else:
                    sel_dir = sel
                # Left-censored events: the ion was already inside on the first
                # frame, so its entrance was never seen. Substituting the exit time
                # would invent an entrance; the count is reported instead.
                n_dropped_entry = int((~sel_dir["entrance_observed"]).sum()) if len(sel_dir) else 0
                entry_times = sel_dir.loc[sel_dir["entrance_observed"], "t_entry_ps"].to_numpy()
                exit_times = sel_dir["t_exit_ps"].to_numpy()
                base = pd.DataFrame({"condition": cond["id"], "replica": rep["id"],
                                     "time_ps": times})
                parts = [base]
                rep_sum: dict = {"n_events_dir": int(len(sel_dir)),
                                 "n_events_total": int(len(sel)),
                                 "n_entry_unobserved_excluded_from_entrance": n_dropped_entry}
                for tau in taus:
                    lab = label_one_trajectory(times, exit_times, entry_times, tau, t_last)
                    parts.append(lab)
                    # The base rate is taken over the valid frames only, the same
                    # denominator the metrics use later. Averaging over all frames
                    # would count the censored tail as negative and quote a rate
                    # below the one the model is actually scored against.
                    rep_sum[f"tau{int(tau)}"] = {
                        "positives": int(lab[f"y_tau{int(tau)}"].sum()),
                        "censored_dropped": int((~lab[f"valid_tau{int(tau)}"]).sum()),
                        "base_rate_valid": float(
                            lab.loc[lab[f"valid_tau{int(tau)}"], f"y_tau{int(tau)}"].mean())
                        if lab[f"valid_tau{int(tau)}"].any() else float("nan"),
                    }
                # Transit frames (permeant already inside) — for the
                # entrance-anchor ablation.
                # A frame counts as transit if it falls between an observed entrance
                # and its exit, inclusive; only events with a seen entrance define an
                # interval at all. Occupancy of the lumen does not depend on
                # which way
                # the ion is going, so this uses all events of the replica
                # rather than
                # the direction-filtered subset. The evaluation step scores the
                # resting
                # and the transit frames separately: on a transit frame the
                # descriptors
                # can see the crossing ion itself, and the study's claim has to
                # hold on
                # the frames where they cannot.
                transit = np.zeros(len(times), dtype=bool)
                obs = sel[sel["entrance_observed"]]
                for te, tx in zip(obs["t_entry_ps"], obs["t_exit_ps"]):
                    transit |= (times >= te) & (times <= tx)
                # base is already inside parts, so the flag travels with the block.
                base["in_transit"] = transit
                frames_list.append(pd.concat(parts, axis=1))
                summary["per_replica"][key] = rep_sum
                log.say(f"[{key}] events in direction: {rep_sum['n_events_dir']}; "
                        + "; ".join(
                            f"τ={int(t)}: +{rep_sum[f'tau{int(t)}']['positives']}, "
                            f"censored dropped {rep_sum[f'tau{int(t)}']['censored_dropped']}, "
                            f"base rate {rep_sum[f'tau{int(t)}']['base_rate_valid']:.3f}"
                            for t in taus)
                        + f"; entrance labels exclude events without entry: {n_dropped_entry}")
        labels = pd.concat(frames_list, ignore_index=True)
        labels.to_parquet(out / "labels.parquet")
        (out / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
        log.close()


def load_labels(cfg: Config) -> pd.DataFrame:
    """Read the labels artifact of this run; refuse if the step has not been run.

    All replicas in one table, keyed by (condition, replica, time_ps), with a
    y/valid pair per horizon named after τ. The error names the step that has to
    be run first, which is otherwise reported as a missing file somewhere inside
    the training step.
    """
    p = run_dir(cfg) / "labels" / "labels.parquet"
    if not p.exists():
        raise FileNotFoundError(f"{p} missing: the labels step has not run")
    return pd.read_parquet(p)
