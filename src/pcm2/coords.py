"""Superimposed C-alpha coordinates of the channel, exported per frame.

This artifact exists for the PLS-FMA benchmark arm (Krivobokova et al. 2012,
Biophys. J. 103:786-796, doi:10.1016/j.bpj.2012.07.022): that method regresses a
functional quantity on least-squares-superimposed Cartesian coordinates, which
is a different input from the descriptor table. Keeping the export here, in the
trajectory-reading layer, preserves the architectural rule that the model layer
never opens a trajectory: models consume this table like any other artifact.

Superposition: for every frame, channel C-alpha positions are centred on their
geometric mean and rotated onto the reference frame (frame 0 of the first
replica) with the Kabsch algorithm. Rigid-body motion is thereby removed, so
the exported columns describe internal conformational motion only.
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
from MDAnalysis.analysis.align import rotation_matrix

from . import io as pio
from .config import Config, load as load_config
from .runtime import StepLog, run_dir, step_output


def reference_ca(cfg: Config) -> np.ndarray:
    """Reference structure for the superposition: frame 0 of the first replica, centred.

    → [n_CA, 3] in A, with the geometric mean at the origin. One reference for the whole
    run, so coordinates from different replicas and conditions live in the same frame and
    can be regressed together. Which structure serves as the reference is arbitrary; that
    it is the same one for every frame is not, since a reference changing between replicas
    would enter the exported columns as conformational variance.
    """
    first = cfg["data.conditions"][0]["replicas"][0]
    u = pio.open_replica(first)
    ca = u.select_atoms(f"({cfg['data.selections.channel']}) and name CA")
    u.trajectory[0]
    ref = ca.positions.copy()
    return ref - ref.mean(axis=0)


def _worker(cfg_path: str, cond_id: str, rep_id: str, ref: np.ndarray) -> pd.DataFrame:
    """One replica, superimposed frame by frame → a table of Cartesian components.

    Runs in its own process, so it reopens the config and the trajectory rather than
    receiving them. Columns are named ca{resid}_{x,y,z} and hold coordinates in A; the
    leading three (condition, replica, time_ps) are the key on which the model layer joins
    this table to the descriptor table.
    """
    cfg = load_config(cfg_path)
    cond = next(c for c in cfg["data.conditions"] if c["id"] == cond_id)
    rep = next(r for r in cond["replicas"] if r["id"] == rep_id)
    u = pio.open_replica(rep)
    ca = u.select_atoms(f"({cfg['data.selections.channel']}) and name CA")
    # A superposition onto a reference of a different size is not defined, and a
    # table
    # whose rows describe different atoms would still look like a valid matrix
    # downstream.
    if len(ca) != len(ref):
        raise RuntimeError(f"{cond_id}/{rep_id}: C-alpha count {len(ca)} differs from "
                           f"the reference ({len(ref)}) — replicas are not comparable")
    # Same stride as the descriptor table: the two artifacts are joined on
    # time_ps, so a
    # frame exported here without a descriptor row would be dropped by the join
    # anyway.
    stride = cfg["data.stride"]
    frames = list(range(0, u.trajectory.n_frames, stride))
    rows = np.empty((len(frames), ref.size), dtype=np.float32)
    times = np.empty(len(frames))
    for k, fi in enumerate(frames):
        ts = u.trajectory[fi]
        times[k] = ts.time
        # Kabsch superposition: translate to the geometric mean (unweighted, so every
        # C-alpha counts the same), then rotate onto the reference.
        # rotation_matrix returns
        # R such that pos @ R.T is the rotated set, and the RMSD it also returns
        # is not
        # kept — the residual is a property of the fit, not of the conformation.
        pos = ca.positions - ca.positions.mean(axis=0)
        rot, _rmsd = rotation_matrix(pos, ref)
        rows[k] = (pos @ rot.T).astype(np.float32).ravel()
    # Column order follows the flattening above: resid-major, then x, y, z.
    # Resids differ
    # between the KcsA mutant topologies, so concatenating conditions can yield
    # a union of
    # column names; the model layer drops the columns not shared by all conditions.
    cols = [f"ca{ix}_{ax}" for ix in ca.resids for ax in "xyz"]
    df = pd.DataFrame(rows, columns=cols)
    df.insert(0, "condition", cond_id)
    df.insert(1, "replica", rep_id)
    df.insert(2, "time_ps", times)
    return df


def run_step(cfg: Config) -> None:
    """Write coords/coords.parquet: every frame of every replica, superimposed.

    The reference is built once in this process and handed to the workers, so all of them
    align to the same structure. The pool is capped by the number of replicas and by
    PCM2_MAX_WORKERS; each worker holds one open trajectory and one table of its own.
    """
    from concurrent.futures import ProcessPoolExecutor

    ref = reference_ca(cfg)
    with step_output(cfg, "coords") as out:
        log = StepLog(out)
        jobs = [(c["id"], r["id"]) for c in cfg["data.conditions"] for r in c["replicas"]]
        dfs = []
        with ProcessPoolExecutor(max_workers=min(len(jobs), int(os.environ.get("PCM2_MAX_WORKERS", "14")))) as ex:
            futs = {ex.submit(_worker, str(cfg.source_path), ci, ri, ref): (ci, ri)
                    for ci, ri in jobs}
            for fut, key in futs.items():
                dfs.append(fut.result())
                log.say(f"[{key[0]}/{key[1]}] {len(dfs[-1])} frames aligned")
        table = pd.concat(dfs, ignore_index=True)
        table.to_parquet(out / "coords.parquet")
        log.say(f"superimposed coordinates: {len(table)} frames x "
                f"{table.shape[1] - 3} Cartesian components")
        log.close()


def load_coords(cfg: Config) -> pd.DataFrame:
    """The coords artifact of this run, or a refusal naming the step that produces it."""
    p = run_dir(cfg) / "coords" / "coords.parquet"
    if not p.exists():
        raise FileNotFoundError(f"{p} is missing: run the coords step first")
    return pd.read_parquet(p)
