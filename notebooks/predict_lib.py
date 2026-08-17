"""Single-point and range prediction on top of the trained pooled monitor.

Three pieces used by the interactive notebook:
  * a pooled model trained on all conducting systems (cached on first use);
  * feature rows for a requested time or time range, read from the system's
    computed table when one exists, otherwise measured directly from the
    trajectory over the backward window the features need;
  * an explanation of each verdict as interpretation axes with weights.

Nothing here modifies the pipeline; every measurement path is imported from
the pipeline modules so a notebook prediction and a pipeline table agree.
"""

from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

PROJ = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJ / "src"))
sys.path.insert(0, str(PROJ / "tools"))
os.chdir(PROJ)  # system configs address their files relative to the project root

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

from pcm2 import config  # noqa: E402
from pcm2 import io as pio  # noqa: E402
from pcm2.features import BLOCK_MODULES, ReplicaComputer, build_schema  # noqa: E402
from pcm2.features._common import ColSpec  # noqa: E402
from pcm2.evaluate import ap_with_base  # noqa: E402
from pcm2.interpret import AXES  # noqa: E402
from pcm2.labels import label_one_trajectory  # noqa: E402
from pcm2.models import Calibrator, f1_threshold  # noqa: E402
from pcm2.models.trees import fit_trees_fold, trees_contributions, trees_predict  # noqa: E402
from pcm2.runtime import run_dir  # noqa: E402

from monitor_generalisation import SYSTEM_CONFIGS, is_conducting, load_all  # noqa: E402

POOLED_SEED = 20260816  # the pooled fit-on-all-conducting experiment seed
CACHE_DIR = PROJ / "runs" / "pooled-predict"
CACHE_VERSION = 3  # bumped when the stored bundle gains fields or the
                   # descriptor schema changes (v3: dimensionless columns)


# --------------------------------------------------------------------------
# Config resolution: a data folder points back to the system that owns it.
# --------------------------------------------------------------------------

def resolve_config(protein: str | Path) -> config.Config:
    """Config for a protein, by system id or by its data folder.

    A bare name ("gramicidin", "cx43") loads configs/<name>.yaml directly and
    is portable across machines; a filesystem path is matched against the
    data_root recorded in each config.
    """
    as_id = PROJ / "configs" / f"{protein}.yaml"
    if isinstance(protein, str) and "/" not in protein and as_id.exists():
        return config.load(as_id)
    protein_dir = Path(protein).resolve()
    candidates = []
    for p in sorted((PROJ / "configs").glob("*.yaml")):
        doc = yaml.safe_load(p.read_text())
        root = Path(doc.get("data", {}).get("data_root", "")).resolve()
        if root == protein_dir or root.is_relative_to(protein_dir):
            candidates.append(p)
    if not candidates:
        known = [str(p.name) for p in sorted((PROJ / "configs").glob("*.yaml"))]
        raise FileNotFoundError(
            f"no system config points at {protein_dir}; known configs: {known}. "
            "A new protein needs its config assembled and the autodetect step "
            "run once (tools/init_configs.py, then `python -m pcm2.run "
            "autodetect --config configs/<system>.yaml`).")
    if len(candidates) > 1:
        raise RuntimeError(f"ambiguous folder: {candidates}")
    return config.load(candidates[0])


def resolve_replica(cfg: config.Config, condition: str | None,
                    replica: str | None) -> tuple[dict, dict]:
    conds = cfg["data.conditions"]
    cond = (next((c for c in conds if c["id"] == condition), None)
            if condition else conds[0])
    if cond is None:
        raise KeyError(f"condition {condition!r} not in "
                       f"{[c['id'] for c in conds]}")
    rep = (next((r for r in cond["replicas"] if r["id"] == replica), None)
           if replica else cond["replicas"][0])
    if rep is None:
        raise KeyError(f"replica {replica!r} not in "
                       f"{[r['id'] for r in cond['replicas']]}")
    return cond, rep


# --------------------------------------------------------------------------
# Pooled model: trained once on every conducting system's rows, then cached.
# --------------------------------------------------------------------------

def _physics_only_columns(pool: pd.DataFrame, model_cols: list[str],
                          conducting: np.ndarray) -> list[str]:
    """Drop columns whose per-system ranges are disjoint (system labels)."""
    sys_label = []
    cpool = pool[conducting]
    for c in model_cols:
        ranges = [(g[c].min(), g[c].max()) for _s, g in cpool.groupby("system")
                  if g[c].notna().any()]
        if len(ranges) >= 2 and max(r[0] for r in ranges) > min(r[1] for r in ranges):
            sys_label.append(c)
    return [c for c in model_cols if c not in sys_label]


def relabel_pool(pool: pd.DataFrame, tau_ps: float) -> tuple[np.ndarray, np.ndarray]:
    """Labels at an arbitrary horizon, rebuilt from the recorded crossings.

    Same rule the pipeline uses: a frame is positive when a crossing in the
    condition's conduction direction completes within the next tau ps, and a
    frame whose window runs past the end of the trajectory without a crossing
    is censored rather than called negative.
    """
    ycol, vcol = f"y_tau{int(tau_ps)}", f"valid_tau{int(tau_ps)}"
    y = np.zeros(len(pool), dtype=int)
    valid = np.zeros(len(pool), dtype=bool)
    for path in SYSTEM_CONFIGS:
        cfg = config.load(path)
        sysid = cfg["system.id"]
        ev = pd.read_parquet(run_dir(cfg) / "events" / "events.parquet")
        for cond in cfg["data.conditions"]:
            direction = cond.get("direction")
            for rep in cond["replicas"]:
                m = ((pool["system"] == sysid) & (pool["condition"] == cond["id"])
                     & (pool["replica"] == rep["id"])).to_numpy()
                if not m.any():
                    continue
                times = pool.loc[m, "time_ps"].to_numpy()
                sel = ev[(ev["condition"] == cond["id"])
                         & (ev["replica"] == rep["id"])]
                if direction is not None:
                    sel = sel[(sel["direction"] == direction) | (sel["direction"] == 0)]
                lab = label_one_trajectory(
                    times, sel["t_exit_ps"].to_numpy(),
                    sel.loc[sel["entrance_observed"], "t_entry_ps"].to_numpy(),
                    tau_ps, float(times[-1]))
                idx = np.flatnonzero(m)
                y[idx] = lab[ycol].to_numpy()
                valid[idx] = lab[vcol].to_numpy()
    return y, valid


def resolve_labels(pool: pd.DataFrame, tau_ps: float,
                   relabel=None) -> tuple[np.ndarray, np.ndarray, str]:
    """Labels at this horizon, from the tables when they are complete.

    A horizon declared by only part of the systems gives a column that is
    missing on the rest of the rows, and missing rows would be read as
    negatives. Only a column complete over every row is used as it stands;
    otherwise the labels are rebuilt from the recorded crossings.
    """
    relabel = relabel or relabel_pool
    ycol, vcol = f"y_tau{int(tau_ps)}", f"valid_tau{int(tau_ps)}"
    present = ycol in pool.columns and vcol in pool.columns
    complete = present and bool(pool[ycol].notna().all()) and bool(pool[vcol].notna().all())
    if complete:
        y = pool[ycol].to_numpy(int)
        valid = pool[vcol].to_numpy(bool)
        source = "pipeline tables"
    else:
        y, valid = relabel(pool, tau_ps)
        source = ("relabelled from events (horizon declared by part of the systems)"
                  if present else "relabelled from events")
    unique = np.unique(y)
    if not np.isin(unique, (0, 1)).all():
        raise ValueError(f"labels at horizon {tau_ps:.0f} ps are not binary "
                         f"(values {unique[:5]}) — refusing to train on them")
    return y, valid.astype(bool), source


def load_pooled_model(tau_ps: float = 1000.0, variant: str = "physics_only",
                      refresh: bool = False) -> dict:
    """Train (or load from cache) the pooled monitor at the given horizon.

    variant "full" uses every model column; "physics_only" removes columns
    whose per-system value ranges do not overlap, so the model cannot lean on
    system identity.
    """
    cache = CACHE_DIR / f"tau{int(tau_ps)}-{variant}-v{CACHE_VERSION}.pkl"
    if cache.exists() and not refresh:
        with open(cache, "rb") as f:
            return pickle.load(f)

    pool, schema = load_all()
    ycol, vcol = f"y_tau{int(tau_ps)}", f"valid_tau{int(tau_ps)}"
    model_cols = [c.name for c in schema if c.to_model]
    y, valid, labels_from = resolve_labels(pool, tau_ps)
    X = pool[model_cols].to_numpy(np.float32)
    units = pool["unit"].to_numpy()
    conducting = np.array([is_conducting(s, c)
                           for s, c in zip(pool["system"], pool["condition"])])

    cols = (model_cols if variant == "full"
            else _physics_only_columns(pool, model_cols, conducting))
    idx = [model_cols.index(c) for c in cols]
    tr = np.flatnonzero(conducting & valid)

    r = fit_trees_fold(X[np.ix_(tr, idx)], y[tr], units[tr], units[tr], cols,
                       [2, 3, 4], [2.0, 4.0, 8.0], [0.7], [0.8],
                       lr=0.05, ceiling=2000, patience=50, monotone_map={},
                       seed=POOLED_SEED, n_threads=1)
    cal = Calibrator("sigmoid", 200).fit(r.inner_oof_scores, r.inner_oof_y)
    oof_prob = cal.predict(r.inner_oof_scores)
    thr = f1_threshold(oof_prob, r.inner_oof_y)
    # Held-out quality: the inner rotation leaves one trajectory unit out at a
    # time, so no score comes from a fit that included that unit.
    quality = ap_with_base(r.inner_oof_y, oof_prob)

    bundle = {
        "booster": r.booster, "calibrator": cal, "threshold": float(thr),
        "columns": cols, "variant": variant, "tau_ps": float(tau_ps),
        "block_of": {c.name: c.block for c in schema},
        "trained_on": sorted(set(pool.loc[conducting & valid, "unit"])),
        "hyperparams": {"chosen": r.chosen, "M": r.M}, "seed": POOLED_SEED,
        "labels_from": labels_from,
        "train_positives": int(y[tr].sum()), "train_rows": int(len(tr)),
        "base_rate": float(y[tr].mean()),
        "oof_quality": quality,
    }
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    with open(cache, "wb") as f:
        pickle.dump(bundle, f)
    return bundle


# --------------------------------------------------------------------------
# Feature rows for a query: computed table when available, else measured.
# --------------------------------------------------------------------------

def _exit_times(cfg, cond, rep) -> np.ndarray:
    """Directional crossing exit times for this replica.

    The delivery block is built from these times, so the events step is
    required: without it time-since-last-crossing reads as never-crossed on
    every frame.
    """
    path = run_dir(cfg) / "events" / "events.parquet"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} is missing: run `python -m pcm2.run events --config "
            f"configs/{cfg['system.id']}.yaml` first. Without it the delivery "
            "descriptors (time since the last crossing, arrival rate) would be "
            "silently wrong rather than missing.")
    ev = pd.read_parquet(path)
    sel = ev[(ev["condition"] == cond["id"]) & (ev["replica"] == rep["id"])]
    d = cond.get("direction")
    if d is not None and "direction" in sel.columns:
        sel = sel[(sel["direction"] == d) | (sel["direction"] == 0)]
    return sel["t_exit_ps"].to_numpy()


def rows_from_table(cfg, cond, rep, t_lo_ps: float,
                    t_hi_ps: float) -> pd.DataFrame | None:
    """Rows straight out of the system's computed feature table, if it exists."""
    path = run_dir(cfg) / "features" / "features.parquet"
    if not path.exists():
        return None
    t = pd.read_parquet(path)
    t = t[(t["condition"] == cond["id"]) & (t["replica"] == rep["id"])]
    out = t[(t["time_ps"] >= t_lo_ps) & (t["time_ps"] <= t_hi_ps)]
    return out.reset_index(drop=True) if len(out) else None


def rows_from_trajectory(cfg, cond, rep, t_lo_ps: float,
                         t_hi_ps: float) -> pd.DataFrame:
    """Measure feature rows for [t_lo, t_hi] directly from the trajectory.

    Only the frames inside the interval plus one backward window before it are
    read, on the same frame grid the pipeline uses. Two caveats: run-length
    ages (dyn_wet_state_age_ps, dyn_occ_change_age_ps) restart at the window
    edge and are capped by the window length, and the published Rao baseline
    columns (to_model=False, never fed to the model) are left unfilled here.
    """
    rc = ReplicaComputer(cfg, cond, rep)
    ctx = rc.ctx
    ctx.params["exit_times_dir"] = _exit_times(cfg, cond, rep)
    ctx.params["direction"] = rc.direction
    schema = build_schema(cfg)

    stride = cfg["data.stride"]
    window = float(cfg.get("features.window_ps") or 0.0)
    t0, dt = pio.frame_times_ps(rc.u)
    grid = np.arange(0, rc.u.trajectory.n_frames, stride)
    grid_t = t0 + grid * dt
    keep = (grid_t >= t_lo_ps - window - 0.5 * dt) & (grid_t <= t_hi_ps + 0.5 * dt)
    frames = grid[keep]
    if len(frames) == 0:
        raise ValueError(f"no grid frames in [{t_lo_ps}, {t_hi_ps}] ps "
                         f"(trajectory spans {grid_t[0]:.0f}–{grid_t[-1]:.0f} ps)")

    T, N = len(frames), len(rc.ions)
    cols = {c.name: np.full(T, np.nan) for c in schema}
    times = np.empty(T)
    ion_z_all = np.full((T, N), np.nan)
    ion_in_pore = np.zeros((T, N), bool)
    innermost = np.full(T, -1, int)
    active = [BLOCK_MODULES[b] for b in ctx.blocks]

    for k, fi in enumerate(frames):
        ts = rc.u.trajectory[int(fi)]
        fd = rc.measure_frame(ts)
        times[k] = fd.time_ps
        ion_z_all[k] = fd.ion_z_rel
        ion_in_pore[k] = fd.ion_in_pore
        innermost[k] = fd.innermost if fd.innermost is not None else -1
        for m in active:
            if hasattr(m, "compute"):
                m.compute(ctx, fd, k, cols)

    ctx.params["_aux_ion_z_all"] = ion_z_all
    ctx.params["_aux_ion_in_pore"] = ion_in_pore
    ctx.params["_aux_innermost"] = innermost
    for m in active:
        if hasattr(m, "post"):
            m.post(ctx, times, cols)

    if rc.direction is not None and rc.direction < 0:
        for c in schema:
            if c.sign_flip:
                cols[c.name] = -cols[c.name]

    df = pd.DataFrame({"condition": cond["id"], "replica": rep["id"],
                       "time_ps": times})
    for c in schema:
        df[c.name] = cols[c.name]
    return df[(df["time_ps"] >= t_lo_ps - 0.5 * dt)
              & (df["time_ps"] <= t_hi_ps + 0.5 * dt)].reset_index(drop=True)


def row_at(cfg, condition: str | None, replica: str | None,
           time_ps: float) -> tuple[pd.DataFrame, dict, str, float]:
    """The measured frame nearest to time_ps, with the offset from it.

    Requested times fall between strided frames, so an exact time filter would
    usually return nothing.
    """
    cond, rep = resolve_replica(cfg, condition, replica)
    rows = rows_from_table(cfg, cond, rep, -np.inf, np.inf)
    source = "computed table"
    if rows is None:
        u = pio.open_replica(rep)
        _, dt = pio.frame_times_ps(u)
        step = dt * cfg["data.stride"]
        rows = rows_from_trajectory(cfg, cond, rep, time_ps - step, time_ps + step)
        source = "measured from trajectory (window mode)"
    i = int((rows["time_ps"] - time_ps).abs().idxmin())
    row = rows.loc[[i]].reset_index(drop=True)
    return row, cond, source, float(abs(row["time_ps"].iloc[0] - time_ps))


def query_rows(cfg, condition: str | None, replica: str | None,
               t_lo_ps: float, t_hi_ps: float) -> tuple[pd.DataFrame, dict, str]:
    cond, rep = resolve_replica(cfg, condition, replica)
    rows = rows_from_table(cfg, cond, rep, t_lo_ps, t_hi_ps)
    source = "computed table"
    if rows is None:
        rows = rows_from_trajectory(cfg, cond, rep, t_lo_ps, t_hi_ps)
        source = "measured from trajectory (window mode)"
    return rows, cond, source


# --------------------------------------------------------------------------
# Prediction and per-verdict explanation.
# --------------------------------------------------------------------------

def predict_rows(bundle: dict, rows: pd.DataFrame) -> pd.DataFrame:
    missing = [c for c in bundle["columns"] if c not in rows.columns]
    if missing:
        raise KeyError(f"feature rows lack model columns: {missing[:5]} …")
    X = rows[bundle["columns"]].to_numpy(np.float32)
    raw = trees_predict(bundle["booster"], X, bundle["columns"])
    prob = bundle["calibrator"].predict(raw)
    out = rows[["condition", "replica", "time_ps"]].copy()
    out["p_ready"] = prob
    out["ready"] = prob >= bundle["threshold"]
    return out


def explain_rows(bundle: dict, rows: pd.DataFrame,
                 top_k: int = 3) -> list[list[dict]]:
    """Top axes behind each verdict, with normalized weights.

    Contributions are TreeSHAP values in log-odds; a column's missing value
    contributes what the tree path assigns to its NaN branch. Columns are
    grouped into interpretation axes by their feature block; each axis gets
    the signed sum of its columns' contributions. For a READY verdict the
    reasons are the axes pushing the score up, for NOT_READY — down; the
    weight is the axis share of the total push in that direction.
    """
    cols = bundle["columns"]
    block_of = bundle["block_of"]
    axis_cols = {ax: [c for c in cols if block_of.get(c) in blocks]
                 for ax, blocks in AXES.items()}
    axis_cols = {ax: cc for ax, cc in axis_cols.items() if cc}

    X = rows[cols].to_numpy(np.float32)
    raw = trees_predict(bundle["booster"], X, cols)
    prob = bundle["calibrator"].predict(raw)
    contrib = trees_contributions(bundle["booster"], X, cols)

    reports = []
    for i in range(len(rows)):
        ready = prob[i] >= bundle["threshold"]
        sign = 1.0 if ready else -1.0
        axes_val = {}
        for ax, cc in axis_cols.items():
            sub = contrib[i, [cols.index(c) for c in cc]]
            axes_val[ax] = float(np.nansum(np.where(np.isnan(sub), 0.0, sub)))
        pushed = {ax: v for ax, v in axes_val.items() if sign * v > 0}
        total = sum(abs(v) for v in pushed.values()) or 1.0
        ranked = sorted(pushed.items(), key=lambda kv: -abs(kv[1]))[:top_k]
        report = []
        for ax, v in ranked:
            cc = axis_cols[ax]
            sub = contrib[i, [cols.index(c) for c in cc]]
            # sign first, then the -inf floor, so a missing contribution cannot
            # win the argmax on a NOT-READY verdict.
            j = int(np.nanargmax(np.where(np.isnan(sub), -np.inf, sign * sub)))
            col = cc[j]
            val = rows.iloc[i][col]
            report.append({
                "axis": ax, "weight": abs(v) / total,
                "logodds": v, "top_column": col,
                "column_value": None if pd.isna(val) else float(val),
            })
        reports.append(report)
    return reports


# --------------------------------------------------------------------------
# Snapshot: the measured state of one moment, no horizon required.
# --------------------------------------------------------------------------

SNAPSHOT_QUANTITIES = [
    ("pore geometry", [
        ("geo_r_constriction_A", "constriction radius", "Å"),
        ("geo_z_constriction_A", "constriction position on the axis", "Å"),
        ("hyd_lining_hydrophobicity", "lining hydrophobicity at the constriction", "WW scale"),
    ]),
    ("hydration", [
        ("hyd_water_pore_n", "water molecules in the pore", ""),
        ("hyd_wet_frac_win", "wet fraction of the recent window", ""),
        ("ww_continuous", "water wire continuous (1 = yes)", ""),
        ("ww_max_gap_A", "largest gap in the water wire", "Å"),
    ]),
    ("ions", [
        ("occ_n_ions_pore", "permeant ions in the pore", ""),
        ("occ_innermost_z_A", "innermost ion position on the axis", "Å"),
        ("occ_desolvation", "innermost ion desolvation", ""),
        ("dlv_n_entry_mouth", "permeant ions at the entry mouth", ""),
        ("dlv_dist_entry_A", "nearest ion distance to the entry", "Å"),
        ("dlv_t_since_cross_ps", "time since the last crossing", "ps"),
    ]),
    ("published criterion", [
        ("bl_rao_E_kJmol", "Rao-2019 dewetting energy (>= 2.6 = dewetted)", "kJ/mol"),
    ]),
]


def snapshot(cfg, condition: str | None, replica: str | None,
             time_ps: float) -> tuple[pd.Series, dict, float]:
    """The nearest measured frame to the requested time, as a state readout."""
    cond, rep = resolve_replica(cfg, condition, replica)
    rows = rows_from_table(cfg, cond, rep, -np.inf, np.inf)
    if rows is None:
        u = pio.open_replica(rep)
        t0, dt = pio.frame_times_ps(u)
        step = dt * cfg["data.stride"]
        rows = rows_from_trajectory(cfg, cond, rep, time_ps - step, time_ps + step)
    i = int((rows["time_ps"] - time_ps).abs().idxmin())
    row = rows.loc[i]
    return row, cond, float(abs(row["time_ps"] - time_ps))


def snapshot_text(row: pd.Series, tau_probs: dict[float, tuple[float, bool]]) -> str:
    t_ns = row["time_ps"] / 1000.0
    lines = [f"state at t = {row['time_ps']:.0f} ps ({t_ns:.1f} ns), "
             f"{row['condition']}/{row['replica']}"]
    for group, quants in SNAPSHOT_QUANTITIES:
        shown = []
        for col, label, unit in quants:
            if col not in row.index:
                continue
            v = row[col]
            txt = "not measured on this frame" if pd.isna(v) else f"{v:.3g}"
            shown.append(f"    {label}: {txt}{(' ' + unit) if unit and not pd.isna(v) else ''}")
        if shown:
            lines.append(f"  {group}:")
            lines.extend(shown)
    if tau_probs:
        lines.append("  readiness at this moment:")
        for tau, (p, ready) in sorted(tau_probs.items()):
            state = "READY" if ready else "NOT READY"
            lines.append(f"    within {tau:.0f} ps: P = {p:.3f} -> {state}")
    return "\n".join(lines)


def cached_horizons(variant: str = "physics_only") -> list[float]:
    """Horizons for which a pooled model is already trained and cached."""
    out = []
    for p in sorted(CACHE_DIR.glob(f"tau*-{variant}-v{CACHE_VERSION}.pkl")):
        out.append(float(p.name.split("-")[0].removeprefix("tau")))
    return out


def verdict_text(pred_row: pd.Series, reasons: list[dict],
                 tau_ps: float) -> str:
    state = "READY" if pred_row["ready"] else "NOT READY"
    lines = [f"t = {pred_row['time_ps']:.0f} ps | horizon τ = {tau_ps:.0f} ps | "
             f"P(crossing within τ) = {pred_row['p_ready']:.3f} → {state}"]
    for k, r in enumerate(reasons, 1):
        v = ("missing" if r["column_value"] is None
             else f"{r['column_value']:.3g}")
        lines.append(f"  {k}. {r['axis']} — weight {r['weight']:.0%} "
                     f"(driver: {r['top_column']} = {v})")
    return "\n".join(lines)
