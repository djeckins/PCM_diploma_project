"""Nested CV, calibration, thresholding, arms, and the final model.

The head model is declared IN ADVANCE in the config (model.head) — no selection
over outer folds is performed, and outer-fold metrics play no part in the choice.

The nesting has two levels and one direction of information flow. The outer level
leaves one independent assembly out; inside a fold, a rotation over the remaining
assemblies selects the hyperparameters, produces the out-of-fold scores the
calibrator and the decision threshold are fitted on, and only then is the held-out
assembly scored. Nothing measured on the held-out assembly re-enters the fit, which
is why the calibrator and the threshold are built per fold rather than once at the
end.

Every arm goes through the same folds, the same calibration and the same artifact
layout, so the arms differ only in what they are allowed to look at: all modelled
columns, two structural columns, three temporal columns, superimposed Cartesian
coordinates, or a published table with no fitting at all.
"""

from __future__ import annotations

import json
import pickle
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import precision_recall_curve

from .. import evaluate as evaluate_mod
from .. import interpret as interpret_mod
from ..baselines.clock import clock_columns
from ..baselines.linear_control import control_columns
from ..config import Config
from ..config import load as load_config
from ..datasets import build_matrix
from ..features import load_features
from ..labels import load_labels
from ..runtime import StepLog, code_fingerprint, derive_seed, pin_threads, run_dir, step_output
from ..splits import check_no_frame_leak, outer_folds
from .linear import fit_linear_fold, linear_contributions
from .plsfma import fit_plsfma_fold, plsfma_predict
from .trees import fit_trees_fold, trees_contributions, trees_predict

# Arms that are fitted here, and arms that only read a column of the table.
# The published arms are listed separately because they are never trained: they
# have no fold, no calibrator and no probability scale.
ML_ARMS = ("trees", "linear", "linear_control_structural", "clock")
PUBLISHED_ARMS = ("published_rao2019", "published_rao2019_win")


def _logit(p):
    # Clipped away from 0 and 1: an exact 0 or 1 probability would map to
    # infinity, and the sigmoid calibrator is fitted on this scale. 1e-7 is far
    # below any probability the arms resolve on samples of this size.
    p = np.clip(p, 1e-7, 1 - 1e-7)
    return np.log(p / (1 - p))


class Calibrator:
    """Map that is monotone BY CONSTRUCTION: Platt sigmoid or isotonic regression.
    The flexible isotonic fit needs many positives: below the threshold it falls
    back to the sigmoid.

    kind="probability": the score is a probability; the sigmoid is fit on its logit.
    kind="raw": the score is an arbitrary real value (PLS regression); the sigmoid
    is fit on the score itself; monotonicity is preserved in both cases.

    Monotone by construction is the requirement that matters here: calibration
    must not change the ranking of frames, or the ranking metrics and the
    calibrated probabilities would describe two different models. Both families
    have that property; a general regression of y on the score would not.
    """

    def __init__(self, family: str, min_positives_isotonic: int,
                 kind: str = "probability"):
        self.family_requested = family
        self.min_pos = min_positives_isotonic
        self.family = family
        self.kind = kind
        self.degenerate = False

    def _design(self, s: np.ndarray) -> np.ndarray:
        """The single regressor the sigmoid is fitted on.

        A probability enters as its logit, which is Platt's original
        parameterization and makes an already well-calibrated arm a fit with
        slope 1 and intercept 0. A raw score enters unchanged, since it has no
        logit to take.
        """
        return _logit(s) if self.kind == "probability" else s

    def fit(self, scores: np.ndarray, y: np.ndarray) -> "Calibrator":
        """Fit the map on held-out scores and their labels; returns self.

        scores must be out-of-fold: fitting on in-sample scores would learn the
        overconfidence of the training fit and then remove it, which looks like
        calibration and is not. Non-finite scores are dropped here, which is how
        rows that no rotation slot could score leave the calibration set.
        """
        m = np.isfinite(scores)
        s, yy = scores[m], y[m]
        # Single-class calibration data (an all-negative training pool of
        # near-non-conducting units): a probability map cannot be fitted.
        # The arm records this state; raw scores stay usable for ranking and
        # calibrated probabilities are absent.
        if len(np.unique(yy)) < 2:
            self.unfit = True
            self.degenerate = True
            return self
        self.unfit = False
        # Isotonic regression is a free step function: it needs enough positives
        # to place its steps, or it reproduces the noise of the calibration set.
        # The number is model.calibration.min_positives_isotonic.
        if self.family == "isotonic" and yy.sum() < self.min_pos:
            self.family = "sigmoid"  # fallback is recorded in the artifact
        if self.family == "isotonic":
            self.iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
            self.iso.fit(s, yy)
        else:
            # Unpenalized on purpose: the fit has two parameters, and shrinking
            # them would bias the calibrated probabilities toward the base rate.
            self.lr = LogisticRegression(penalty=None, solver="lbfgs", max_iter=1000)
            self.lr.fit(self._design(s)[:, None], yy)
            # A non-positive slope means the calibration set says the score points
            # the wrong way. It is recorded, not corrected: flipping the sign here
            # would hide a broken arm behind a working calibrator.
            self.degenerate = bool(self.lr.coef_[0][0] <= 0)
        return self

    def predict(self, scores: np.ndarray) -> np.ndarray:
        """Calibrated P(y=1) per score; NaN where the score is absent or the map is unfit.

        NaN rather than a substituted value, so a missing probability travels into
        the OOF table as missing and no metric silently treats it as 0.5.
        """
        out = np.full(len(scores), np.nan)
        if getattr(self, "unfit", False):
            return out
        m = np.isfinite(scores)
        if self.family == "isotonic":
            out[m] = self.iso.predict(scores[m])
        else:
            out[m] = self.lr.predict_proba(self._design(scores[m])[:, None])[:, 1]
        return out

    def describe(self) -> dict:
        """What the fold artifact records: family asked for, family used, and the flags.

        The requested and applied families are both kept because the isotonic
        fallback happens per fold, so a run can be part isotonic and part sigmoid.
        """
        return {"requested": self.family_requested,
                "applied": "unfit" if getattr(self, "unfit", False) else self.family,
                "kind": self.kind, "degenerate_fit": self.degenerate}


def f1_threshold(prob: np.ndarray, y: np.ndarray) -> float:
    """Threshold is chosen inside the training fold; the rule is maximum F1.

    Returns the probability above which a frame is called ready. Chosen on the
    inner out-of-fold scores, never on the test rows of the outer fold: a
    threshold tuned on the test side would turn a reported operating point into a
    quantity fitted to the answer. 0.5 is returned when the calibration set has no
    positives and there is nothing to maximize.
    """
    m = np.isfinite(prob)
    if y[m].sum() == 0:
        return 0.5
    prec, rec, thr = precision_recall_curve(y[m], prob[m])
    # 1e-12 guards the point where precision and recall are both zero.
    f1 = 2 * prec * rec / np.maximum(prec + rec, 1e-12)
    # precision_recall_curve returns one more (precision, recall) point than
    # thresholds: the last one is the degenerate (recall 0, precision 1) end of
    # the curve and has no threshold, so it is excluded from the argmax.
    best = int(np.argmax(f1[:-1])) if len(thr) else 0
    return float(thr[best]) if len(thr) else 0.5


def _fit_ml_arm(arm: str, cfg_tree: dict, X: np.ndarray, y: np.ndarray,
                assemblies: np.ndarray, traj: np.ndarray,
                names: list[str], seed: int, times: np.ndarray | None = None,
                block_embargo_ps: float | None = None) -> tuple[object, dict, np.ndarray, np.ndarray]:
    """→ (model, hyperparameters/diagnostics, inner-OOF scores, inner-OOF y).

    arm is the branch, "trees" or "linear", and not the reported arm name: the
    structural control and the clock arm are linear fits on a restricted column
    set, so they arrive here as "linear" with fewer names. Everything downstream
    of this call treats the two branches alike, which is why both return the same
    four things. Every argument is already restricted to the training rows of the
    fold; nothing here knows the test rows exist.
    """
    mc = cfg_tree["model"]
    if arm == "trees":
        r = fit_trees_fold(X, y, assemblies, traj, names,
                           mc["trees"]["depth_grid"], mc["trees"]["k_events_grid"],
                           mc["trees"]["colsample_grid"], mc["trees"]["subsample_grid"],
                           mc["trees"]["learning_rate"], mc["trees"]["rounds_ceiling"],
                           mc["trees"]["patience"], cfg_tree["model"]["monotone"],
                           seed, mc["n_threads"], times_tr=times,
                           block_embargo_ps=block_embargo_ps)
        hp = {"chosen": r.chosen, "M_rounds": r.M, "grid_scores": r.grid_scores,
              "edges": r.edges, "leaf_math": r.leaf_math, "hit_ceiling": r.hit_ceiling,
              "learning_rate_fixed": mc["trees"]["learning_rate"],
              "blocked_inner_rotation": r.blocked_inner}
        return r.booster, hp, r.inner_oof_scores, r.inner_oof_y
    lin = mc["linear"]
    r = fit_linear_fold(X, y, assemblies, names, lin["C_grid"], lin["tol"],
                        lin["max_iter"], lin["class_weight_mode"],
                        cfg_tree["evaluate"]["corr_threshold"])
    hp = {"chosen_C": r.chosen_C, "curve": r.curve, "edge": r.edge,
          "df_effective": r.df_effective, "df_at_grid_ends": r.df_at_grid_ends,
          "n_iter": r.n_iter, "class_weight": str(r.class_weight),
          "weight_measurement": r.weight_measurement,
          "dup_pairs": r.pipeline.named_steps["dup"].report_}
    return r.pipeline, hp, r.inner_oof_scores, r.inner_oof_y


def _predict_ml_arm(arm: str, model, X: np.ndarray, names: list[str]) -> np.ndarray:
    """Uncalibrated P(y=1) per row for either branch.

    The linear branch goes through the pipeline, so the missing-value encoding
    and the scaling fitted on the training rows are reapplied to X here; the
    boosting branch consumes the NaN directly.
    """
    if arm == "trees":
        return trees_predict(model, X, names)
    return model.predict_proba(X)[:, 1]


def _fold_worker(payload: dict) -> dict:
    """One (τ, outer fold): fit the arms, calibrate, predict on the test set.

    Runs in a worker process: the config is reloaded from its path instead of being
    shipped in the payload, and the thread pinning is applied again inside the
    child, so every fit made here runs at the fixed thread count that keeps
    floating-point summations reproducible (see runtime.pin_threads).

    Returns the rows destined for oof.parquet, the fold artifact, and — only for
    the head arm at the primary horizon — the per-frame contributions.
    """
    cfg = load_config(payload["cfg_path"])
    pin_threads(cfg)
    tau = payload["tau"]
    fold_id = payload["fold_id"]
    data = pickle.loads(payload["blob"])
    X, names = data["X"], data["names"]
    y, valid = data["y"], data["valid"]
    assemblies, traj = data["assemblies"], data["traj"]
    meta = data["meta"]
    published = data["published"]
    tr_idx, te_idx = data["tr_idx"], data["te_idx"]
    check_no_frame_leak(tr_idx, te_idx, traj)
    seed = cfg["model.seed"]
    cal_cfg = (cfg["model.calibration.family"],
               cfg["model.calibration.min_positives_isotonic"])
    # Training uses only rows with a decided label: a frame whose horizon runs off
    # the end of the trajectory is unknown, not negative (see labels). The test
    # side keeps every row and carries `valid` into the OOF table, so the metric
    # layer applies the same censoring rule to every arm at once.
    tr_fit = tr_idx[valid[tr_idx]]
    oof_rows = []
    fold_art: dict = {"fold": fold_id, "tau": tau, "arms": {}}
    contrib_out = None

    # Which columns each arm is allowed to see. The two comparison arms are the
    # linear branch on a restricted subset — same fitting code, same folds, same
    # calibration — so their distance from the full linear arm is attributable to
    # the columns they were denied and to nothing else.
    arm_cols = {
        "trees": names, "linear": names,
        "linear_control_structural": control_columns(names),
        "clock": clock_columns(names),
    }
    # The two comparison arms always run; the config only chooses among the full
    # ones. A result without them would have nothing to be compared against.
    which = list(cfg["model.which"]) + ["linear_control_structural", "clock"]
    for arm in which:
        cols = arm_cols[arm]
        sub = [names.index(c) for c in cols]
        if len(sub) == 0 or len(np.unique(y[tr_fit])) < 2:
            fold_art["arms"][arm] = {"skipped": "no columns or a single class in training"}
            continue
        kind = "trees" if arm == "trees" else "linear"
        try:
            model, hp, oof_s, oof_y = _fit_ml_arm(
                kind, cfg.tree, X[np.ix_(tr_fit, sub)], y[tr_fit],
                assemblies[tr_fit], traj[tr_fit], cols,
                # One recorded seed; the per-fit seed is derived from it and from
                # the names of arm, horizon and fold, so two fits never share a
                # random stream and the whole run still reproduces from one number.
                derive_seed(seed, arm, str(tau), fold_id),
                times=np.asarray(meta["time_ps"])[tr_fit],
                # Two-sided embargo for the time-block fallback: the feature
                # window reaches backward and the label horizon forward, so a
                # block must be separated from its training rows by the sum.
                block_embargo_ps=cfg["features.window_ps"] + tau)
        except (RuntimeError, ValueError) as e:
            # A degenerate fold is handled explicitly: the arm is skipped BY
            # NAME and stays visible in the fold artifact.
            fold_art["arms"][arm] = {"skipped": f"degenerate training: {e}"}
            continue
        # Calibrator and threshold come from the INNER out-of-fold scores of this
        # fold; the test rows have not been touched at this point.
        cal = Calibrator(*cal_cfg).fit(oof_s, oof_y)
        thr = f1_threshold(cal.predict(oof_s), oof_y)
        raw_te = _predict_ml_arm(kind, model, X[np.ix_(te_idx, sub)], cols)
        prob_te = cal.predict(raw_te)
        fold_art["arms"][arm] = {"hyperparams": hp, "calibration": cal.describe(),
                                 "threshold": thr, "branch": "native" if kind == "trees"
                                 else "encoded", "n_features_seen": len(cols)}
        for k, i in enumerate(te_idx):
            oof_rows.append((meta["condition"][i], meta["replica"][i],
                             meta["time_ps"][i], tau, arm, float(raw_te[k]),
                             float(prob_te[k]), int(y[i]), bool(valid[i]),
                             bool(meta["in_transit"][i]), fold_id))
        # Mechanism diagnosis runs on the declared head arm at the primary horizon
        # only. Attributions of an arm that was not declared in advance would be a
        # second, unrecorded choice on top of the first.
        if arm == cfg["model.head"] and tau == cfg["labels.primary_tau_ps"]:
            if kind == "trees":
                contrib = trees_contributions(model, X[np.ix_(te_idx, sub)], cols)
            else:
                contrib = linear_contributions(model, X[np.ix_(te_idx, sub)])
            contrib_out = {"columns": cols, "rows_te": te_idx.tolist(),
                           "contrib": contrib.tolist(),
                           "attribution_kind": "TreeSHAP (xgboost pred_contribs, exact)"
                           if kind == "trees" else
                           "coef·z of the linear model in the post-preprocessing space"}
            # Block ablation: one block at a time, fold hyperparameters, no
            # re-search.
            # Re-selecting hyperparameters after dropping a block would mix two
            # effects — the loss of the block and a new selection — and the
            # difference could no longer be read as the block's contribution.
            # Only the raw scores are stored; evaluate turns them into per-fold
            # average precision, so the aggregation rule lives in one place.
            fold_art["block_ablation"] = {}
            for block, bcols in data["block_cols"].items():
                keep = [c for c in cols if c not in bcols]
                # Nothing dropped (the block is not in this arm's columns) or
                # nothing left: neither case is an ablation of this block.
                if len(keep) == len(cols) or not keep:
                    continue
                ksub = [names.index(c) for c in keep]
                m2, _hp2, s2, y2 = _fit_ml_arm(kind, cfg.tree, X[np.ix_(tr_fit, ksub)],
                                               y[tr_fit], assemblies[tr_fit],
                                               traj[tr_fit], keep,
                                               derive_seed(seed, arm, str(tau), fold_id, block),
                                               times=np.asarray(meta["time_ps"])[tr_fit],
                                               block_embargo_ps=cfg["features.window_ps"] + tau)
                raw2 = _predict_ml_arm(kind, m2, X[np.ix_(te_idx, ksub)], keep)
                fold_art["block_ablation"][block] = {
                    "scores": raw2.tolist(),
                }
    # PLS-FMA benchmark arm (Krivobokova et al. 2012): readiness regressed on
    # superimposed Cartesian C-alpha coordinates from the coords artifact.
    coords_art = _aligned_coords(cfg, meta)
    C, coords_dropped = coords_art if coords_art is not None else (None, 0)
    if C is not None and len(np.unique(y[tr_fit])) >= 2:
        try:
            r = fit_plsfma_fold(C[tr_fit], y[tr_fit], assemblies[tr_fit],
                                cfg["model.plsfma.components_grid"])
            # kind="raw": the PLS score is an unbounded regression value, not a
            # probability, so the sigmoid is fitted on the score itself.
            cal = Calibrator(cal_cfg[0], cal_cfg[1], kind="raw").fit(
                r.inner_oof_scores, r.inner_oof_y)
            thr = f1_threshold(cal.predict(r.inner_oof_scores), r.inner_oof_y)
            raw_te = plsfma_predict(r.model, C[te_idx])
            prob_te = cal.predict(raw_te)
            fold_art["arms"]["plsfma_coords"] = {
                "hyperparams": {"chosen_components": r.chosen_components,
                                "curve": r.curve, "edge": r.edge,
                                "mode_vector_norm": float(np.linalg.norm(r.mode_vector)),
                                "ewmcm_vector_norm": float(np.linalg.norm(r.ewmcm_vector))},
                "calibration": cal.describe(), "threshold": thr,
                "branch": "cartesian_coords", "n_features_seen": int(C.shape[1]),
                "n_union_columns_dropped": coords_dropped}
            for k, i in enumerate(te_idx):
                oof_rows.append((meta["condition"][i], meta["replica"][i],
                                 meta["time_ps"][i], tau, "plsfma_coords",
                                 float(raw_te[k]), float(prob_te[k]), int(y[i]),
                                 bool(valid[i]), bool(meta["in_transit"][i]), fold_id))
        except (RuntimeError, ValueError) as e:
            fold_art["arms"]["plsfma_coords"] = {"skipped": str(e)}
    else:
        fold_art["arms"]["plsfma_coords"] = {
            "skipped": "coords artifact missing or single-class training fold"}

    # Published criterion: score comes from the table, no training; no
    # probability scale.
    # It appears on the test rows of every fold so that it is scored on exactly the
    # frames the trained arms are scored on, even though it has no training side.
    for arm, colname in zip(PUBLISHED_ARMS, ("bl_rao_E_kJmol", "bl_rao_E_win_kJmol")):
        # Sign: the column is a dewetting barrier in kJ/mol, so a HIGH energy means
        # a dry, non-conducting constriction. Negating puts the arm on the same
        # orientation as the trained ones — larger score, readier frame — which is
        # what the ranking metrics assume. prob_cal below is NaN by construction:
        # an energy has no probability scale, and none is invented for it.
        s = -published[colname]
        for i in te_idx:
            oof_rows.append((meta["condition"][i], meta["replica"][i],
                             meta["time_ps"][i], tau, arm, float(s[i]),
                             float("nan"), int(y[i]), bool(valid[i]),
                             bool(meta["in_transit"][i]), fold_id))
    return {"oof_rows": oof_rows, "fold_art": fold_art, "contrib": contrib_out,
            "te_idx": te_idx.tolist()}


# Column order of oof.parquet. One row per (frame, horizon, arm): the frame keys
# identify the frame, score_raw is the arm's own scale and prob_cal the calibrated
# probability where one exists, y and valid carry the label and its censoring, and
# in_transit marks a frame with a permeant already inside, on which the descriptors
# can see the crossing ion itself — evaluate scores resting and transit frames apart.
OOF_COLUMNS = ["condition", "replica", "time_ps", "tau", "arm", "score_raw",
               "prob_cal", "y", "valid", "in_transit", "fold"]


def run_step(cfg: Config) -> None:
    """The train step: nested CV over all horizons, then the final model and reports.

    Writes oof.parquet, folds.json, contrib_head.parquet, attribution.json,
    condition_transfer.json and the final model to the step's output directory,
    then hands the same in-memory results to the evaluate and interpret steps so
    that reporting reads exactly the numbers that were just produced.
    """
    pin_threads(cfg)
    table, schema, verdict = load_features(cfg)
    labels = load_labels(cfg)
    # one_to_one: features and labels are both one row per frame, and a merge that
    # duplicated or dropped rows would change the sample without saying so. The
    # assertion then confirms every feature row found its label.
    merged = table.merge(labels, on=["condition", "replica", "time_ps"],
                         validate="one_to_one")
    assert len(merged) == len(table), "labels and features are on different frame grids"
    with step_output(cfg, "train") as out:
        log = StepLog(out)
        mat = build_matrix(cfg, merged, schema, verdict, log=log.say)
        names = mat["columns"]
        X = mat["X"].to_numpy(float)
        assemblies, traj = mat["assembly"], mat["trajectory"]
        # Block membership of the columns that survived into the matrix; the
        # ablation drops one block at a time and needs the names, not the schema.
        block_cols = {}
        for c in schema:
            if c.name in names:
                block_cols.setdefault(c.block, []).append(c.name)
        folds = outer_folds(assemblies)
        log.say(f"outer folds (over independent assemblies): {len(folds)}; "
                f"the effective sample size is the number of trajectories")
        taus = cfg["labels.tau_ps"]
        jobs = []
        for tau in taus:
            y = merged[f"y_tau{int(tau)}"].to_numpy(int)
            valid = merged[f"valid_tau{int(tau)}"].to_numpy(bool)
            meta = {"condition": merged["condition"].to_numpy(),
                    "replica": merged["replica"].to_numpy(),
                    "time_ps": merged["time_ps"].to_numpy(),
                    "in_transit": merged["in_transit"].to_numpy()}
            published = {c: merged[c].to_numpy(float)
                         for c in ("bl_rao_E_kJmol", "bl_rao_E_win_kJmol")}
            # Every (horizon, fold) pair is an independent job: the folds share no
            # fitted object, so they are dispatched to worker processes. The matrix
            # is serialized once per job rather than shared, which costs memory and
            # keeps each worker's data immutable.
            for fi, (tr_idx, te_idx, gname) in enumerate(folds):
                blob = pickle.dumps({"X": X, "names": names, "y": y, "valid": valid,
                                     "assemblies": assemblies, "traj": traj,
                                     "meta": meta, "published": published,
                                     "tr_idx": tr_idx, "te_idx": te_idx,
                                     "block_cols": block_cols})
                jobs.append({"cfg_path": str(cfg.source_path), "tau": tau,
                             "fold_id": f"{fi}:{gname}", "blob": blob})
        results = []
        with ProcessPoolExecutor(max_workers=min(len(jobs), 14)) as ex:
            # ex.map preserves the order of jobs, so the artifacts are written in a
            # deterministic order however the workers finish.
            for res in ex.map(_fold_worker_entry, jobs):
                results.append(res)
                log.say(f"fold {res['fold_art']['fold']} τ={res['fold_art']['tau']} done")
        oof = pd.DataFrame([r for res in results for r in res["oof_rows"]],
                           columns=OOF_COLUMNS)
        oof.to_parquet(out / "oof.parquet")
        fold_arts = [res["fold_art"] for res in results]
        (out / "folds.json").write_text(json.dumps(fold_arts, indent=1,
                                                   ensure_ascii=False, default=str))
        # Head-arm contributions at the primary τ (for mechanism diagnosis).
        contribs = [res["contrib"] for res in results if res["contrib"]]
        contrib_df = _assemble_contrib(contribs, merged)
        contrib_df.to_parquet(out / "contrib_head.parquet")
        (out / "attribution.json").write_text(json.dumps(
            {"kind": contribs[0]["attribution_kind"] if contribs else "none",
             "head_arm": cfg["model.head"]}, ensure_ascii=False, indent=1))

        # Block ablation: per-fold AP aggregation is computed by evaluate.
        # Transfer between conditions.
        cond_transfer = _condition_transfer(cfg, X, names, merged, assemblies, traj, log)
        (out / "condition_transfer.json").write_text(
            json.dumps(cond_transfer, indent=1, ensure_ascii=False, default=str))

        # Final deliverable model: fit on all valid rows; parameters follow the
        # recorded selection rule; the calibrator is fit on the pooled OOF scores.
        final = _final_model(cfg, X, names, merged, assemblies, traj, oof, fold_arts)
        (out / "final_model.pkl").write_bytes(pickle.dumps(final))
        # The same record without the fitted object: readable, diffable, and
        # independent of the library version that could unpickle the model.
        (out / "final_model.json").write_text(json.dumps(
            {k: v for k, v in final.items() if k not in ("model",)},
            indent=1, ensure_ascii=False, default=str))

        evaluate_mod.evaluate_run(cfg, oof, merged, fold_arts, results, out, log)
        interpret_mod.interpret_run(cfg, oof, contrib_df, schema, names, merged, out, log)
        log.close()


def _fold_worker_entry(payload):
    """Target handed to the process pool; the work itself is in _fold_worker."""
    return _fold_worker(payload)


def _aligned_coords(cfg, meta) -> tuple[np.ndarray, int] | None:
    """Superimposed C-alpha coordinates aligned to the model-matrix row order.

    Conditions with differently numbered topologies (the KcsA mutants) yield a
    union column set in the coords artifact: a column absent from one condition
    arrives as missing there. Those columns are a representation artifact, not
    a measurement gap, and are dropped — the analysis runs on the C-alphas all
    conditions share, and the number dropped is reported alongside. A residual
    non-finite value after that IS a data gap and refuses as before.
    Returns None when the coords artifact has not been produced or is unusable.
    """
    path = run_dir(cfg) / "coords" / "coords.parquet"
    if not path.exists():
        return None
    cdf = pd.read_parquet(path)
    # Left merge onto the frame keys of the model matrix, in that order: the row
    # index of the returned matrix has to match X row for row, since the caller
    # indexes both with the same fold indices.
    key = pd.DataFrame({"condition": meta["condition"], "replica": meta["replica"],
                        "time_ps": meta["time_ps"]})
    merged = key.merge(cdf, on=["condition", "replica", "time_ps"], how="left",
                       validate="one_to_one")
    # float32 is the artifact's own precision and halves a matrix of 3 numbers per
    # C-alpha per frame; PLS runs in float64 internally (see pls_denham).
    mat = merged.drop(columns=["condition", "replica", "time_ps"]).to_numpy(np.float32)
    # Column-wise: a coordinate column is kept only if it is measured in EVERY
    # frame. This is what removes the union columns of the differently numbered
    # topologies, since those are missing across whole conditions at a time.
    keep = np.isfinite(mat).all(axis=0)
    n_dropped = int((~keep).sum())
    mat = mat[:, keep]
    if mat.shape[1] == 0 or np.any(~np.isfinite(mat)):
        return None
    return mat, n_dropped


def _assemble_contrib(contribs: list[dict], merged: pd.DataFrame) -> pd.DataFrame:
    """One table of per-frame contributions from the per-fold pieces.

    Every fold contributed its own test rows, so the pieces are disjoint and their
    concatenation covers each frame once. Returns an empty frame when no fold
    produced contributions, which happens when the head arm was skipped in every
    fold. The frame keys are recovered from the stored row indices into `merged`,
    the same frame order the matrix was built in.
    """
    rows, datas = [], []
    for c in contribs:
        rows.extend(c["rows_te"])
        datas.append(np.asarray(c["contrib"]))
    if not contribs:
        return pd.DataFrame()
    cols = contribs[0]["columns"]
    data = np.vstack(datas)
    df = pd.DataFrame(data, columns=cols)
    idx = np.asarray(rows)
    # insert(0, ...) in this order leaves the keys as condition, replica, time_ps
    # ahead of the contribution columns.
    for m in ("condition", "replica", "time_ps"):
        df.insert(0, m, merged[m].to_numpy()[idx])
    return df


def _condition_transfer(cfg, X, names, merged, assemblies, traj, log) -> dict:
    """Train on one condition, test on another: was the condition label learned?

    Returns a dict keyed "train->test" holding raw scores, labels and the validity
    mask of the test condition, or a named reason for skipping. Only the head arm
    at the primary horizon is used, and only raw scores are stored: the metrics are
    computed downstream by evaluate, so this control is scored by the same code as
    everything else. A control that transfers no better than chance means the arm
    had learned to recognize the condition rather than the conformation.
    """
    conds = [c["id"] for c in cfg["data.conditions"]]
    if len(conds) < 2:
        return {"skipped": "single condition"}
    tau = cfg["labels.primary_tau_ps"]
    y = merged[f"y_tau{int(tau)}"].to_numpy(int)
    valid = merged[f"valid_tau{int(tau)}"].to_numpy(bool)
    cond_arr = merged["condition"].to_numpy()
    out = {}
    head = cfg["model.head"]
    for train_c in conds:
        for test_c in conds:
            if train_c == test_c:
                continue
            # Training on the decided rows of one condition, testing on all rows of
            # another. No embargo is needed between the two: they are different
            # simulations, so no feature window or label horizon spans them.
            tr = np.flatnonzero((cond_arr == train_c) & valid)
            te = np.flatnonzero(cond_arr == test_c)
            if len(np.unique(y[tr])) < 2:
                out[f"{train_c}->{test_c}"] = {"skipped": "single class"}
                continue
            try:
                model, hp, oof_s, oof_y = _fit_ml_arm(
                    "trees" if head == "trees" else "linear", cfg.tree,
                    X[tr], y[tr], assemblies[tr], traj[tr], names,
                    # Offset from the recorded seed so this control never reuses
                    # the random stream of a cross-validation fold.
                    cfg["model.seed"] + 777)
            except (RuntimeError, ValueError) as e:
                # A single-replica or event-starved condition cannot sustain the
                # rotation-based fit; the control is skipped by name.
                out[f"{train_c}->{test_c}"] = {"skipped": f"degenerate training: {e}"}
                continue
            raw = _predict_ml_arm("trees" if head == "trees" else "linear",
                                  model, X[te], names)
            out[f"{train_c}->{test_c}"] = {
                "scores": raw.tolist(), "y": y[te].tolist(),
                "valid": valid[te].tolist(), "n_train": int(len(tr)),
            }
            log.say(f"condition transfer {train_c}->{test_c}: trained on {len(tr)} rows")
    return out


def _final_model(cfg, X, names, merged, assemblies, traj, oof, fold_arts) -> dict:
    """The deliverable model: head arm at the primary horizon, fitted on all valid rows.

    Returns the record written to final_model.pkl (and, minus the fitted object, to
    final_model.json): the model, the column order it expects, its hyperparameters,
    its calibrator and its threshold, plus the code fingerprint. It carries no
    performance number of its own — this fit has no held-out data left, and the
    reported metrics are the cross-validated ones.
    """
    tau = cfg["labels.primary_tau_ps"]
    head = cfg["model.head"]
    y = merged[f"y_tau{int(tau)}"].to_numpy(int)
    valid = merged[f"valid_tau{int(tau)}"].to_numpy(bool)
    rows = np.flatnonzero(valid)
    kind = "trees" if head == "trees" else "linear"
    try:
        model, hp, oof_s, oof_y = _fit_ml_arm(kind, cfg.tree, X[rows], y[rows],
                                              assemblies[rows], traj[rows], names,
                                              # Again offset, for the same reason as
                                              # in the transfer control above.
                                              cfg["model.seed"] + 999,
                                              times=merged["time_ps"].to_numpy()[rows],
                                              block_embargo_ps=cfg["features.window_ps"] + tau)
    except (RuntimeError, ValueError) as e:
        # An event-starved system cannot ship a model; the absence is recorded,
        # and applying this system as a transfer source will be refused.
        return {"model": None, "arm": head, "tau_ps": tau,
                "unavailable": f"degenerate training on the full sample: {e}",
                "columns_in_order": names, "code": code_fingerprint()["tree_sha256"]}
    # The calibrator of the shipped model is fitted on the POOLED cross-validated
    # scores of the head arm — every frame scored by a model that had not seen its
    # assembly. Fitting it on the final model's own predictions would calibrate
    # against in-sample scores and promise a sharpness the model does not have.
    pool = oof[(oof["arm"] == head) & (oof["tau"] == tau) & oof["valid"]]
    cal = Calibrator(cfg["model.calibration.family"],
                     cfg["model.calibration.min_positives_isotonic"]).fit(
        pool["score_raw"].to_numpy(), pool["y"].to_numpy())
    thr = f1_threshold(cal.predict(pool["score_raw"].to_numpy()), pool["y"].to_numpy())
    return {"model": model, "arm": head, "tau_ps": tau,
            "columns_in_order": names,
            "hyperparams": hp, "calibration": cal.describe(),
            "calibrator": cal, "threshold": thr,
            "hyperparam_rule": "refit on all valid rows with parameters re-selected "
                               "by the nested procedure on the full sample; final-model "
                               "metrics are taken from the cross-validation",
            "code": code_fingerprint()["tree_sha256"],
            "note": "the model is unusable without the schema: columns are carried by name"}
