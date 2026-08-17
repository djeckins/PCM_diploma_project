"""Applying the model to a channel it has never seen.

1. Columns are matched BY NAME, never by position.
2. Columns absent in the new protein arrive as missing values and stay missing.
3. A transfer plan is printed: which columns are absent and what share of the
   explanatory power they carried.
4. Refusal if more than the declared share is lost: the output of a model
   missing that much of its input is indistinguishable from a valid one.
"""

from __future__ import annotations

import json
import pickle
import shutil

import numpy as np
import pandas as pd

from .config import Config
from .evaluate import ap_with_base
from .models import _predict_ml_arm
from .runtime import PROJECT_ROOT, run_dir, write_provenance


def run_transfer(source_cfg: Config, target_cfg: Config) -> dict:
    """Apply the source system's final model to the target system's frames, or refuse.

    Returns the transfer plan: which of the model's columns the target protein
    provides, what share of the source model's explanatory power the absent ones
    carried, the verdict, and — where the target has labels — the average precision
    of the ranking on the target. Nothing is fitted on the target: the model, its
    calibrator and the horizon all come from the source run.
    """
    src_train = run_dir(source_cfg) / "train"
    final = pickle.loads((src_train / "final_model.pkl").read_bytes())
    # An event-starved source ships a record with model=None instead of a model. Its
    # own reason is quoted rather than replaced, so the refusal stays traceable to
    # the refit that failed on the source's full sample.
    if final.get("model") is None:
        return {"source": source_cfg["system.id"], "target": target_cfg["system.id"],
                "verdict": "REFUSED",
                "reason": "the source system shipped no model: "
                          + str(final.get("unavailable", "unknown"))}
    # Mean absolute contribution per column, as the interpretation stage recorded it
    # in the source's train artifact. It is the weight in which the loss below is
    # measured: what matters is not how many columns the target is missing, but how
    # much of the source model rested on them.
    imp = json.loads((src_train / "importances.json").read_text())["per_column"]
    cols = final["columns_in_order"]

    tgt_feats = pd.read_parquet(run_dir(target_cfg) / "features" / "features.parquet")
    n = len(tgt_feats)
    # Prefilled with NaN, then filled by name: a column the target does not have
    # stays missing and is left to the model's own missing-value policy. Absence is
    # normally structural — the target has no such site, no such subunit — so there
    # is nothing to impute, and a value borrowed from the target's other columns
    # would enter the model as a measurement that was never made.
    X = np.full((n, len(cols)), np.nan)
    missing, present = [], []
    for j, c in enumerate(cols):
        if c in tgt_feats.columns:
            X[:, j] = tgt_feats[c].to_numpy(float)
            present.append(c)
        else:
            missing.append(c)
    # The 1e-12 floor only guards the division: a model whose recorded importances
    # are all zero would otherwise raise here instead of reaching its verdict.
    total_imp = sum(abs(v) for v in imp.values()) or 1e-12
    lost = sum(abs(imp.get(c, 0.0)) for c in missing)
    lost_frac = lost / total_imp
    plan = {"source": source_cfg["system.id"], "target": target_cfg["system.id"],
            "n_columns_model": len(cols), "n_present": len(present),
            "missing_columns": missing,
            "importance_lost_frac": lost_frac,
            "max_importance_loss_frac": source_cfg["transfer.max_importance_loss_frac"]}
    # Past the declared share the answer is withheld. A model run on a matrix whose
    # informative half is missing still returns confident-looking probabilities, and
    # there is nothing in the output that would mark them as unsupported.
    if lost_frac > source_cfg["transfer.max_importance_loss_frac"]:
        plan["verdict"] = "REFUSED"
        plan["reason"] = (f"lost {lost_frac:.2f} of explanatory power > "
                          f"threshold {source_cfg['transfer.max_importance_loss_frac']}")
        return plan
    plan["verdict"] = "APPLIED"
    kind = "trees" if final["arm"] == "trees" else "linear"
    raw = _predict_ml_arm(kind, final["model"], X, cols)
    # The source's calibrator, unchanged: refitting it on the target would need the
    # target's labels and would turn the transfer into a fit on the system it claims
    # to be transferring to. The probabilities therefore live on the source's scale,
    # which is why only the ranking is scored below.
    prob = final["calibrator"].predict(raw)
    plan["scores_summary"] = {"prob_mean": float(np.nanmean(prob)),
                              "prob_p95": float(np.nanquantile(prob, 0.95))}
    # Check against target labels when they exist: ranking only, no fitting.
    lab_path = run_dir(target_cfg) / "labels" / "labels.parquet"
    if lab_path.exists():
        lab = pd.read_parquet(lab_path)
        # The horizon of the SOURCE model: it was trained to forecast that τ, and a
        # label column for another horizon would be a different question. The target
        # labels only carry it if the two configs sweep the same horizons.
        tau = source_cfg["labels.primary_tau_ps"]
        ycol, vcol = f"y_tau{int(tau)}", f"valid_tau{int(tau)}"
        if ycol in lab.columns:
            # The one-to-one merge keeps the row order of the feature table, which is
            # the order of raw, so labels and scores line up row by row. Had the
            # label table not covered every frame, the shorter mask would raise on
            # raw[v] rather than quietly score a shifted pairing.
            mg = tgt_feats[["condition", "replica", "time_ps"]].merge(
                lab, on=["condition", "replica", "time_ps"], validate="one_to_one")
            v = mg[vcol].to_numpy(bool)
            # Scored on the raw score, not the probability, and only on frames whose
            # horizon is complete: what transfers is the ordering of frames.
            plan["target_eval"] = ap_with_base(mg[ycol].to_numpy(int)[v], raw[v])
    return plan


def run_both_directions(cfg_a: Config, cfg_b: Config) -> dict:
    """Transfer a→b and b→a, written to one run directory as transfer.json.

    The two directions are separate experiments and neither implies the other: each
    uses a different training sample, a different column set and its own source
    config, down to the horizon and the refusal threshold, which both come from
    whichever system is the source. Returns the two plans keyed a_to_b and b_to_a.
    """
    out = {"a_to_b": run_transfer(cfg_a, cfg_b), "b_to_a": run_transfer(cfg_b, cfg_a)}
    # The pair and both strides are in the directory name: the same two systems read
    # at different strides are different results and must not overwrite each other.
    name = (f"transfer-{cfg_a['system.id']}-stride{cfg_a['data.stride']}"
            f"--{cfg_b['system.id']}-stride{cfg_b['data.stride']}")
    final = PROJECT_ROOT / "runs" / name
    # Write into a temporary directory next to the final one and rename it into
    # place: the same wholesale replacement the pipeline steps get from step_output.
    # Transfer is not one of the pipeline steps, so that helper refuses its name and
    # the sequence is spelled out here. An interrupted run leaves the previous result
    # untouched, with the temporary directory left behind as the sign of the failure.
    tmp = PROJECT_ROOT / "runs" / f".{name}.tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    (tmp / "transfer.json").write_text(json.dumps(out, indent=1, ensure_ascii=False,
                                                  default=str))
    # Provenance is recorded against cfg_a with the counterpart named, so the pair
    # can be identified from the artifact without reading the directory name.
    write_provenance(tmp, cfg_a, "transfer", extra={"counterpart": cfg_b["system.id"]})
    if final.exists():
        shutil.rmtree(final)
    tmp.rename(final)
    return out
