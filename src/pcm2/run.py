"""Entry point: python -m pcm2.run <step> --config configs/<system>.yaml.

Pipeline steps: each one recomputes from scratch and replaces its own output
directory wholesale.

ORDER is the dependency order, and it is what makes running a single step
meaningful: a step reads the finished directories of the steps before it under
runs/<system>-stride<N>/ and never their intermediate state, because
runtime.step_output publishes a directory only after the step has completed. So
a lone step assumes its predecessors have already run on the same config, and
`all` walks the whole chain. Since nothing is cached, a rerun is a full
recomputation; there is no partial state to reconcile and no stale file to
inherit.
"""

from __future__ import annotations

import argparse
import sys

from . import autodetect, benchmark, config, coords, events, report, viz
from . import features as features_mod
from . import labels as labels_mod
from .runtime import pin_threads

# CLI verb -> the function that writes that step's directory. "train" is
# deliberately missing: it is dispatched through _train so that this module can
# be imported without the model layer.
STEP_FUNCS = {
    "autodetect": autodetect.run_step,
    "events": events.run_step,
    "features": features_mod.run_step,
    "coords": coords.run_step,
    "labels": labels_mod.run_step,
    "figures": viz.run_step,
    "report": report.run_step,
    "benchmark": benchmark.run_step,
}

# What each step needs from the ones before it: autodetect writes the
# selections, the pore axis and the axial offsets into the config, and events,
# features and coords all read the trajectories through those selections;
# labels needs the crossing times from events; train needs the feature table
# and the labels on the same frame grid; figures, report and benchmark only read
# what train wrote. coords is the one soft dependency — its Cartesian C-alpha
# table feeds the PLS-FMA arm, and train skips that arm if the table is absent.
ORDER = ["autodetect", "events", "features", "coords", "labels", "train",
         "figures", "report", "benchmark"]


def _train(cfg):
    """Run the training step."""
    # Imported here rather than at module level: the model layer pulls in
    # scikit-learn and xgboost, and a run of the reading-layer steps should not
    # have to load them.
    from . import models
    models.run_step(cfg)


def main(argv=None) -> int:
    """Run the requested step (or the whole chain) for one system; returns the exit code."""
    p = argparse.ArgumentParser(prog="pcm2")
    p.add_argument("step", choices=ORDER + ["all", "transfer"])
    p.add_argument("--config", required=True)
    p.add_argument("--source-config", help="for the transfer step: source-system config")
    args = p.parse_args(argv)
    cfg = config.load(args.config)
    # Threads are pinned before any step imports a linear-algebra backend: the
    # OMP/BLAS variables are read when a library builds its thread pool, so a
    # later assignment would have no effect and the summation order — hence the
    # last digits of every score — would stop being reproducible.
    pin_threads(cfg)
    if args.step == "transfer":
        # transfer is not a pipeline step: it consumes two finished per-system
        # runs and writes one shared artifact of its own under runs/, so it is
        # neither in ORDER nor in STEP_FUNCS. Both directions are computed in the
        # one call and land in the one file, so the pair is always reported
        # together rather than one direction at a time.
        from . import transfer
        src = config.load(args.source_config)
        out = transfer.run_both_directions(src, cfg)
        print(out)
        return 0
    steps = ORDER if args.step == "all" else [args.step]
    for s in steps:
        print(f"=== step {s} ({cfg['system.id']}) ===", flush=True)
        if s == "train":
            _train(cfg)
        else:
            STEP_FUNCS[s](cfg)
        if s == "autodetect":
            # The autodetect step saves its detected keys back into the YAML
            # file in place; the cfg object held here was loaded before that
            # write, so every later step in an `all` run would otherwise see the
            # config as it was before detection.
            cfg = config.load(args.config)  # autodetect appended to the config — reload
    return 0


if __name__ == "__main__":
    sys.exit(main())
