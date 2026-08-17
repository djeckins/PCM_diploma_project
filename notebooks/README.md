# Interactive notebooks

Three notebooks and the library they share. The pipeline in `src/pcm2` measures
and trains; the notebooks ask questions of what it wrote — a verdict for one
frame of one trajectory, the per-frame output of every surveyed method at one
moment, or the state of a run that has just finished.

None of them recomputes a pipeline step. They read artefacts under `../runs/`
and print or plot; the only thing ever written back is a pooled-model cache in
`../runs/pooled-predict/`. Kernel: the `chem` environment of
`../environment.yml`. Either working directory works, because `predict_lib.py`
changes to the project root on import — the system configs address their files
relative to it.

## Which one to open

| notebook | the question it answers | open it when |
|---|---|---|
| `run_my_trajectory.ipynb` | did my run produce what it should, and what did autodetect decide or refuse? | straight after `python -m pcm2.run all --config configs/<system>.yaml` on a system of your own |
| `permeation_predict.ipynb` | is this channel about to conduct — at this frame, at this horizon — and what carries the verdict? | when applying the frozen monitor to a trajectory, one of the three training proteins or a new one |
| `methods_showcase.ipynb` | what does each surveyed method report on the same frames, in its own units? | when placing this work beside the published criteria, or rebuilding the one-moment table of the write-up |

All three run on the artefacts shipped in this repository, with no trajectory
data present. Only `permeation_predict.ipynb` ever opens a trajectory, and only
when the query falls outside the feature table of the system it is pointed at.

## `run_my_trajectory.ipynb` — read a finished run

The entry notebook: point it at a system and read what the run produced. No
value is typed in — no temperature, no atom selection, no number of sites. The
config path is the only editable field, and even that is overridden by the
`PCM2_CONFIG` environment variable, so the notebook can be executed
non-interactively without touching a cell. It recomputes nothing by design: the
features step *writes* its output, so a cell that looked harmless would
overwrite a published run.

The cells print, in order: the origin and basis recorded for every
autodetected key, and every key autodetect refused with its reason; per-replica
event counts (detected here against provided) with the left-censored count and
the label summary per horizon; the first rows of `train/answers.parquet` with
the READY fraction and the mechanism-axis counts; per-arm average precision at
the primary horizon with its ratio to chance; and the precision-recall figure of
that run. A missing artefact prints its own name and the command that would
produce it instead of raising, so a partially finished run is still readable —
on Cx43, where every fold was refused, the answers cell prints that line rather
than a table. The last cell is empty on purpose: it is where a manual override
goes, and it stays empty unless autodetect refused a key.

Requires one run directory and nothing else. On the shipped `gramicidin-stride4`
run it executes end to end in about a second.

## `permeation_predict.ipynb` — apply the monitor

The pooled monitor, trained once on the conducting systems and then frozen, is
applied to a trajectory and answers three questions on the same frame set:
what is measured at the pore at time t (no model, no horizon); whether a
crossing completes in (t, t + τ] together with the mechanism behind that verdict
and the same frame at the other cached horizons; and how readiness runs over a
stretch, against the crossings that actually happened.

Edit the first code cell only:

```python
PROTEIN = "gramicidin"    # system id from configs/ — or a path to the data folder
CONDITION = None          # None -> first condition of the system
REPLICA = None            # None -> first replica of the condition

TAU_PS = 1000.0           # prediction horizon: "a crossing completes within τ ps"
VARIANT = "physics_only"  # "physics_only" (recommended off-family) or "full"

TIME_PS = 100000.0              # the single moment for questions 1 and 2
RANGE_PS = (50000.0, 150000.0)  # the stretch for question 3
HORIZONS_PS = (500.0, 1000.0, 2000.0, 4000.0)  # horizon sweep shown in question 2
```

A system id carries no machine paths, so it is the portable way to name a
protein; a folder path also works and is matched against the `data_root`
recorded in the configs.

Run all cells. The point query prints one verdict with its ranked mechanism
axes and the descriptor driving each:

```
t = 100000 ps | horizon τ = 1000 ps | P(crossing within τ) = 0.022 → NOT READY
  1. constriction dehydration — weight 57% (driver: flu_r_constr_var_win = 0.0908)
  2. pathway/lining rearrangement — weight 19% (driver: geo_bin7_nsearch = 2)
  3. permeant delivery — weight 17% (driver: dlv_t_since_cross_ps = 2.12e+04)
```

The range query prints the READY fraction over the stretch and draws the
probability timeline with the calibrated threshold and the observed crossings
marked. Nothing is saved: the figure is inline, and the only file the notebook
can create is a pooled cache for a horizon that had none.

Requires a config, one `events` run, and either the system's feature table or
the trajectory itself; the five configured systems all ship their tables, so
they can be queried as they are. At a cached (τ, variant) the whole notebook
runs in about three seconds on the shipped gramicidin tables. A horizon that is
not cached costs one pooled fit before the first verdict — roughly a minute,
once per (τ, variant).

## `methods_showcase.ipynb` — every method on the same frames

The comparison notebook: each benchmarked method's own per-frame output on the
same frames, in its native units, with its verdict against its own declared
threshold and its measured discrimination there. Every number is read
or computed from the run artefacts when the cell executes; nothing is
transcribed. Change `PROTEIN` and `RANGE_NS` in the first cell to move the
comparison to another system.

It produces four things: one trace panel per method that emitted output over the
chosen stretch, with completed crossings marked; the full readout at one
frame, including this work's calibrated probability, verdict and mechanism axes;
a capability matrix of which output forms each method emits over the trajectory;
and the headline metrics of the run's benchmark step. `figH_methods_showcase.png`
and the one-moment table in `../runs/FINAL/` were exported from here by hand —
`../docs/REPRODUCTION.md` records them as the two write-up artefacts that no
script regenerates.

Requires a `benchmark` step next to the `features` and `events` of the same run;
all five systems ship theirs. About three seconds end to end on the shipped
gramicidin artefacts.

## Helper modules

| file | purpose |
|---|---|
| `predict_lib.py` | everything the notebooks call, kept out of the cells so it can be tested: config resolution from an id or a folder, pooled-model training and caching, feature rows for a query, prediction, per-verdict explanation, and the one-moment state readout |
| `_train_pooled.py` | the pooled fit as a command (`python _train_pooled.py <tau_ps> <variant>`), kept separate so that importing the notebook library never triggers a fit; use it to warm the cache for a horizon before opening a notebook |
| `verify_predict.py` | checks on the prediction path: nine horizons cross-checked against the pipeline's own labels, seven time-range edge cases, and window-mode rows compared against the pipeline tables. Run it after changing `predict_lib.py` |

The first two blocks of `verify_predict.py` read tables only and finish in a few
seconds; the third measures rows from a trajectory and therefore stops with
`IoError: missing topology …` when the trajectories are not present. That is the
expected outcome in a checkout of this repository alone.

## Implementation

- The model is one pooled classifier trained on every conducting trajectory
  unit — 15 units of gramicidin A, MthK and KcsA-E71A, 26 672 rows at τ = 1000
  ps, 12.2 % of them positive. The first request for a given (τ, variant) fits
  it and caches it in `../runs/pooled-predict/`; later requests read the cache.
  `physics_only` removes columns whose per-system value ranges do not overlap
  (79 of 100 model columns survive at τ = 1000 ps), so the model cannot lean on
  system identity.
- Horizons already carried by the pipeline tables load as they are; any other τ
  is relabelled from the recorded crossings by the same rule, and is
  bit-identical where both exist — `verify_predict.py` asserts that agreement
  on every row that carries a pipeline label.
- Held-out quality falls smoothly with the horizon, so short horizons are the
  informative ones. Nine horizons are cached in both variants; for
  `physics_only` the leave-one-unit-out average precision runs from ×5.9 over
  chance at 250 ps to ×1.7 at 6 ns, read from the `oof_quality` field of the
  cached bundles and printed by `verify_predict.py`.
- Feature rows come from the protein's pipeline feature table when that run
  exists. Otherwise the notebook measures only the frames inside the query plus
  one backward window: minutes, against hours for a full pass. Both paths
  produce identical model columns; the window-mode caveat is two run-length ages
  capped by the window length.
- A new protein needs its config (`../tools/init_configs.py`), one `autodetect`
  run and one `events` run (`python -m pcm2.run <step> --config
  configs/<system>.yaml`); after that it can be queried like any other, and the
  frozen pooled model is applied to it by feature name. Labels are never needed
  for the target system. The training pool stays fixed at the three conducting
  proteins, so a new protein is scored without being added to it.

## Limitations

- On a protein family absent from training the verdicts carry zero-shot
  transfer quality (measured in `../runs/monitor-generalisation-tau1000/`),
  not in-system quality.
- The explanation weights are TreeSHAP contributions in log-odds, grouped
  into the interpretation axes of `../src/pcm2/interpret.py`; they report what
  moved this verdict and carry no causal claim.
