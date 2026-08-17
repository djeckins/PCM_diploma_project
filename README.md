# PCM 2 — permeation-readiness forecasting for ion channels

A per-frame monitor for ion-channel MD trajectories. For every frame it answers
two questions: whether the configuration is ready to conduct an ion within the
next τ ps, and if not, which mechanism is holding it back.

This repository is the code, the configurations and the run artefacts behind an
MSc dissertation (module CHE701P, *Artificial Intelligence for Drug Discovery*,
Queen Mary University of London, 2026). It is the first public release of the
work.

The trajectories themselves are not in it, and most of it works without them:
the feature and label tables ship, so the notebooks, the reports, every figure
and the whole test suite run on a clean clone.

## What to install

Python 3.11. Nine libraries do the work, and each earns its place:

| library | floor | role; what breaks without it |
|---|---|---|
| MDAnalysis | ≥ 2.8 | reading layer: streaming xtc/tpr input, PBC distances in triclinic cells. Without it neither autodetection nor features exist. The floor is 2.8 because that release relicensed the package to the LGPL |
| numpy, scipy | — | all numerical work; scipy supplies the density peaks that become sites. No optimiser is used |
| numba | — | inscribed-sphere kernel of the pore profile; without it the profile is prohibitively slow |
| pandas, pyarrow | — | feature, label and out-of-fold score tables on disk (parquet) |
| scikit-learn | — | linear arm (ridge logistic regression, newton-cholesky), calibration, metrics |
| xgboost | — | boosting with native missing-value handling, monotone constraints and a readable base_score |
| matplotlib | — | the figure catalogue |
| PyYAML | — | configs with a recorded basis per key |
| pytest, ruff | — | the test suite and the lint configuration (`[dev]` extra) |

Not dependencies, but worth knowing about:

- Jupyter is not declared anywhere here, because the pipeline does not import
  it. Install it separately if you want the notebooks (step 4 below).
- GROMACS is not required to run any code in this repository — MDAnalysis
  parses the `.tpr` itself. Its only role was one reading step outside the
  pipeline: `gmx dump` on each run's `.tpr` to recover that run's thermostat
  reference temperature, from an installation kept in a separate `gmx` conda
  environment.

Version floors live in `pyproject.toml` and `environment.yml`; the exact versions
the shipped artefacts were produced with are tabulated in
`THIRD_PARTY_NOTICES.md`.

## Installation

Four steps, on a clean clone, with no simulation data present.

1. Clone and enter the repository.

```bash
git clone https://github.com/djeckins/PCM_diploma_project.git
cd PCM_diploma_project
```

2. Create the environment. Conda is what the published numbers were produced
with; a plain virtual environment works too and is what CI uses.

```bash
conda env create -f environment.yml     # environment name: chem
conda activate chem
```

Without conda: `python -m venv .venv && source .venv/bin/activate`, then
`python -m pip install -e ".[dev]"`, which pulls the same libraries at the same
declared floors and finishes step 3 as well.

3. Install the package from *this* clone.

```bash
pip install -e . --no-deps
```

Expect `Successfully installed pcm2-0.1.0`. The install must be editable and must
be made from the clone: `external/` is not packaged into a wheel and the feature
layer reads a file from it, so a non-editable install cannot find its vendored
data. If another checkout of `pcm2` is already installed in this environment it
shadows this one for `import pcm2` and for `python -m pcm2.run` — installing from
here fixes that, and `PYTHONPATH=src` in front of a command does too. `pytest` is
unaffected, because `pyproject.toml` puts `src` first for it. To see which copy
is live:

```bash
python -c "import pcm2, pathlib; print(pathlib.Path(pcm2.__file__).parent)"
```

4. Check the installation.

```bash
python -m pytest
```

Expect `84 passed, 39 warnings in 2.3s`. The warnings are library deprecations,
not failures. The suite builds its own synthetic inputs or reads the committed
tables — nothing in `tests/` opens a trajectory, so it passes with
`PCM2_DATA_ROOT` pointing at a path that does not exist. (`make setup && make
test` does steps 3 and 4 together, but the `Makefile` hard-codes a `~/miniforge3`
conda hook.)

For the notebooks, once:

```bash
conda install -c conda-forge jupyterlab   # or: pip install jupyterlab
jupyter lab                               # start it from the repository root
```

Choose the `chem` environment as the kernel.

## Using the monitor

Three notebooks, three questions. All are started from the repository root, and
none of them needs a trajectory: they read the committed tables under `runs/`.

| notebook | the question it answers | when to open it |
|---|---|---|
| `notebooks/permeation_predict.ipynb` | For a chosen protein, moment and horizon: is this channel ready to conduct, and if not, which mechanism is holding it back? | whenever you want a verdict rather than a summary — on one of the five shipped systems, or on a protein of your own once it has a config, an `autodetect` run and an `events` run |
| `notebooks/methods_showcase.ipynb` | At one moment, what does each benchmarked method report in its own native units, and what does each claim to predict? | when placing this monitor beside the published alternatives, or to regenerate the showcase panel |
| `notebooks/run_my_trajectory.ipynb` | What did one pipeline run actually produce, step by step? | straight after running the pipeline on a system of your own. It recomputes nothing, so it cannot overwrite a finished run |

In `permeation_predict.ipynb` you edit only the first code cell — protein,
condition, replica, the `full` or `physics_only` variant, the time point, the
time range, and the horizons to sweep. One run answers all three in order: the
measured pore state at *t*, the verdict at τ with its mechanism ranking and the
same frame at the other cached horizons, then the READY fraction over the range
with its probability timeline. What the first cell prints as it ships
(gramicidin, `298K_2Msalt`, replica `01`, t = 100 000 ps, τ = 1000 ps,
`physics_only`):

```
t = 100000 ps | horizon τ = 1000 ps | P(crossing within τ) = 0.022 → NOT READY
  1. constriction dehydration — weight 57% (driver: flu_r_constr_var_win = 0.0908)
  2. pathway/lining rearrangement — weight 19% (driver: geo_bin7_nsearch = 2)
  3. permeant delivery — weight 17% (driver: dlv_t_since_cross_ps = 2.12e+04)
```

The first request for a (τ, variant) pair that is not already cached in
`runs/pooled-predict/` trains the pooled model, which takes roughly a minute.
`notebooks/predict_lib.py` holds everything the cells call, so it can be tested:
`python notebooks/verify_predict.py` re-checks the prediction path against the
pipeline's own labels, printing one line per horizon and one per time-range case
marked `OK` or `MISMATCH`, and asserting on the mismatches. Run it after any
change to that path. `notebooks/README.md` describes each notebook in full.

## Running the pipeline on your own trajectories

The simulation data is needed only for the measuring steps.

1. Point `PCM2_DATA_ROOT` at the directory holding the datasets
   (`docs/DATA.md` describes the layout each config expects; see also
   [Where the data came from](#where-the-data-came-from)).

```bash
export PCM2_DATA_ROOT=/absolute/path/to/your/data
```

2. Regenerate the configs from the files on disk. Discovery walks the data root,
writes a fresh `configs/<system>.yaml`, then loads it under strict validation and
prints one line per system.

```bash
python tools/init_configs.py
```

3. Run the pipeline for one system.

```bash
python -m pcm2.run all --config configs/mthk.yaml
```

The steps run in dependency order, each printing its own banner (`=== step
autodetect (mthk) ===`):

```
autodetect → events → features → coords → labels → train → figures → report → benchmark
```

Any single step can be run on its own with the same command. `autodetect` must
come first, because nothing can be measured before it has fixed the atom
selections, the pore axis and extent, the coordination cutoff and the
architecture profile. Each step recomputes and replaces its directory
`runs/<system>-stride<N>/<step>/` wholesale — there is no caching, so stale
intermediate state cannot leak into a run — and leaves a `PROVENANCE.json` and a
log beside its output.

`docs/REPRODUCTION.md` gives the full sequence that produced everything under
`runs/`: the five per-system pipelines, the pairwise transfer runs, the
cross-protein layer, then the figures and tables.

### Why the shipped paths read `/path/to/data`

The configs and the `runs/**/PROVENANCE.json` records were written on the machine
the runs were made on, and the absolute dataset paths in them were replaced by
the placeholder `/path/to/data` before publication. They are kept as records
rather than deleted: a provenance file states which input file produced which
artefact, and the directory names below the placeholder are part of that record.
They will not resolve on your machine and are not meant to — step 2 above writes
real paths. `configs/templates/` holds the same five configs with
`${PCM2_DATA_ROOT}` in place of the data root; nothing in `src/pcm2` expands
environment variables inside a config, so a template used as-is will load,
validate and then fail to find its trajectories. Regenerating with
`tools/init_configs.py` is the preferred route;
`configs/templates/README.md` gives the manual one.

## Training, scoring a new protein, and retraining

There are two models, and they are trained in different places.

The per-system model is the `train` step of the pipeline above. It fits the
linear and boosted arms for one system with leave-one-unit-out folds, writes
`runs/<system>-stride<N>/train/` (out-of-fold scores, calibration, folds,
evaluation, the fitted model, the alarm and the mechanism specificity), and is
what the in-system numbers report.

The pooled monitor is what the notebooks apply and what every cross-protein
result is measured on. It is trained on all conducting units at once and cached:

```bash
cd notebooks
python _train_pooled.py 2000 physics_only     # <tau_ps> <variant>
```

It prints the fitted threshold, the positive count and base rate, and the
out-of-fold lift, and writes `runs/pooled-predict/tau2000-physics_only-v3.pkl`.
`predict_lib.load_pooled_model(..., refresh=True)` forces a refit; the notebooks
read the current cache generation only, `v3`, which is the generation shipped
here. The `physics_only` variant drops the columns whose per-system value ranges
do not overlap, since a non-overlapping column can act as a hidden system label.

### Scoring a protein the model has never seen

The prediction path is generic: descriptors are measured from the trajectory and
the frozen pooled model is applied by column name. Adding a protein *as a
prediction target* needs no code change — it is a YAML config plus two steps:

1. `python tools/init_configs.py` — or write `configs/<id>.yaml` by hand:
   `system.id`, `data.data_root`, the conditions with their replica files and
   lineages, and the horizon keys.
2. `python -m pcm2.run autodetect --config configs/<id>.yaml` — fills the atom
   selections, pore axis and extent, coordination cutoff and architecture
   profile. Nothing can be measured before this.
3. `python -m pcm2.run events --config configs/<id>.yaml` — required, because the
   delivery descriptors are built from detected crossings.
4. Optional: `python -m pcm2.run features --config configs/<id>.yaml`, which lets
   the notebook read rows from a table instead of measuring them on demand.

The target needs no labels: the monitor does not label what it predicts. An
off-family verdict carries the zero-shot transfer quality measured in
`runs/monitor-generalisation-tau*/`, not in-system accuracy — and, on the
evidence below, ranking above chance off family does not imply a usable alarm.

### Extending the training pool

The pool is deliberately frozen at the three conducting proteins, so widening it
*is* a code change — two constants in `tools/monitor_generalisation.py`:

```python
SYSTEM_CONFIGS = ["configs/gramicidin.yaml", …]        # which configs are loaded
CONDUCTING = {("gramicidin", None), ("mthk", None), ("kcsa", "E71A")}
```

A `(system, None)` entry admits every condition of that system; a
`(system, condition)` entry admits one condition, which is how the conducting
E71A arm of KcsA is separated from its non-conducting siblings in the same
config. Add your system's config to the first and its conducting units to the
second, run the pipeline for it through at least `labels`, then refit the pooled
model. The cross-protein layer is re-run with:

```bash
PCM2_GEN_TAU=2000 python tools/monitor_generalisation.py
```

`PCM2_POOL_ADD_CX43=1` reproduces the recorded pool-inclusion ablation into a
`-cx43pool` directory, and `PCM2_INNER=protein` the rejected
protein-inner-selection check, leaving the canonical runs untouched. Only
conducting systems can teach the positive class, and three architectures is the
minimum that makes leave-one-protein-out meaningful. Before adding a fourth,
read the pool-inclusion result under "Measured and then rejected" below: adding
Cx43 helped nothing and cost gramicidin.

## What this repository contains

| path | what it is | size |
|---|---|---|
| `src/pcm2/` | the pipeline: reading layer, autodetection, event detection, descriptors, labels, models, interpretation, figures, report. 42 modules | 0.50 MB |
| `tools/` | 13 scripts. Eleven turn run artefacts into further artefacts: the cross-protein layer, the assembled results, the per-unit reports, the write-up figures and tables, the adjudicated side experiments, and `genblocks.py`, which rewrites the generated blocks in `docs/`. The other two read no run artefact — `init_configs.py` walks the data root to write `configs/`, and `workflow_figure.py` draws the schematic from its own source. See `tools/README.md` | 0.18 MB |
| `tests/` | 19 test modules, 84 tests. They run without any trajectory data | 46 kB |
| `configs/` | one YAML per system (5), each key carrying a recorded basis for its value, plus `configs/templates/` — the same five in portable `${PCM2_DATA_ROOT}` form | 0.17 MB |
| `docs/` | `METHODS.md` (the canonical description of the method), `RESULTS.md`, `REPRODUCTION.md`, `SCHEMA.md`, `DATA.md`. The numeric blocks are generated by `tools/genblocks.py`, not typed | 0.08 MB |
| `notebooks/` | three notebooks and the library they call | 0.55 MB |
| `external/` | two data files vendored byte-for-byte from the CHAP project, with checksums, provenance and the upstream licence notice. Used by the Rao-2019 comparison arm *and* by the feature layer | 0.60 MB |
| `runs/` | 39 directories of run artefacts: per-system pipeline runs, the cross-protein experiments, the negative controls, per-unit reports, and the assembled results. 437 figures, 72 parquet tables | 135 MB |

1 039 files and 137 MB in total (decimal: 1 MB is 10⁶ bytes), of which `runs/` is
135 MB. There are no trajectories and no coordinate files of any kind; a CI step
fails the build if any appear.

Where to look for what:

- `runs/RESULTS/<role>/` — the assembled results, one folder per narrative role
  (`training_*`, `control_negative_*`, `test_*`, `cross_protein`), each with a
  `README.md` carrying its numbers. Rebuilt with
  `python tools/collect_results.py`.
- `runs/FINAL/` — the figures and tables of the write-up, regenerated by
  `python tools/final_results.py` (the workflow schematic and the worked SHAP
  example have their own scripts; `cx43_arms.py` must run before
  `final_results.py`).
- `runs/<system>-stride<N>/` — one subdirectory per pipeline step, each with its
  own `PROVENANCE.json` and log. `benchmark/benchmark.parquet` is the
  authoritative comparison table; the `benchmark` step also writes a CSV of the
  same rows, which is not kept here because nothing reads it.
- `runs/monitor-generalisation-tau{1000,2000,4000}[-cx43pool|-pinner]/` — the
  cross-protein experiments and their `GENERALISATION.md`.
- `runs/postproc-check/`, `runs/hyperflat-check/`, `runs/plsfma-selfcheck/`,
  `runs/transfer-*/`, `runs/pooled-loco-tau1000/`, `runs/unit-report-*/` — the
  adjudicated side experiments, the published-method self-check and the per-unit
  reports.

The `*.pkl` files under `runs/` are scikit-learn and XGBoost artefacts. Python
pickles execute code on load, so load them only from a checkout you trust.

## Systems studied

One config per system, in two roles.

- Training pool — three conducting proteins: `gramicidin` (gA, KCl, 500 mV),
  `mthk` (MthK, ±200 mV), and from the KcsA family only the conducting K⁺ unit
  `kcsa:E71A` (300 mV), where crossings are frequent enough to learn from.
- Never in the pool — five KcsA filter-mutant arms (K⁺: G77A·E71A,
  T75A·E71A; Na⁺: E71A, G77A·E71A, T75A·E71A) and `cx43` (connexin-43 gap
  junction, 200 mV), an unrelated architecture with a wide aqueous pore, two
  membranes and rare permeation. Of the mutant arms, `kcsa_na:G77AE71A` leaks
  enough to rank, with seven complete crossings; the K⁺ G77A·E71A arm carries
  three isolated crossings and the Na⁺ T75A·E71A arm one, and the remaining two
  carry none. The two K⁺ arms are held out of the *pool* only — they are units of
  `configs/kcsa.yaml`, so they do appear inside the per-system KcsA model.

The horizon is part of the question. τ is not a hyperparameter to be tuned.
To be answerable it must be commensurate with the channel's own kinetics, so the
declared sweep τ = 1, 2 and 4 ns brackets the measured median transit times of
three of the four systems — gramicidin A 1.36 ns, MthK 1.95 ns, Cx43 4.38 ns
(`runs/FINAL/tableA_proteins.md`); τ = 2 ns was fixed in advance as primary.
KcsA is the exception: the median transit of its conducting E71A arm is 36.5 ns,
far longer than any horizon it is scored at. Changing τ changes the answer for a
given frame, not just its precision, so the horizon must always travel with the
number.

## Results at a glance

The headline metric is average precision divided by the base rate ("×chance"),
because absolute precision is not comparable across systems whose base rates run
from 0.2 % (the KcsA frames the tree head is scored on) to 50 % (MthK).

1. In-system, one model per protein (τ = 2000 ps,
`runs/FINAL/tableD_method_matrix.md`):

| method | gramicidin A | MthK | KcsA-E71A |
|---|---|---|---|
| Rao 2019 Σ_d (published structure rule) | ×0.97 | ×1.05 | ×0.64 |
| PLS-FMA (independent implementation) | ×1.04 | ×0.95 | ×0.78 |
| Linear (ridge-logistic) | ×2.66 | ×1.16 | ×0.57 |
| GBT (this work) | ×2.59 | ×1.37 | ×1.98 |
| Clock (structure-free temporal null) | ×1.53 | ×1.01 | ×0.55 |
| GBT zero-shot (trained on the other conducting proteins) | ×1.60 | ×1.14 | ×1.21 |
| Base rate (chance) | ×1.00 | ×1.00 | ×1.00 |

Gramicidin A rests on the most evidence: of 12 510 frames analysed, 12 385 carry
a defined label at τ = 2000 ps and are the frames every arm is scored on, with
1 357 positives, base rate 0.1096, AP 0.284, ROC-AUC 0.798, MCC 0.302 and all 10
folds fitted (`runs/FINAL/tableE_ml_quality.md`,
`runs/RESULTS/training_gramicidin/`). The table above is read at one shared
horizon, τ = 2000 ps; `docs/RESULTS.md` and the `BENCHMARK.md` files instead read
each system at the horizon its own config declares primary — gramicidin A and
MthK at 2000 ps, both KcsA arms at 1000 ps, Cx43 at 4000 ps — so the same KcsA
tree arm reads ×1.98 here and ×3.14 there, differing only by horizon. Three
cautions belong with this table, and `docs/RESULTS.md` carries them in full: on
gramicidin the tree head is not
distinguishable from the linear arm (paired bootstrap over trajectories
Δ(AP) = −0.0083, CI [−0.042, 0.022]), while the trees do win on MthK (+0.1031,
CI [0.068, 0.127]); the KcsA column rests on 37 positive frames at a base rate of
0.00198 with only 2 of 3 folds fitted, because no training rotation for the E71A
fold produced both classes; and that refusal means the GBT figure is not scored
on the same frames as the rest of its column. "KcsA-E71A" names the whole K⁺
config — a leave-one-arm-out fit over its three arms — not a model of E71A alone.

2. Zero-shot: leave one protein out of the pool. Median, over the held-out
protein's trajectories, of within-trajectory AP divided by that trajectory's own
base rate (`runs/FINAL/tableB_lopo4.md`):

| held-out protein | τ = 1000 | τ = 2000 | τ = 4000 | physics-only at τ = 2000 |
|---|---|---|---|---|
| gramicidin A (10 units) | ×1.46 | ×1.60 | ×1.46 | ×1.73 |
| MthK (4 units) | ×1.39 | ×1.14 | ×1.06 | ×1.13 |
| KcsA-E71A family (1 unit) | ×1.22 | ×1.21 | ×1.16 | ×1.31 |

The per-unit spread matters more than the median: at τ = 2000 the ten gramicidin
units span ×1.02 to ×2.11 and the four MthK units ×1.08 to ×1.36, so the MthK
range nearly touches chance, and KcsA-E71A contributes exactly one unit. The
physics-only variant is the anti-leakage control — for this pool it removes 21 of
the 100 columns the model is fed and still matches or beats the full set in four
of nine cells, so the transfer lift is not simply the model recognising which
protein it is looking at.

3. Negative controls: silence, and one ranked leak. The frozen three-protein
monitor applied to five KcsA arms it never saw
(`runs/FINAL/tableC_negative_controls.md`). Over 43 755 frames the READY fraction
is 0.0000 on all five arms at τ = 1000 and τ = 4000 and 0.0000–0.0001 at
τ = 2000, with mean P(ready) far below the fitted threshold at every horizon. The
one arm that does leak, `kcsa_na:G77AE71A/1` — 7 crossings in 12 501 frames — is
ranked ×17.66 over its own chance at τ = 1000 (×13.99 physics-only), ×4.41 at
τ = 2000 and ×3.30 at τ = 4000. Both halves are needed: silence alone would be
satisfied by a model that always says NOT READY, and the ranking shows the score
still tracks something physical.

4. Off-family: connexin-43. The same frozen monitor on an unrelated
architecture reads ×1.17 at τ = 1000, ×2.01 at τ = 2000 and ×2.8–×3.9 at
τ = 4000. The range at that horizon is the spread between two different
experiments on the same unit: ×3.91 when the three-protein monitor is applied to
Cx43 as a negative control (`runs/monitor-generalisation-tau4000/`), and ×2.76
when Cx43 is instead put inside a four-protein pool and held out from it
(`runs/monitor-generalisation-tau4000-cx43pool/`, the recorded pool-inclusion
ablation). Two things must be read with it: it
rests on one trajectory unit, `cx43:200mV/R2`, with 1 001 frames and two
crossings detected system-wide; and on that unit the baselines beat the monitor —
the clock null reaches ×5.2 at τ = 2000 against the GBT's ×2.01, and at τ = 4000
the linear arm reaches ×7.93 (`runs/FINAL/cx43_arms.json`).

5. Event-level forewarning, declared in advance as the primary quantity: the
fraction of real crossings preceded by an alarm under hysteresis, against a
structure-matched random alarm of the same duty (`runs/*/train/alarm.json`).

| protein | crossings warned | random alarm of same duty | empty alarm episodes |
|---|---|---|---|
| gramicidin A | 0.350 | 0.111 | 12 of 40 |
| MthK | 0.068 | 0.039 | 5 of 26 |
| KcsA-E71A | 1.000 | 1.000 | 0 of 1 |

Gramicidin is the result. MthK's margin is weak. The KcsA row rests on one
always-on alarm episode against a null that also scores 1.000; it is
uninformative and must not be quoted as a success.

6. Mechanism diagnosis. The verdict is the mechanism axis contributing most
toward "not ready", from TreeSHAP contributions of the already-trained model
aggregated onto seven physical axes. Axis specificity is published as an odds
ratio — fires on not-ready frames against ready frames — *before* the winner is
chosen (`runs/*/train/mechanism_specificity.json`). On gramicidin A: ion
placement and desolvation 3.82, column immobility 3.09, constriction dehydration
2.63, water-wire rupture 2.51, permeant delivery 2.17, pathway/lining
rearrangement 1.83. On MthK the axes fire on 81–100 % of *all* frames with odds
ratios between 0.58 and 2.17, several below 1, so specificity is essentially
absent there. Two wording rules follow, and they are not stylistic:
explanations are associational, so a mechanism is described as *enriched among
non-conducting frames*, never as the cause of the block, and a mechanism share is
never quoted without its firing rate on conducting frames.

### Measured and then rejected

The rejections are results too, and the evidence for each is archived.

- Adding Cx43 to the training pool. Held-out medians move by at most 0.04 in
  ratio for MthK and KcsA and cost gramicidin 0.07–0.16 depending on horizon.
  Nothing gains, so Cx43 stays outside the pool
  (`runs/FINAL/tableC2_pool_choice.md`).
- Causal score smoothing. It inflates per-frame AP (×3.91 → ×5.72 on Cx43)
  but at matched event-level duty it is neutral on Cx43 and harmful on
  gramicidin (0.327 → 0.286 at 5 % duty) and KcsA (0.367 → 0.327). Raw
  per-frame AP stays canonical
  (`runs/postproc-check/event_level_lopo4.json`, written by
  `tools/event_level_check.py`; the adjudication in
  `runs/postproc-check/VERDICT.md` quotes an earlier run of the same test, kept
  as evidence, whose figures differ slightly and whose conclusion does not).
- Transfer-aligned hyperparameter selection. The inner criterion is not flat
  and is misaligned with transfer: its best points give the worst gramicidin
  transfer (`runs/hyperflat-check/FLATNESS.md`). The grids are frozen.
- Protein-level inner selection. Rejected on the negative controls, where it
  fired READY on 0.877 and 0.999 of two non-conducting arms
  (`runs/monitor-generalisation-tau2000-pinner/GENERALISATION.md`).

## Limitations

State these alongside any result taken from this repository.

- Sample size. Four proteins in five configs, 22 trajectory units, of which
  15 are in the training pool (10 gramicidin, 4 MthK, one KcsA-E71A); the pooled
  training matrix at τ = 1000 physics-only is 26 672 rows with 3 243 positives
  (12.2 %). Only gramicidin and MthK are genuinely replicated. The KcsA
  in-system figure rests on 37 positive frames, the Cx43 transfer figure on one
  trajectory unit with one labelled crossing.
- Ranking, not thresholded decision. On the rare-positive systems the
  calibrated verdict is not usable (KcsA: MCC −0.001, F1 0.0). AP/chance is the
  metric to quote, and it is a statement about ranking.
- The alarm does not transfer. Off family the monitor ranks pre-permeation
  frames above chance but the calibrated alarm never fires: 0 of 7 crossings on
  the leaking Na⁺ arm, 0 of 1 on the Cx43 unit.
- The negative controls cannot prove specificity. Two of the five arms have
  no positive frames, so AP is undefined; the Na⁺ T75A·E71A arm ranks *below*
  chance (×0.65). Four of five arms are too event-poor to score. The silence is
  also a property of this fit rather than of the method: in the physics-only
  variant at τ = 4000 the READY fraction rises to 0.027 on the Na⁺ G77A arm.
- Baselines are not always beaten. On Cx43 the clock null and the linear arm
  score higher; on gramicidin the tree head is indistinguishable from the linear
  arm.
- Explanations are associational. The mechanism axis is what moved this
  verdict, not a demonstrated cause, and its specificity is protein-dependent.
- Not every acceptance check passes, by design. Cx43 fails four checks
  because it has too few crossings for an in-system model; the K⁺ KcsA run fails
  one and the Na⁺ run two, for the same reason. Only gramicidin and MthK are
  fully green (`docs/REPRODUCTION.md`).
- The readiness target is a horizon label, not a committor. `y(t) = 1` iff at
  least one complete crossing finishes in (t, t+τ]; a frame whose horizon window
  is incomplete is *unknown*, not negative. No Markov-state-model or committor
  estimate is computed anywhere here, and a principled data-driven definition of
  readiness remains open work.
- Descriptor counts differ by context and all three are correct: the schema
  documents 104 columns, of which 100 are fed to the model (79 physics-only,
  after dropping 21), and the shipped feature table adds condition, replica and
  time_ps for 107.

## What is not here, and what that costs

- The trajectories and topologies. No `.xtc`, `.trr`, `.tpr`, `.top`,
  `.itp`, `.gro` or `.pdb` file is present; the configs reference them by path
  only. They amount to tens of gigabytes and come from three sources, not all of
  which permit redistribution — see [Where the data came from](#where-the-data-came-from).
- `runs/*/coords/coords.parquet` — the per-frame superimposed Cα coordinate
  tables, 468 MB across the five systems as measured on the source runs before
  they were excluded (170 MB for each KcsA arm, 77 MB Cx43, 45 MB MthK, 7 MB
  gramicidin A; the frame and column counts behind those sizes are recorded in
  each `runs/*/coords/log.txt`). These are closer to a re-encoding
  of the input trajectory than to a derived measurement, so they are
  `.gitignore`d. The `coords/` directories themselves ship, with their
  `PROVENANCE.json` and log, and the `coords` step regenerates the table locally
  once the trajectories are in place.

Two things therefore cannot be reproduced from the repository as shipped.
`tools/plsfma_selfcheck.py` reads `coords/coords.parquet` directly, so it cannot
be re-run, although its output directory is committed and readable. And the
`train` step records its `plsfma_coords` comparison arm as skipped rather than
failing when the coordinate table is absent, so a pipeline re-run without
trajectories produces a quietly smaller method comparison than `docs/RESULTS.md`
shows. Separately, `runs/hyperflat-check/` is evidence rather than a reproducible
step: its generating script was not kept.

Everything else works without the trajectories, deliberately: the feature, label
and event tables ship for all five systems, which is why the per-protein results,
the reports and every figure can be read as they are, why the prediction notebook
runs in table mode, and why the full test suite passes with no data at all.

## Reproducibility

- Generated numbers. Every numeric block in `docs/` is written by
  `python tools/genblocks.py` from run artefacts, never typed by hand.
  `--check` verifies that documentation and artefacts still agree, and that check
  is part of the test suite.
- Vendored data integrity. `tests/test_vendored_checksums.py` recomputes the
  checksums of the vendored CHAP files on every run, so a silent substitution
  fails the suite.
- Continuous integration. `.github/workflows/tests.yml` installs the package,
  runs the suite with `PCM2_DATA_ROOT` set to a non-existent path, and fails the
  build if any trajectory or topology file has been committed.
- Command sequence. `docs/REPRODUCTION.md` gives the main line and states
  what it does not regenerate; `tools/README.md` documents every script and its
  output.

`CONTRIBUTING.md` describes what a useful change looks like and the invariants
the test suite exists to protect.

## Licence

The code — `src/`, `tools/`, `tests/`, `notebooks/`, `configs/`, the `Makefile`
and the build metadata — is released under the MIT licence (`LICENSE`). The
figures, tables, reports and derived data under `runs/` and `docs/` are released
under CC-BY-4.0 (`LICENSE-DATA.md`), which fits scientific artefacts better than
a software licence and expressly covers database rights in structured tables.

Neither licence covers `external/rao2019_heuristic/*.json`, which is third-party
material under its own licence, nor the input trajectories, which are not
contained here and over which no licence is claimed.

## Citation

`CITATION.cff` carries the machine-readable form. In text:

> Ivanova, A. (2026). *PCM 2 — permeation-readiness forecasting for ion
> channels* (version 0.1.0) [software and derived data].
> <https://github.com/djeckins/PCM_diploma_project>

The work accompanies an MSc dissertation at Queen Mary University of London;
please cite the dissertation as well if you draw on its findings. If you use the
trajectory datasets, cite them at their own DOIs as listed under
[Where the data came from](#where-the-data-came-from): a derived analysis does not replace the
deposit's own citation.

## Third-party components

`THIRD_PARTY_NOTICES.md` is the full attribution file. In brief:

Redistributed. `external/rao2019_heuristic/` holds two data files taken
byte-for-byte from CHAP (<https://github.com/channotation/chap>) — the heuristic
grid of the Rao et al. 2019 hydrophobic-gate criterion and CHAP's encoding of the
Wimley–White 1996 interfacial hydrophobicity scale. CHAP is MIT-licensed; the
upstream notice is reproduced verbatim in
`external/rao2019_heuristic/LICENSE.CHAP.md`, with checksums and the reason for
vendoring in its `PROVENANCE.md`. These files are load-bearing for the feature
layer as well as for the comparison arm, so removing them breaks the pipeline.

Reimplemented from publications, with citation. The Rao-2019 criterion, the
CHAP pore-radius and pore-facing-residue definitions, and PLS-FMA. The PLS-FMA
arm is written independently in NumPy from the published description and verified
against scikit-learn's PLS1 to numerical precision; no source from the reference
tool is copied. Its self-check reaches a best canonical correlation of 0.503 on
KcsA, 0.420 on MthK, 0.389 on the KcsA Na⁺ arm, 0.138 on Cx43 and 0.051 on
gramicidin, so the implementation is validated for the tetramers and not for
gramicidin.

Imported, not redistributed. The libraries listed at the top. MDAnalysis is
the only copyleft dependency and is weak copyleft (LGPL) from version 2.8 onward,
which is why that is the declared floor.

## Where the data came from

The molecular-dynamics trajectories analysed here are not contained in this
repository and are not redistributed by it. The five configs cover four proteins
— the KcsA filter mutants were simulated twice, once in K⁺ and once in Na⁺, and
each ionic condition is its own config — and they come from three sources:

| system | protocol as recorded | source | may be redistributed? |
|---|---|---|---|
| MthK, 2 replicas at +200 mV and 2 at −200 mV | CHARMM36, charge scaling 0.78, 323.15 K | the project supervisor at Queen Mary University of London | no |
| Gramicidin A, 5 replicas at 298.15 K / 2 M KCl and 5 at 330 K / 1 M KCl, 500 mV | charge scaling 0.70 and 0.75 respectively; force field not recorded in the config | the project supervisor at Queen Mary University of London | no |
| KcsA filter mutants, K⁺ — E71A, G77A·E71A, T75A·E71A, one replica each | CHARMM36m, no charge scaling, 300 mV, 290 K | Zenodo [10.5281/zenodo.12623504](https://doi.org/10.5281/zenodo.12623504) | yes, CC-BY-4.0 at the deposit |
| KcsA filter mutants, Na⁺ — the same three variants | CHARMM36m, no charge scaling, 300 mV, 290 K | the same deposit | yes, CC-BY-4.0 at the deposit |
| Connexin-43, replicas R2 and R4 at 200 mV | CHARMM36, no charge scaling, 303.15 K | Zenodo [10.5281/zenodo.8191584](https://doi.org/10.5281/zenodo.8191584) | yes, CC-BY-4.0 at the deposit |

Gramicidin A and MthK were provided by the project supervisor at Queen Mary
University of London and are used with permission. Permission to use is not
permission to redistribute: they are not published, linked or re-hosted here, and
no licence over them is granted or implied. Enquiries about access should go to
the supervising laboratory. Unlike the two deposits below, this entry rests on
the project record rather than on anything a reader can check here. For
gramicidin the `*-transit-times.dat` files that came with the trajectories are a
pre-existing crossing annotation, used only as a cross-check: the config records
`events.source: own`, so every label in this work comes from this project's own
detector, and `runs/gramicidin-stride4/events/summary.json` records `n_own` and
`n_provided` side by side.

The KcsA filter-mutant arms come from the Zenodo deposit *Data from Selectivity
filter mutations shift ion permeation mechanism in potassium channels*,
A. Mironenko, B. L. de Groot and W. Kopec, 2024, CC-BY-4.0; the accompanying
paper is *PNAS Nexus* 2024, 3(7), pgae272,
[10.1093/pnasnexus/pgae272](https://doi.org/10.1093/pnasnexus/pgae272). The three
variants are the mutants that paper studies, which is why they enter this work as
engineered non-conducting controls rather than as training material.

Connexin-43 comes from the Zenodo deposit *Molecular dynamics simulation data 1:
Structure of the connexin-43 gap junction channel in a putative closed state*,
S. Acosta-Gutierrez and F. L. Gervasio, 2023, CC-BY-4.0; the accompanying paper
is C. Qi, S. Acosta Gutierrez, P. Lavriha *et al.*, *eLife* 2023, 12, RP87616,
[10.7554/eLife.87616](https://doi.org/10.7554/eLife.87616). Two things about this
set are recorded in the `origins` block of `configs/cx43.yaml`. The topology used
here is a GROMACS `.top` assembled from the deposit's own files, because the
deposited `topol.top` water count predates the final solvent trim. And the
deposit ships three production replicas of which the third is not independent:
`protein_center_200mV_R3.xtc` is a byte-identical copy of the R2 file — same
size, same MD5, identical detected crossings — so it is excluded and the two
genuine replicas, R2 and R4, are analysed.

Cite the deposits and their papers as the source of the simulations, and this
repository for the descriptors, labels, models and results computed from them.
Both deposits are CC-BY-4.0, which is what makes publishing derived tables
legitimate as long as they are credited; titles, authorship and licence terms
were checked against the Zenodo records and the papers against PubMed on
2026-08-17.

`docs/DATA.md` carries the technical side: trajectory spans and strides, the
deposit archives, the directory layout each config expects, how to run the
pipeline against a local copy, and what works with no trajectories at all.
