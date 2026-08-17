# REPRODUCTION — how to repeat

```bash
conda activate chem            # environment from environment.yml
make setup                     # pip install -e . --no-deps
make test                      # the whole test suite must be green
python tools/init_configs.py   # configs from files on disk + declared keys

# 1. per-system pipelines
for S in mthk gramicidin kcsa kcsa_na cx43; do
    python -m pcm2.run all --config configs/$S.yaml
done

# 2. pairwise transfer (source model applied to a target system)
python -m pcm2.run transfer --config configs/gramicidin.yaml \
       --source-config configs/mthk.yaml
python -m pcm2.run transfer --config configs/kcsa.yaml \
       --source-config configs/gramicidin.yaml
python -m pcm2.run transfer --config configs/kcsa.yaml \
       --source-config configs/mthk.yaml
python -m pcm2.run transfer --config configs/kcsa_na.yaml \
       --source-config configs/kcsa.yaml

# 3. cross-protein layer: leave-one-protein-out, negative controls, ablations
for T in 1000 2000 4000; do
    PCM2_GEN_TAU=$T python tools/monitor_generalisation.py
    PCM2_GEN_TAU=$T PCM2_POOL_ADD_CX43=1 python tools/monitor_generalisation.py
done

# 4. write-up figures and tables
python tools/shap_example_figure.py
python tools/workflow_figure.py
python tools/cx43_arms.py         # -> runs/FINAL/cx43_arms.json
python tools/final_results.py     # -> runs/FINAL
python tools/collect_results.py   # -> runs/RESULTS

python tools/genblocks.py      # last, once the artifacts are final
```

`tools/cx43_arms.py` has to run before `tools/final_results.py`: the method
matrices read `runs/FINAL/cx43_arms.json` and leave the Cx43 group empty when
that file is absent, and the arms are skipped rather than estimated.

`runs/FINAL/figH_methods_showcase.png` and the one-moment table come from
`notebooks/methods_showcase.ipynb`, which is run interactively and its panel
exported by hand; every other figure and table above is written by a script.
`runs/FINAL/figG2_shap_mthk.png` is written by a script but by no command line
above: it is `tools/shap_example_figure.py` with its unit and output path passed
to `main()` as overrides. `runs/hyperflat-check/` records a measured decision
whose generating script was not kept, so that one directory is evidence rather
than a reproducible step.

Each step recomputes and replaces its folder under `runs/<system>-stride<N>/`
wholesale; there is no cache and no partial state. Each step's provenance is
the `PROVENANCE.json` next to its output. A run and source edits never overlap
in time; tests are not run in parallel with a run.

## What the sequence above does and does not regenerate

It regenerates the main line: the five per-system pipelines, the pairwise
transfer runs, the cross-protein layer at three horizons with its pool-inclusion
ablation, and the write-up figures and tables. It does not regenerate the
adjudicated side experiments, which have their own entry points and are listed in
`tools/README.md`: the per-unit reports (`unit_report.py`), the post-processing
and event-level checks (`postproc_check.py`, `event_level_check.py`), the PLS-FMA
self-check (`plsfma_selfcheck.py`), the earlier unit-level leave-one-channel-out
run (`pooled_loco.py`), and the rejected protein-inner-selection variant of the
cross-protein layer (`PCM2_INNER=protein python tools/monitor_generalisation.py`,
which writes the `-pinner` directories). So a reader who follows the sequence and
then compares directory listings will find those extra directories in `runs/`;
they are evidence that was kept, produced by the scripts named above, not output
the main line forgot to write.

## Without the trajectories

Nothing above works without the trajectory datasets, which are not part of this
repository — see `docs/DATA.md` for what they are, where they came from and how to
point the pipeline at a local copy. Three things can still be checked on the
shipped tree, with no data at all, and they are worth doing first:

```bash
conda env create -f environment.yml && conda activate chem
pip install -e . --no-deps            # what `make setup` runs
python -m pytest                      # 84 tests, a few seconds
python tools/genblocks.py --check     # documents still agree with the artifacts
```

Use these plain commands rather than `make setup && make test` if you do not have
conda installed at `~/miniforge3`, which is the path the `Makefile` activates.

The test suite needs neither the datasets nor `PCM2_DATA_ROOT`; it exercises the
measurement, model and reporting layers on synthetic inputs and on the shipped
artifacts. `genblocks.py --check` recomputes every generated block from the
artifacts under `runs/` and reports a mismatch instead of rewriting, which is the
quickest way to confirm that the numbers in `docs/` are the numbers in the run
directories. Beyond that, every run's report, figures and tables read as shipped,
and the prediction notebook runs in table mode; `docs/DATA.md` §4 lists what is
and is not possible in that state.

## Acceptance table (generated from artifacts)

<!-- BEGIN GENERATED: acceptance -->
### cx43-stride2
- edge_exception2_recorded: PASS
- no_ceiling_hits: PASS
- penalty_in_df_units: PASS
- calibration_recorded_per_fold: FAIL
- events_source_declared: PASS
- vendored_checksums: PASS
- figure_catalog: FAIL
- answers_with_mechanism: FAIL
- resting_frames_defined: FAIL
- applicability_mask_written: PASS

### gramicidin-stride4
- edge_exception2_recorded: PASS
- no_ceiling_hits: PASS
- penalty_in_df_units: PASS
- calibration_recorded_per_fold: PASS
- events_source_declared: PASS
- vendored_checksums: PASS
- figure_catalog: PASS
- answers_with_mechanism: PASS
- resting_frames_defined: PASS
- applicability_mask_written: PASS

### kcsa-stride4
- edge_exception2_recorded: PASS
- no_ceiling_hits: PASS
- penalty_in_df_units: PASS
- calibration_recorded_per_fold: PASS
- events_source_declared: PASS
- vendored_checksums: PASS
- figure_catalog: PASS
- answers_with_mechanism: PASS
- resting_frames_defined: FAIL
- applicability_mask_written: PASS

### kcsa_na-stride4
- edge_exception2_recorded: PASS
- no_ceiling_hits: PASS
- penalty_in_df_units: PASS
- calibration_recorded_per_fold: FAIL
- events_source_declared: PASS
- vendored_checksums: PASS
- figure_catalog: PASS
- answers_with_mechanism: PASS
- resting_frames_defined: FAIL
- applicability_mask_written: PASS

### mthk-stride5
- edge_exception2_recorded: PASS
- no_ceiling_hits: PASS
- penalty_in_df_units: PASS
- calibration_recorded_per_fold: PASS
- events_source_declared: PASS
- vendored_checksums: PASS
- figure_catalog: PASS
- answers_with_mechanism: PASS
- resting_frames_defined: PASS
- applicability_mask_written: PASS
<!-- END GENERATED -->
