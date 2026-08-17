# tools

Eleven of these thirteen scripts read run artefacts and write further
artefacts: most write a directory of their own under `runs/`, a few write one
named file into an existing one, and `genblocks.py` rewrites the generated
blocks inside the `docs/` files. The other two read no run artefact at all —
`init_configs.py` walks the data root under `PCM2_DATA_ROOT` to write
`configs/`, and `workflow_figure.py` draws the schematic from its own source
text. The pipeline itself lives in `src/pcm2`.

| script | writes | role |
|---|---|---|
| `init_configs.py` | `configs/*.yaml` | discovers the datasets on disk and assembles one config per system, with a recorded basis for every declared key |
| `genblocks.py` | the generated blocks in `docs/` | turns artefacts into the numeric tables of the documentation; `--check` mode is enforced by `tests/test_docs_generated_blocks.py` |
| `monitor_generalisation.py` | `runs/monitor-generalisation-tau*/` | the cross-protein layer: leave-one-protein-out over the conducting pool, the final model on the negative controls, and the pool-inclusion ablation. Switches: `PCM2_GEN_TAU`, `PCM2_POOL_ADD_CX43`, `PCM2_INNER` |
| `final_results.py` | `runs/FINAL/` | every figure and table of the write-up, in one command; refits nothing |
| `workflow_figure.py` | `runs/FINAL/figure_2_1_workflow.png` | the workflow schematic, drawn from a script so it can be regenerated |
| `shap_example_figure.py` | `runs/FINAL/figG_shap_example.png` | the worked example of a single verdict decomposed into mechanism axes |
| `cx43_arms.py` | `runs/FINAL/cx43_arms.json` | the comparison arms for connexin-43, which has no in-system model: the published rules frame by frame, and the fitted arms applied zero-shot from the conducting pool. Read by `final_results.py`, so run it first |
| `collect_results.py` | `runs/RESULTS/` | the curated per-role view of the raw run directories |
| `unit_report.py` | `runs/unit-report-*/` | a per-trajectory report for one held-out unit |
| `postproc_check.py` | `runs/postproc-check/results.json` | measures the post-processing candidates (score smoothing, arm ensembling, threshold headroom) |
| `event_level_check.py` | `runs/postproc-check/event_level_lopo4.json` | the event-level test at matched time-in-warning, which rejected smoothing |
| `plsfma_selfcheck.py` | `runs/plsfma-selfcheck/` | validates the PLS-FMA port in the regime its authors designed it for |
| `pooled_loco.py` | `runs/pooled-loco-tau1000/` | the earlier unit-level leave-one-channel-out experiment, superseded by `monitor_generalisation.py`; kept because its output directory is committed |

Beside `results.json` and `event_level_lopo4.json`, `runs/postproc-check/` holds
three files that neither script writes; they are kept as evidence rather than
regenerated. `VERDICT.md` is the note recording what was decided,
`event_level.json` an earlier run of the event-level test — the file that
`VERDICT.md` quotes — and `descriptor_change.json` the Cx43 zero-shot ratios
under the 83-column and the 104-column canon.

Not every directory under `runs/` is written from here. The per-system
`runs/<system>-stride<N>/` and the pairwise `runs/transfer-*/` come from the
pipeline itself (`python -m pcm2.run <step>` and `python -m pcm2.run transfer`),
and `runs/pooled-predict/` is the model cache of `notebooks/predict_lib.py`.
`runs/hyperflat-check/` records a measured decision (the hyperparameter
selection criterion is not flat and is misaligned with transfer); its
generating script was not kept, so the directory cannot be regenerated.
