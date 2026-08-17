# FINAL: the three-block thesis experiment

Block A — the four proteins on their own terms: `tableA_proteins`,
`figA_methods_per_protein.png` (six-way method comparison per protein).

Block B — train on three, zero-shot on the fourth (all four rotations):
`tableB_lopo4`, `figB1_zeroshot_timelines.png` (prediction vs reality at each
protein's own horizon), `figB2_lopo_medians.png` (per-unit spread).

Block C — the final three-protein monitor on the non-conducting mutant arms:
`tableC_negative_controls`, `figC_controls.png` (silence + the one real leak),
and `tableC2_pool_choice` (measured three-vs-four pool comparison behind the
canonical-pool decision).

Blocks D, E and F were added later and come from the same command. D is the
method-by-protein matrix: `figD_method_matrix.png` with each protein at its own
horizon, `figD_method_matrix_tau{1000,2000,4000}.png` at the three common
horizons, and `tableD_method_matrix`. E is `tableE_ml_quality`, the training
passport of each in-system model. F is `figF_training_quality.png`, the same
models in the standard ML views (ROC with AUC, precision-recall with PR-AUC,
reliability with the Brier score).

All numbers are read from canonical run artifacts (104-column canon); sources:
runs/monitor-generalisation-tau*{,-cx43pool}, per-system runs/<id>/train and
events. Every `table*` above is written twice, as `.csv` and as `.md`, from the
same rows.

## What writes what

| artifact | written by |
|---|---|
| `tableA_proteins`, `figA_methods_per_protein.png` | `tools/final_results.py`, block A |
| `tableB_lopo4`, `figB1_zeroshot_timelines.png`, `figB2_lopo_medians.png` | `tools/final_results.py`, block B |
| `tableC_negative_controls`, `figC_controls.png`, `tableC2_pool_choice` | `tools/final_results.py`, block C |
| `tableD_method_matrix`, `figD_method_matrix.png`, `figD_method_matrix_tau{1000,2000,4000}.png` | `tools/final_results.py`, block D |
| `tableE_ml_quality` | `tools/final_results.py`, block E |
| `figF_training_quality.png` | `tools/final_results.py`, block F |
| `cx43_arms.json` | `tools/cx43_arms.py`, which has to run first: `final_results.py` reads it and leaves the Cx43 group of the matrices empty when it is missing |
| `figG_shap_example.png` | `tools/shap_example_figure.py` |
| `figG2_shap_mthk.png` | the same script with the system, condition, replica, time window and output path passed to its `main()`; no command line in the repository reproduces it |
| `figure_2_1_workflow.png` | `tools/workflow_figure.py`, whose panel text is literal and does not follow the pipeline by itself |
| `figH_methods_showcase.png` | exported by hand from `notebooks/methods_showcase.ipynb` |
| `tableT5_one_moment.txt` | exported by hand from the same notebook |
