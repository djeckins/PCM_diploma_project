# Method benchmark — cx43 (horizon tau = 4000 ps)

Each method is judged against what it claims to predict, on the frames where both its own inputs and the label are defined — which is why the frames column is not the same for every row. `verdict_accuracy` is the fraction of valid frames where the method's READY/NOT_READY reading (or its stated proxy) matched whether a crossing actually completed within tau; `ap_over_chance` is ranking quality relative to the base rate.

| method | claims | frames | verdict acc. | balanced acc. | AP / chance |
|---|---|---|---|---|---|
| clock | temporal null: event seriality only (time since last crossing, event rate, pore ion count); no structural input | 0 | nan | nan | nan |
| linear | this work (linear branch): same probabilistic claim, ridge-logistic on all descriptors | 0 | nan | nan | nan |
| linear_control_structural | fitted linear control on the published criterion's inputs (constriction radius + lining hydrophobicity) | 0 | nan | nan | nan |
| plsfma_coords | PLS-FMA (Krivobokova et al. 2012; the authors' g_fma PLS core, Helland/Denham): collective mode maximally covarying with the readiness label, regressed on superimposed C-alpha Cartesian coordinates — the published method's machinery applied to a binary functional quantity | 0 | nan | nan | nan |
| published_rao2019 | Rao et al. 2019 heuristic: local DEWETTING of the constriction (frame-wise application); readiness proxy = wetted | 1922 | 0.037 | 0.500 | 0.78 |
| published_rao2019_win | Rao et al. 2019 heuristic on window-averaged inputs (closer to the authors' time-averaged definition) | 1922 | 0.037 | 0.500 | 0.65 |
| published_rao2019_sigma | Rao et al. 2019 STRUCTURE classifier as published: sum of contour distances (sigma_d) over all pore-facing residue points past the 1RT line; non-conductive when sigma_d > 0.55 — the rule the paper's AUC 0.91 belongs to | 1922 | 0.037 | 0.500 | 1.00 |

This system additionally reports, per frame: the calibrated probability of readiness, the READY/NOT_READY verdict, and the diagnosed pore state (mechanism axis with the measured value it rests on) — columns `pcm_prob_ready`, `pcm_verdict`, `pcm_pore_state`, `pcm_reason` in `benchmark.parquet`/`.csv`.