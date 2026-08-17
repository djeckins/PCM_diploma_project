# Method benchmark — kcsa (horizon tau = 1000 ps)

Each method is judged against what it claims to predict, on the frames where both its own inputs and the label are defined — which is why the frames column is not the same for every row. `verdict_accuracy` is the fraction of valid frames where the method's READY/NOT_READY reading (or its stated proxy) matched whether a crossing actually completed within tau; `ap_over_chance` is ranking quality relative to the base rate.

| method | claims | frames | verdict acc. | balanced acc. | AP / chance |
|---|---|---|---|---|---|
| trees | this work (head model): P(crossing completes within tau), READY/NOT_READY verdict and mechanism axis (pore state) per frame | 18738 | 0.998 | 0.500 | 3.14 |
| clock | temporal null: event seriality only (time since last crossing, event rate, pore ion count); no structural input | 6244 | 1.000 | nan | 0.55 |
| linear | this work (linear branch): same probabilistic claim, ridge-logistic on all descriptors | 6244 | 0.999 | nan | 0.57 |
| linear_control_structural | fitted linear control on the published criterion's inputs (constriction radius + lining hydrophobicity) | 6244 | 1.000 | nan | 0.55 |
| plsfma_coords | PLS-FMA (Krivobokova et al. 2012; the authors' g_fma PLS core, Helland/Denham): collective mode maximally covarying with the readiness label, regressed on superimposed C-alpha Cartesian coordinates — the published method's machinery applied to a binary functional quantity | 6244 | 0.000 | nan | 0.79 |
| published_rao2019 | Rao et al. 2019 heuristic: local DEWETTING of the constriction (frame-wise application); readiness proxy = wetted | 24982 | 0.861 | 0.550 | 1.84 |
| published_rao2019_win | Rao et al. 2019 heuristic on window-averaged inputs (closer to the authors' time-averaged definition) | 24982 | 0.980 | 0.500 | 1.43 |
| published_rao2019_sigma | Rao et al. 2019 STRUCTURE classifier as published: sum of contour distances (sigma_d) over all pore-facing residue points past the 1RT line; non-conductive when sigma_d > 0.55 — the rule the paper's AUC 0.91 belongs to | 24982 | 0.986 | 0.499 | 0.63 |

This system additionally reports, per frame: the calibrated probability of readiness, the READY/NOT_READY verdict, and the diagnosed pore state (mechanism axis with the measured value it rests on) — columns `pcm_prob_ready`, `pcm_verdict`, `pcm_pore_state`, `pcm_reason` in `benchmark.parquet`/`.csv`.

Event-level monitor (this work only): 1.00 of crossings warned in advance vs 1.00 for a structure-preserving random alarm.