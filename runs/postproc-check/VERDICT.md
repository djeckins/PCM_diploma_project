# Post-processing candidates: measured verdict

Per-frame AP gains from causal score smoothing (x3.91 -> x5.72 on cx43)
looked spectacular but FAILED the decisive event-level test at matched
time-in-warning (event_level.json): smoothing is neutral on cx43 (raw and
smoothed both warn 2/2 crossings from 5% duty vs 0.28 random null) and
HARMFUL on gramicidin (0.406 -> 0.348 warned at 5% duty) and KcsA
(0.388 -> 0.347). The per-frame gain is within-block credit multiplication
over autocorrelated labels — the artifact predicted by the literature
(point-adjust inflation, Kim et al. AAAI 2022; sample-vs-alarm gap in
seizure prediction; overlapping-window bias, Hammerla & Plotz 2015).

Decisions:
- Smoothing is NOT adopted for any ranking/skill claim. It already lives
  where it belongs: the event-level alarm's hysteresis smoothing
  (firing-power-style post-processing; Teixeira et al. 2012).
- Raw per-frame AP stays the canonical ranking metric everywhere.
- Positive finding kept for the thesis: the zero-shot cx43 monitor warns
  BOTH crossings at 5% time-in-warning vs 0.28 for a structure-preserving
  random alarm of the same duty; gramicidin cross-trained: 0.41 vs 0.21.
- The GBT+Linear ensemble (x2.59 -> x2.74 in-system gA; x2.84 only with the rejected smoothing on top) and max-MCC
  thresholds (MthK 0.06 -> 0.25) remain measured-but-unadopted options
  (secondary rows at most; results in results.json).
