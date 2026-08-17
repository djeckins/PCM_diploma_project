# Unit report — gramicidin:330K_1Msalt/02 (cross-trained monitor, tau = 1000 ps)

Every prediction here comes from a model trained ONLY on: kcsa:E71A/1, mthk:200MV/1, mthk:200MV/2, mthk:neg-200MV/1, mthk:neg-200MV/2.

## Headline
- frames: 1251; base rate: 0.0297
- AP = 0.0379, x1.27 over this unit's own chance
- READY fraction at the training threshold: 0.5068
- monitor: 1/6 crossings warned; duty cycle 0.1735; random-null warned fraction 0.200

## Mechanism axes named on this unit's frames
- ion placement and desolvation: 625 frames
- permeant delivery: 234 frames
- column immobility: 228 frames
- constriction dehydration: 164 frames

## Reality checks
- training pool (leakage impossible by construction): family 'gramicidin' fully excluded
- event counts: full-membrane collector 6; label positive runs 6
- reproducibility vs stored experiment scores: {'max_prob_deviation': 0.0, 'ok': True}
- authors'-counter divergence explained by measurement: their state machine advances only on strict +-1 region steps over ~3 A half-filter regions, while Na+ in the broken TVAYG filter moves in bursts of up to 3.0 A per 40 ps frame (ion traces on file), so fast passes stall their counter; on the normal regime (E71A, K+) the two counters agree (51 vs 49)