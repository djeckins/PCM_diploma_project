# Unit report — gramicidin:298K_2Msalt/01 (cross-trained monitor, tau = 1000 ps)

Every prediction here comes from a model trained ONLY on: kcsa:E71A/1, mthk:200MV/1, mthk:200MV/2, mthk:neg-200MV/1, mthk:neg-200MV/2.

## Headline
- frames: 1251; base rate: 0.0611
- AP = 0.0908, x1.49 over this unit's own chance
- READY fraction at the training threshold: 0.3405
- monitor: 0/12 crossings warned; duty cycle 0.0671; random-null warned fraction 0.065

## Mechanism axes named on this unit's frames
- column immobility: 433 frames
- constriction dehydration: 348 frames
- ion placement and desolvation: 257 frames
- permeant delivery: 213 frames

## Reality checks
- training pool (leakage impossible by construction): family 'gramicidin' fully excluded
- event counts: full-membrane collector 13; label positive runs 12
- reproducibility vs stored experiment scores: {'max_prob_deviation': 0.0, 'ok': True}
- authors'-counter divergence explained by measurement: their state machine advances only on strict +-1 region steps over ~3 A half-filter regions, while Na+ in the broken TVAYG filter moves in bursts of up to 3.0 A per 40 ps frame (ion traces on file), so fast passes stall their counter; on the normal regime (E71A, K+) the two counters agree (51 vs 49)