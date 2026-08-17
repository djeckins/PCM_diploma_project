# Unit report — gramicidin:330K_1Msalt/04 (cross-trained monitor, tau = 1000 ps)

Every prediction here comes from a model trained ONLY on: kcsa:E71A/1, mthk:200MV/1, mthk:200MV/2, mthk:neg-200MV/1, mthk:neg-200MV/2.

## Headline
- frames: 1251; base rate: 0.0490
- AP = 0.1006, x2.05 over this unit's own chance
- READY fraction at the training threshold: 0.5172
- monitor: 1/10 crossings warned; duty cycle 0.0440; random-null warned fraction 0.084

## Mechanism axes named on this unit's frames
- ion placement and desolvation: 743 frames
- permeant delivery: 289 frames
- column immobility: 132 frames
- constriction dehydration: 87 frames

## Reality checks
- training pool (leakage impossible by construction): family 'gramicidin' fully excluded
- event counts: full-membrane collector 10; label positive runs 10
- reproducibility vs stored experiment scores: {'max_prob_deviation': 0.0, 'ok': True}
- authors'-counter divergence explained by measurement: their state machine advances only on strict +-1 region steps over ~3 A half-filter regions, while Na+ in the broken TVAYG filter moves in bursts of up to 3.0 A per 40 ps frame (ion traces on file), so fast passes stall their counter; on the normal regime (E71A, K+) the two counters agree (51 vs 49)