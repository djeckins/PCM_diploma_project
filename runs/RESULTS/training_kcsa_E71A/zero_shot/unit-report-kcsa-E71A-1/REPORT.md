# Unit report — kcsa:E71A/1 (cross-trained monitor, tau = 1000 ps)

Every prediction here comes from a model trained ONLY on: gramicidin:298K_2Msalt/01, gramicidin:298K_2Msalt/02, gramicidin:298K_2Msalt/03, gramicidin:298K_2Msalt/04, gramicidin:298K_2Msalt/05, gramicidin:330K_1Msalt/01, gramicidin:330K_1Msalt/02, gramicidin:330K_1Msalt/03, gramicidin:330K_1Msalt/04, gramicidin:330K_1Msalt/05, mthk:200MV/1, mthk:200MV/2, mthk:neg-200MV/1, mthk:neg-200MV/2.

## Headline
- frames: 6251; base rate: 0.0488
- AP = 0.0596, x1.22 over this unit's own chance
- READY fraction at the training threshold: 0.0710
- monitor: 0/49 crossings warned; duty cycle 0.0117; random-null warned fraction 0.025

## Mechanism axes named on this unit's frames
- column immobility: 5348 frames
- permeant delivery: 587 frames
- constriction dehydration: 297 frames
- unclear: 19 frames

## Reality checks
- training pool (leakage impossible by construction): family 'kcsa_family' fully excluded
- event counts: full-membrane collector 49; label positive runs 49
- reproducibility vs stored experiment scores: {'max_prob_deviation': 0.0, 'ok': True}
- authors'-counter divergence explained by measurement: their state machine advances only on strict +-1 region steps over ~3 A half-filter regions, while Na+ in the broken TVAYG filter moves in bursts of up to 3.0 A per 40 ps frame (ion traces on file), so fast passes stall their counter; on the normal regime (E71A, K+) the two counters agree (51 vs 49)