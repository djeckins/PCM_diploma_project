# Unit report — mthk:200MV/2 (cross-trained monitor, tau = 1000 ps)

Every prediction here comes from a model trained ONLY on: gramicidin:298K_2Msalt/01, gramicidin:298K_2Msalt/02, gramicidin:298K_2Msalt/03, gramicidin:298K_2Msalt/04, gramicidin:298K_2Msalt/05, gramicidin:330K_1Msalt/01, gramicidin:330K_1Msalt/02, gramicidin:330K_1Msalt/03, gramicidin:330K_1Msalt/04, gramicidin:330K_1Msalt/05, kcsa:E71A/1.

## Headline
- frames: 2001; base rate: 0.1983
- AP = 0.3950, x1.99 over this unit's own chance
- READY fraction at the training threshold: 0.1239
- monitor: 0/89 crossings warned; duty cycle 0.0000; random-null warned fraction nan

## Mechanism axes named on this unit's frames
- ion placement and desolvation: 1668 frames
- column immobility: 154 frames
- permeant delivery: 121 frames
- pathway/lining rearrangement: 54 frames
- constriction dehydration: 4 frames

## Reality checks
- training pool (leakage impossible by construction): family 'mthk' fully excluded
- event counts: full-membrane collector 104; label positive runs 89
- reproducibility vs stored experiment scores: {'max_prob_deviation': 0.0, 'ok': True}
- authors'-counter divergence explained by measurement: their state machine advances only on strict +-1 region steps over ~3 A half-filter regions, while Na+ in the broken TVAYG filter moves in bursts of up to 3.0 A per 40 ps frame (ion traces on file), so fast passes stall their counter; on the normal regime (E71A, K+) the two counters agree (51 vs 49)