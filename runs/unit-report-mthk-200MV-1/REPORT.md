# Unit report — mthk:200MV/1 (cross-trained monitor, tau = 1000 ps)

Every prediction here comes from a model trained ONLY on: gramicidin:298K_2Msalt/01, gramicidin:298K_2Msalt/02, gramicidin:298K_2Msalt/03, gramicidin:298K_2Msalt/04, gramicidin:298K_2Msalt/05, gramicidin:330K_1Msalt/01, gramicidin:330K_1Msalt/02, gramicidin:330K_1Msalt/03, gramicidin:330K_1Msalt/04, gramicidin:330K_1Msalt/05, kcsa:E71A/1.

## Headline
- frames: 2001; base rate: 0.3245
- AP = 0.5051, x1.56 over this unit's own chance
- READY fraction at the training threshold: 0.1494
- monitor: 0/141 crossings warned; duty cycle 0.0000; random-null warned fraction nan

## Mechanism axes named on this unit's frames
- ion placement and desolvation: 1613 frames
- permeant delivery: 178 frames
- column immobility: 145 frames
- pathway/lining rearrangement: 64 frames
- constriction dehydration: 1 frames

## Reality checks
- training pool (leakage impossible by construction): family 'mthk' fully excluded
- event counts: full-membrane collector 175; label positive runs 141
- reproducibility vs stored experiment scores: {'max_prob_deviation': 0.0, 'ok': True}
- authors'-counter divergence explained by measurement: their state machine advances only on strict +-1 region steps over ~3 A half-filter regions, while Na+ in the broken TVAYG filter moves in bursts of up to 3.0 A per 40 ps frame (ion traces on file), so fast passes stall their counter; on the normal regime (E71A, K+) the two counters agree (51 vs 49)