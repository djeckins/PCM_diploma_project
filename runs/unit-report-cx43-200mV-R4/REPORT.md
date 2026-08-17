# Unit report — cx43:200mV/R4 (cross-trained monitor, tau = 4000 ps)

Every prediction here comes from a model trained ONLY on: gramicidin:298K_2Msalt/01, gramicidin:298K_2Msalt/02, gramicidin:298K_2Msalt/03, gramicidin:298K_2Msalt/04, gramicidin:298K_2Msalt/05, gramicidin:330K_1Msalt/01, gramicidin:330K_1Msalt/02, gramicidin:330K_1Msalt/03, gramicidin:330K_1Msalt/04, gramicidin:330K_1Msalt/05, kcsa:E71A/1, mthk:200MV/1, mthk:200MV/2, mthk:neg-200MV/1, mthk:neg-200MV/2.

## Headline
- frames: 1001; base rate: 0.0000
- AP undefined (no positive frames)
- READY fraction at the training threshold: 0.0789
- monitor: 0/0 crossings warned; duty cycle 0.0000; random-null warned fraction nan

## Mechanism axes named on this unit's frames
- unclear: 784 frames
- permeant delivery: 180 frames
- constriction dehydration: 37 frames

## Reality checks
- training pool (leakage impossible by construction): family 'cx43' fully excluded
- event counts: full-membrane collector 0; label positive runs 0
- reproducibility vs stored experiment scores: experiment scores not found
- authors'-counter divergence explained by measurement: their state machine advances only on strict +-1 region steps over ~3 A half-filter regions, while Na+ in the broken TVAYG filter moves in bursts of up to 3.0 A per 40 ps frame (ion traces on file), so fast passes stall their counter; on the normal regime (E71A, K+) the two counters agree (51 vs 49)