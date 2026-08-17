# Conducting-channels generalisation (tau = 4000 ps)

Training pool: gramicidin + MthK + KcsA-E71A (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.47 | x1.6800000000000002 |
| mthk | x1.0350000000000001 | x1.0550000000000002 |
| kcsa_family | x1.19 | x1.18 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| cx43:200mV/R2 | 1001 | 0.0749 | 0.000 | 0.346 | 1.43 |
| cx43:200mV/R4 | 1001 | 0.0000 | 0.000 | 0.347 | undefined |
| kcsa:G77AE71A/1 | 12501 | 0.0060 | 0.000 | 0.380 | 0.96 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.000 | 0.388 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.000 | 0.343 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0140 | 0.000 | 0.345 | 1.13 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0040 | 0.000 | 0.361 | 1.73 |

Columns dropped in physics-only: 21 (listed in results.json).