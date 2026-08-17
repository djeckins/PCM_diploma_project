# Conducting-channels generalisation (tau = 2000 ps)

Training pool: gramicidin + MthK + KcsA-E71A (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.605 | x1.7349999999999999 |
| mthk | x1.1400000000000001 | x1.13 |
| kcsa_family | x1.21 | x1.31 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| cx43:200mV/R2 | 1001 | 0.0408 | 0.000 | 0.100 | 2.01 |
| cx43:200mV/R4 | 1001 | 0.0000 | 0.000 | 0.100 | undefined |
| kcsa:G77AE71A/1 | 12501 | 0.0030 | 0.000 | 0.027 | 2.04 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.000 | 0.031 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.000 | 0.050 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0069 | 0.000 | 0.049 | 4.41 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0021 | 0.000 | 0.076 | 0.65 |

Columns dropped in physics-only: 21 (listed in results.json).