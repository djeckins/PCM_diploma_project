# Conducting-channels generalisation (tau = 2000 ps)

Training pool: gramicidin + MthK + KcsA-E71A + Cx43 (conducting channels only).

## Leave-one-protein-out (median AP over the held protein's own chance)

| held-out protein | full columns | physics-only |
|---|---|---|
| gramicidin | x1.53 | x1.16 |
| mthk | x1.155 | x1.165 |
| kcsa_family | x1.2 | x1.17 |
| cx43 | x2.08 | x2.7 |

## Negative control: monitor applied to unseen mutant arms

| unit | frames | base rate | READY fraction | mean P(ready) | AP/chance |
|---|---|---|---|---|---|
| kcsa:G77AE71A/1 | 12501 | 0.0030 | 0.000 | 0.027 | 2.02 |
| kcsa:T75AE71A/1 | 6251 | 0.0000 | 0.000 | 0.020 | undefined |
| kcsa_na:E71A/1 | 6251 | 0.0000 | 0.000 | 0.041 | undefined |
| kcsa_na:G77AE71A/1 | 12501 | 0.0069 | 0.000 | 0.037 | 5.19 |
| kcsa_na:T75AE71A/1 | 6251 | 0.0021 | 0.000 | 0.034 | 0.92 |

Columns dropped in physics-only: 38 (listed in results.json).